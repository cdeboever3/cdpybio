import pandas as pd

def goseq_gene_enrichment(genes, sig, plot_fn=None, length_correct=True):
    """
    Perform goseq enrichment for an Ensembl gene set.

    Parameters
    ----------
    genes : list
        List of all genes as Ensembl IDs.
        
    sig : list
        List of boolean values indicating whether each gene is significant or
        not.

    plot_fn : str
        Path to save length bias plot to. If not provided, the plot is deleted.

    length_correct : bool
        Correct for length bias.

    Returns
    -------
    go_results : pandas.DataFrame
        Dataframe with goseq results as well as Benjamini-Hochberg correct
        p-values.
    """
    import readline
    import statsmodels.stats.multitest as smm
    import rpy2.robjects as r
    genes = list(genes)
    sig = [bool(x) for x in sig]
    r.r('library(goseq)')
    r.globalenv['genes'] = list(genes)
    r.globalenv['group'] = list(sig)
    r.r('group = as.logical(group)')
    r.r('names(group) = genes')
    r.r('pwf = nullp(group, "hg19", "ensGene")')
    if length_correct:
        r.r('wall = goseq(pwf, "hg19", "ensGene")')
    else:
        r.r('wall = goseq(pwf, "hg19", "ensGene", method="Hypergeometric")')
    r.r('t = as.data.frame(wall)')
    t = r.globalenv['t']
    go_results = pd.DataFrame(columns=list(t.colnames))
    for i, c in enumerate(go_results.columns):
        go_results[c] = list(t[i])
    r, c, ask, abf = smm.multipletests(
        go_results.over_represented_pvalue, alpha=0.05, method='fdr_i')
    go_results['over_represented_pvalue_bh'] = c 
    r, c, ask, abf = smm.multipletests(
        go_results.under_represented_pvalue, alpha=0.05, method='fdr_i')
    go_results['under_represented_pvalue_bh'] = c
    go_results.index = go_results.category
    go_results = go_results.drop('category', axis=1)
    if plot_fn:
        from os import rename
        rename('Rplots.pdf', plot_fn)
    else:
        from os import remove
        remove('Rplots.pdf')
    return go_results

class SVD:
    def __init__(self, df, mean_center=True, scale_variance=False):
        """
        Perform SVD for data matrix using scipy. Note that this is currently
        inefficient for large matrices due to some of the pandas operations.

        Parameters
        ----------
        df : pandas.DataFrame
            Pandas data frame with data.

        mean_center : bool
            If True, mean center the rows. This should be done if not already
            done.

        scale_variance : boole
            If True, scale the variance of each row to be one. Combined with
            mean centering, this will transform your data into z-scores.
        
        """
        import copy
        self.data_orig = copy.deepcopy(df)
        self.data = copy.deepcopy(df)
        if mean_center:
            self.data = (self.data.T - self.data.mean(axis=1)).T
        if scale_variance:
            self.data = (self.data.T / self.data.std(axis=1)).T
        self._perform_svd()

    def _perform_svd(self):
        from scipy.linalg import svd
        u, s, vh = svd(self.data, full_matrices=False)
        self.u_orig = u
        self.s_orig = s
        self.vh_orig = vh

        cols = ['PC{}'.format(x) for x in range(1, self.data.shape[1] + 1)]
        self.u = pd.DataFrame(u, index=self.data.index, columns=cols)
        self.v = pd.DataFrame(vh.T, index=self.data.columns, columns=cols)
        index = ['PC{}'.format(x) for x in range(1, len(s) + 1)]
        self.s_norm = pd.Series(s / s.sum(), index=index)

    def plot_variance_explained(self, xtick_start=1, xtick_spacing=1,
                                num_pc=None):
        """
        Plot amount of variance explained by each principal component.

        Parameters
        ----------
        num_pc : int
            Number of principal components to plot. If None, plot all.

        xtick_start : int
            The first principal component to label on the x-axis.

        xtick_spacing : int
            The spacing between labels on the x-axis.

        """
        import matplotlib.pyplot as plt
        from numpy import arange
        if num_pc:
            s_norm = self.s_norm[0:num_pc]
        else:
            s_norm = self.s_norm
        plt.bar(range(s_norm.shape[0]), s_norm.values)
        plt.xlabel('PC')
        plt.ylabel('Proportion variance explained')
        plt.xlim(0, s_norm.shape[0])
        tick_locs = arange(xtick_start - 1, s_norm.shape[0],
                              step=xtick_spacing)
        # 0.8 is the width of the bars.
        tick_locs = tick_locs + 0.4 
        plt.xticks(tick_locs, 
                   arange(xtick_start, s_norm.shape[0] + 1, xtick_spacing))

    def plot_pc_scatter(self, pc1, pc2, v=True, subset=None, ax=None):
        """
        Make a scatter plot of two principal components. 

        Parameters
        ----------
        pc1 : str
            String of form PCX where X is the number of the principal component
            you want to plot on the x-axis.
        
        pc2 : str
            String of form PCX where X is the number of the principal component
            you want to plot on the y-axis.

        v : bool
            If True, use the v matrix for plotting the principal components
            (typical if input data was genes as rows and samples as columns).
            If False, use the u matrix.

        subset : list
            Make the scatter plot using only a subset of the rows of u or v.

        ax : matplotlib.axes
            Plot the scatter plot on this axis.

        Returns
        -------
        ax : matplotlib.axes._subplots.AxesSubplot
            Scatter plot axis.

        TODO: Add ability to have facets that control point size, color, and
        shape. Add legends for these. Add ability to label points. 
        """
        import matplotlib.pyplot as plt
        if v:
            df = self.v
        else:
            df = self.u
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        ax.scatter(df[pc1], df[pc2])
        ax.set_title('{} vs. {}'.format(pc1, pc2))
        ax.set_xlabel(pc1)
        ax.set_ylabel(pc2)
        return ax
    
    def pc_correlation(self, covariates, num_pc=5):
        """
        Calculate the correlation between the first num_pc prinicipal components
        and known covariates. The size and index of covariates determines
        whether u or v is used.

        Parameters
        ----------
        covariates : pandas.DataFrame
            Dataframe of covariates whose index corresponds to the index of
            either u or v. 
        
        num_pc : int
            Number of principal components to correlate with.

        Returns
        -------
        corr : pandas.Panel
            Panel with correlation values and p-values.

        """
        from scipy.stats import spearmanr
        if (covariates.shape[0] == self.u.shape[0] and 
            len(set(covariates.index) & set(self.u.index)) == self.u.shape[0]):
            mat = self.u
        elif (covariates.shape[0] == self.v.shape[0] and 
            len(set(covariates.index) & set(self.v.index)) == self.v.shape[0]):
            mat = self.v
        corr = pd.Panel(items=['rho', 'pvalue'],
                        major_axis=covariates.columns,
                        minor_axis=mat.columns[0:num_pc])
        for i in corr.major_axis:
            for j in corr.minor_axis:
                rho, p = spearmanr(covariates[i], mat[j])
                corr.ix['rho', i, j] = rho
                corr.ix['pvalue', i, j] = p
        return corr
    
    def pc_anova(self, covariates, num_pc=5):
        """ 
        Calculate one-way ANOVA between the first num_pc prinicipal components
        and known covariates. The size and index of covariates determines
        whether u or v is used.
    
        Parameters
        ----------
        covariates : pandas.DataFrame
            Dataframe of covariates whose index corresponds to the index of
            either u or v. 
    
        num_pc : int
            Number of principal components to correlate with.
    
        Returns
        -------
        anova : pandas.Panel
            Panel with F-values and p-values.
    
        """
        from scipy.stats import f_oneway
        if (covariates.shape[0] == self.u.shape[0] and 
            len(set(covariates.index) & set(self.u.index)) == self.u.shape[0]):
            mat = self.u
        elif (covariates.shape[0] == self.v.shape[0] and 
            len(set(covariates.index) & set(self.v.index)) == self.v.shape[0]):
            mat = self.v
        anova = pd.Panel(items=['fvalue', 'pvalue'],
                         major_axis=covariates.columns,
                         minor_axis=mat.columns[0:num_pc])
        for i in anova.major_axis:
            for j in anova.minor_axis:
                t = [mat[j][covariates[i] == x] for x in set(covariates[i])]
                f, p = f_oneway(*t)
                anova.ix['fvalue', i, j] = f 
                anova.ix['pvalue', i, j] = p 
        return anova
