import pandas as pd

# tableau20 from # http://www.randalolson.com/2014/06/28/how-to-make-beautiful-data-visualizations-in-python-with-matplotlib/
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150), 
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib
# accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)   

def liftover_bed(bed, chain, liftOver=None, mapped=None, unmapped=None):
    """
    Lift over a bed file using a given chain file. 

    Parameters
    ----------
    bed : str or pybedtools.BedTool
        Coordinates to lift over.
        
    chain : str
        Path to chain file to use for lift over.

    liftOver : str
        Path to liftOver executable if not in path.

    mapped : str
        Path for bed file with coordinates that are lifted over correctly.

    unmapped : str
        Path for text file to store coordinates that did not lift over
        correctly. If this is not provided, these are discarded.

    Returns
    -------
    new_coords : pandas.DataFrame
        Pandas data frame with lift over results. Index is old coordinates in
        the form chrom:start-end and columns are chrom, start, end and loc
        (chrom:start-end) in new coordinate system.
    """
    import pybedtools as pbt
    if mapped == None:
        import tempfile
        mapped = tempfile.NamedTemporaryFile()
        mname = mapped.name
    else:
        mname = mapped
    if unmapped == None:
        import tempfile
        unmapped = tempfile.NamedTemporaryFile()
        uname = unmapped.name
    else:
        uname = unmapped
    if liftOver == None:
        liftOver = 'liftOver'
    if type(bed) == str:
        bt = pbt.BedTool(bed)
    bt = bt.sort()
    c = '{} {} {} {} {}'.format(liftOver, bt.fn, chain, mname, uname)
    subprocess.check_call(c, shell=True)
    with open(uname) as f:
        missing = pbt.BedTool(''.join([x for x in f.readlines()[1::2]]),
                              from_string=True)
    bt = bt.subtract(missing)
    bt_mapped = pbt.BedTool(mname)
    old_loc = []
    for r in bt:
        old_loc.append('{}:{}-{}'.format(r.chrom, r.start, r.end))
    new_loc = []
    new_chrom = []
    new_start = []
    new_end = []
    for r in bt_mapped:
        new_loc.append('{}:{}-{}'.format(r.chrom, r.start, r.end))
        new_chrom.append(r.chrom)
        new_start.append(r.start)
        new_end.append(r.end)
    new_coords = pd.DataFrame({'loc':new_loc, 'chrom': new_chrom, 
                               'start': new_start, 'end': new_end},
                              index=old_loc)
    for f in [mapped, unmapped]:
        try:
            f.close()
        except AttributeError:
            continue
    return new_coords

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
    import os
    import readline
    import statsmodels.stats.multitest as smm
    import rpy2.robjects as r
    genes = list(genes)
    sig = [bool(x) for x in sig]
    r.r('suppressMessages(library(goseq))')
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
    if plot_fn and os.path.exists('Rplots.pdf'):
        from os import rename
        rename('Rplots.pdf', plot_fn)
    elif os.path.exists('Rplots.pdf'):
        from os import remove
        remove('Rplots.pdf')
    return go_results

def make_color_legend_rects(colors, labels=None):
    """ 
    Make list of rectangles and labels for making legends.

    Parameters
    ----------
    colors : pandas.Series or list
        Pandas series whose values are colors and index is labels.
        Alternatively, you can provide a list with colors and provide the labels
        as a list.

    labels : list
        If colors is a list, this should be the list of corresponding labels.

    Returns
    -------
    out : pd.Series
        Pandas series whose values are matplotlib rectangles and whose index are
        the legend labels for those rectangles.

    """
    from matplotlib.pyplot import Rectangle
    if labels:
        d = dict(zip(labels, colors))
        se = pd.Series(d)
    else:
        se = colors
    rects = []
    for i in se.index:
        r = Rectangle((0, 0), 0, 0, fc=se[i])
        rects.append(r)
    out = pd.Series(rects, index=se.index)
    return out 

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

    def plot_variance_explained(self, cumulative=False, xtick_start=1,
                                xtick_spacing=1, num_pc=None):
        """ 
        Plot amount of variance explained by each principal component.
    
        Parameters
        ----------
        num_pc : int
            Number of principal components to plot. If None, plot all.
    
        cumulative : bool
            If True, include cumulative variance.
    
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
        if cumulative:
            s_cumsum = s_norm.cumsum()
            plt.bar(range(s_cumsum.shape[0]), s_cumsum.values,
                    label='Cumulative', color=(0.17254901960784313,
                                               0.6274509803921569,
                                               0.17254901960784313))
            plt.bar(range(s_norm.shape[0]), s_norm.values, label='Per PC',
                    color=(0.12156862745098039, 0.4666666666666667,
                           0.7058823529411765))
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.ylabel('Variance')
        else:
            plt.bar(range(s_norm.shape[0]), s_norm.values,
                    color=(0.12156862745098039, 0.4666666666666667,
                           0.7058823529411765))
            plt.ylabel('Proportion variance explained')
        plt.xlabel('PC')
        plt.xlim(0, s_norm.shape[0])
        tick_locs = arange(xtick_start - 1, s_norm.shape[0],
                              step=xtick_spacing)
        # 0.8 is the width of the bars.
        tick_locs = tick_locs + 0.4 
        plt.xticks(tick_locs, 
                   arange(xtick_start, s_norm.shape[0] + 1, xtick_spacing))

    def plot_pc_scatter(self, pc1, pc2, v=True, subset=None, ax=None,
                        color=None, s=None, marker=None, color_name=None,
                        s_name=None, marker_name=None):
        """
        Make a scatter plot of two principal components. You can create
        differently colored, sized, or marked scatter points.

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

        color : pandas.Series
            Pandas series containing a categorical variable to color the scatter
            points. Currently limited to 10 distinct values (colors).

        s : pandas.Series
            Pandas series containing a categorical variable to size the scatter
            points. Currently limited to 7 distinct values (sizes).

        marker : pandas.Series
            Pandas series containing a categorical variable to choose the marker
            type for the scatter points. Currently limited to 21 distinct values
            (marker styles).

        color_name : str
            Name for the color legend if a categorical variable for color is
            provided.

        s_name : str
            Name for the size legend if a categorical variable for size is
            provided.

        marker_name : str
            Name for the marker legend if a categorical variable for marker type
            is provided.

        Returns
        -------
        ax : matplotlib.axes._subplots.AxesSubplot
            Scatter plot axis.

        TODO: Add ability to label points. 
        """
        import matplotlib.pyplot as plt
        if v:
            df = self.v
        else:
            df = self.u
        if color is not None:
            colormap = pd.Series(dict(zip(set(color.values),
                                          tableau20[0:len(set(color)) + 1:2])))
            color = pd.Series([colormap[x] for x in color.values],
                              index=color.index)
            color_legend = True
            if not color_name:
                color_name = color.index.name
        else:
            color = pd.Series([tableau20[0]] * df.shape[0], index=df.index)
            color_legend = False
        if s is not None:
            smap = pd.Series(dict(zip(
                set(s.values), range(30, 351)[0::50][0:len(set(s)) + 1])))
            s = pd.Series([smap[x] for x in s.values],
                          index=s.index)
            s_legend = True
            if not s_name:
                s_name = s.index.name
        else:
            s = pd.Series(30, index=df.index)
            s_legend = False
        markers = ['o', '*', 's', 'v', '+', 'x', 'd', 
                   'p', '2', '<', '|', '>', '_', 'h', 
                   '1', '2', '3', '4', '8', '^', 'D']
        if marker is not None:
            markermap = pd.Series(dict(zip(set(marker.values), markers)))
            marker = pd.Series([markermap[x] for x in marker.values],
                               index=marker.index)
            marker_legend = True
            if not marker_name:
                marker_name = marker.index.name
        else:
            marker = pd.Series('o', index=df.index)
            marker_legend = False
        if ax is None:
            fig, ax = plt.subplots(1, 1)
        for m in set(marker.values):
            mse = marker[marker == m]
            cse = color[mse.index]
            sse = s[mse.index]
            ax.scatter(df.ix[mse.index, pc1], df.ix[mse.index, pc2],
                       s=sse.values, color=list(cse.values), marker=m, 
                       alpha=0.8)
        
        ax.set_title('{} vs. {}'.format(pc1, pc2))
        ax.set_xlabel(pc1)
        ax.set_ylabel(pc2)
    
        if color_legend:
            legend_rects = make_color_legend_rects(colormap)
            for r in legend_rects:
                ax.add_patch(r)
            lgd = ax.legend(legend_rects.values, labels=legend_rects.index, 
                             title=color_name,
                             loc='upper left',
                             bbox_to_anchor=(1, 1))
        
        if s_legend:
            if lgd:
                lgd = ax.add_artist(lgd)
            xa, xb = ax.get_xlim()
            ya, yb = ax.get_ylim()
            for i in smap.index:
                ax.scatter([xb + 1], [yb + 1], marker='o', 
                           s=smap[i], color='black', label=i)
            lgd = ax.legend(title=s_name, loc='center left', 
                            bbox_to_anchor=(1, 0.5))
            ax.set_xlim(xa, xb)
            ax.set_ylim(ya, yb)
            
        if marker_legend:
            if lgd:
                lgd = ax.add_artist(lgd)
            xa, xb = ax.get_xlim()
            ya, yb = ax.get_ylim()
            for i in markermap.index:
                t = ax.scatter([xb + 1], [yb + 1], marker=markermap[i], 
                               s=smap.min(), color='black', label=i)
                
            handles, labels = ax.get_legend_handles_labels()
            if s_legend:
                handles = handles[len(smap):]
                labels = labels[len(smap):]
            lgd = ax.legend(handles, labels, title=marker_name, 
                            loc='lower left', bbox_to_anchor=(1, 0))
            ax.set_xlim(xa, xb)
            ax.set_ylim(ya, yb)
        fig.tight_layout()
        return fig, ax
    
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
