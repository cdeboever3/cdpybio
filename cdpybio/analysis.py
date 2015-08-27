import pandas as pd

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

    def plot_variance_explained(self):
        import matplotlib.pyplot as plt
        plt.bar(range(self.s_norm.shape[0]), self.s_norm.values)
        plt.xlabel('PC')
        plt.ylabel('Proportion variance explained')

    def plot_pc_scatter(self, pc1, pc2, v=True, subset=None):
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

        TODO: Add ability to have facets that control point size, color, and
        shape. Add legends for these. Add ability to label points. Return
        axis object or something useful.
        """
        import matplotlib.pyplot as plt
        from numpy import arange
        if v:
            df = self.v
        else:
            df = self.u
        plt.scatter(df[pc1], df[pc2])
        plt.title('{} vs. {}'.format(pc1, pc2))
        plt.xlabel(pc1)
        plt.ylabel(pc2)
        plt.xlim(0, self.s_norm.shape[0])
        tick_locs = arange(xtick_start - 1, self.s_norm.shape[0],
                              step=xtick_spacing)
        # 0.8 is the width of the bars.
        tick_locs += 0.4 
        plt.xticks(tick_locs, 
                   arange(xtick_start, self.s_norm.shape[0], xtick_spacing))
