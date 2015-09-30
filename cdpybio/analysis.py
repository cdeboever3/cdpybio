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
        if v:
            df = self.v
        else:
            df = self.u
        plt.scatter(df[pc1], df[pc2])
        plt.title('{} vs. {}'.format(pc1, pc2))
        plt.xlabel(pc1)
        plt.ylabel(pc2)
