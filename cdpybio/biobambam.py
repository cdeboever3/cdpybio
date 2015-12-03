import numpy as np
import pandas as pd

def parse_bammarkduplicates(fn):
    """
    Parse the output from biobambam2's bammarkduplicates and return as pandas
    Series.

    Parameters
    ----------
    fn : str
        Path to the output file to parse.

    Returns
    -------
    metrics : pandas.Series
        Duplicate metrics.

    hist : pandas.Series
        Duplicate histogram.

    """
    with open(fn) as f:
        lines = [x.strip().split('\t') for x in f.readlines()]
    metrics = pd.Series(lines[4], lines[3])
    m = pd.to_numeric(metrics[metrics.index[1:]])
    metrics[m.index] = m.values

    vals = np.array(lines[8:-1])
    hist = pd.Series(vals[:, 1], index=[int(float(x)) for x in vals[:, 0]])
    hist = pd.to_numeric(hist)
    return metrics, hist
