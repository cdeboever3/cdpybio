import numpy as np
import pandas as pd

def parse_mark_duplicate_metrics(fn):
    """
    Parse the output from Picard's MarkDuplicates and return as pandas
    Series.

    Parameters
    ----------
    filename : str of filename or file handle
        Filename of the Picard output you want to parse.

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
    metrics = pd.to_numeric(metrics)

    vals = np.array(lines[8:-1])
    hist = pd.Series(vals[:, 1], index=[int(float(x)) for x in vals[:, 0]])
    hist = pd.to_numeric(hist)
    return metrics, hist
