import numpy as np
import pandas as pd

def parse_alignment_summary_metrics(fn):
    """
    Parse the output from Picard's CollectAlignmentSummaryMetrics and return as
    pandas Dataframe.

    Parameters
    ----------
    filename : str of filename or file handle
        Filename of the Picard output you want to parse.

    Returns
    -------
    df : pandas.DataFrame
        Data from output file.

    """
    df = pd.read_table(fn, index_col=0, skiprows=range(6) + [10, 11]).T
    return df

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
    metrics = pd.Series(lines[7], lines[6])
    metrics = metrics.convert_objects(convert_numeric=True)

    vals = np.array(lines[11:-1])
    hist = pd.Series(vals[:, 1], index=[int(float(x)) for x in vals[:, 0]])
    hist = hist.convert_objects(convert_numeric=True)
    return metrics, hist

def parse_insert_metrics(fn):
    """
    Parse the output from Picard's CollectInsertSizeMetrics and return as pandas
    Series.

    Parameters
    ----------
    filename : str of filename or file handle
        Filename of the Picard output you want to parse.

    Returns
    -------
    metrics : pandas.Series
        Insert size metrics.

    hist : pandas.Series
        Insert size histogram.

    """
    with open(fn) as f:
        lines = [x.strip().split('\t') for x in f.readlines()]

    index = lines[6]
    vals = lines[7]
    for i in range(len(index) - len(vals)):
        vals.append(np.nan)
    for i, v in enumerate(vals):
        if type(v) == str:
            try:
                vals[i] = int(v)
            except ValueError:
                try: 
                    vals[i] = float(v)
                except ValueError:
                    continue
    metrics = pd.Series(vals, index=index)
    
    vals = np.array(lines[11:-1])
    hist = pd.Series(vals[:, 1], index=[int(float(x)) for x in vals[:, 0]])
    hist = hist.convert_objects(convert_numeric=True)
    return metrics, hist

def parse_rna_seq_metrics(fn):
    """
    Parse the output from Picard's CollectRnaSeqMetrics and return as pandas
    Series.

    Parameters
    ----------
    filename : str of filename or file handle
        Filename of the Picard output you want to parse.

    Returns
    -------
    metrics : pandas.Series
        Insert size metrics.

    hist : pandas.Series
        5'/3' bias histogram.

    """
    metrics = pd.read_table(fn, skiprows=6, index_col=0)
    with open(fn) as f:
        lines = [x.strip().split('\t') for x in f.readlines()]

    index = lines[6]
    vals = lines[7]
    for i in range(len(index) - len(vals)):
        vals.append(np.nan)
    for i, v in enumerate(vals):
        if type(v) == str:
            try:
                vals[i] = int(v)
            except ValueError:
                try:
                    vals[i] = float(v)
                except ValueError:
                    continue
    metrics = pd.Series(vals, index=index)

    vals = np.array(lines[11:-1])
    hist = pd.Series(vals[:, 1], index=[int(float(x)) for x in vals[:, 0]])
    hist = hist.convert_objects(convert_numeric=True)
    return metrics, hist
