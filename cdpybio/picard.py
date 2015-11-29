import numpy as np
import pandas as pd

def parse_bam_index_stats(fn):
    """
    Parse the output from Picard's BamIndexStast and return as pandas Dataframe.

    Parameters
    ----------
    filename : str of filename or file handle
        Filename of the Picard output you want to parse.

    Returns
    -------
    df : pandas.DataFrame
        Data from output file.

    """
    with open(fn) as f:
        lines = [x.strip().split() for x in f.readlines()]
    no_counts = int(lines[-1][-1])
    lines = lines[:-1]
    chrom = [x[0] for x in lines]
    length = [int(x[2]) for x in lines]
    aligned = [int(x[4]) for x in lines]
    unaligned = [int(x[6]) for x in lines]
    df = pd.DataFrame([length, aligned, unaligned], columns=chrom, 
                      index=['length', 'aligned', 'unaligned']).T
    df = df.ix[sorted(df.index)]
    return df

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
    m = pd.to_numeric(metrics[metrics.index[1:]])
    metrics[m.index] = m.values

    vals = np.array(lines[11:-1])
    hist = pd.Series(vals[:, 1], index=[int(float(x)) for x in vals[:, 0]])
    hist = pd.to_numeric(hist)
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
    hist = pd.to_numeric(hist)
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
    hist = pd.to_numeric(hist)
    return metrics, hist
