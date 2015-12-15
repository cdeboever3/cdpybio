import pandas as pd

def combine_counts(
    fns,
    define_sample_name=None,
):
    """
    Combine featureCounts output files for multiple samples.

    Parameters
    ----------
    fns : list of strings 
        Filenames of featureCounts output files to combine.

    define_sample_name : function
        A function mapping the featureCounts output filenames to sample names.
        If this is not provided, the header of the last column in the
        featureCounts output will be used as the sample name.

    Returns
    -------
    combined_counts :  pandas.DataFrame
        Combined featureCount counts.
    
    """
    counts = []
    for fn in fns:
        df = pd.read_table(fn, skiprows=1, index_col=0)
        counts.append(df[df.columns[-1]])
    combined_counts = pd.DataFrame(counts).T
    if define_sample_name:
        names = [define_sample_name(x) for x in fns]
        combined_counts.columns = names
    combined_counts.index.name = ''
    return combined_counts
