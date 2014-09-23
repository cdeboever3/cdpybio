import pandas as pd

def strip_chr(bt):
    """Strip 'chr' from chromosomes for BedTool object

    Parameters
    ----------
    bt : pybedtools.BedTool
        BedTool to strip 'chr' from.

    Returns
    -------
    out : pybedtools.BedTool
        New BedTool with 'chr' stripped from chromosome names.

    """
    try:
        df = pd.read_table(exons.fn, header=None, dtype=str)
    # If the try fails, I assume that's because the file has a trackline. Note
    # that I don't preserve the trackline (I'm not sure how pybedtools keeps
    # track of it anyway).
    except pd.parser.CParserError:
        pass
    df = pd.read_table(exons.fn, header=None, skiprows=1, dtype=str)
    out = pbt.BedTool('\n'.join(df[0].apply(lambda x: x[3:]) + '\t' + 
                                df[1] + '\t' + df[2]) + '\n', 
                      from_string=True)
    return out

def add_chr(bt):
    """Add 'chr' to chromosomes for BedTool object

    Parameters
    ----------
    bt : pybedtools.BedTool
        BedTool to add 'chr' to.

    Returns
    -------
    out : pybedtools.BedTool
        New BedTool with 'chr' added to chromosome names.

    """
    try:
        df = pd.read_table(exons.fn, header=None, dtype=str)
    # If the try fails, I assume that's because the file has a trackline. Note
    # that I don't preserve the trackline (I'm not sure how pybedtools keeps
    # track of it anyway).
    except pd.parser.CParserError:
        pass
    df = pd.read_table(exons.fn, header=None, skiprows=1, dtype=str)
    out = pbt.BedTool('\n'.join('chr' + df[0] + '\t' + 
                                df[1] + '\t' + df[2]) + '\n', 
                      from_string=True)
    return out
