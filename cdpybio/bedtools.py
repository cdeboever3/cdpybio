import pandas as pd
import pybedtools as pbt

from general import _sample_names

def beds_to_boolean(beds, ref=None, **kwargs):
    """
    Compare a list of bed files or BedTool objects to a reference bed file and
    create a boolean matrix where each row is an interval and each column is a 1
    if that file has an interval that overlaps the row interval and a 0
    otherwise. If no reference bed is provided, the provided bed files will be
    merged into a single bed and compared to that.

    Parameters
    ----------
    beds : list
        List of paths to bed files or BedTool objects.

    ref : str or BedTool
        Reference bed file to compare against. If no reference bed is provided,
        the provided bed files will be merged into a single bed and compared to
        that.

    names : list of strings
        Names to use for columns of output files. Overrides define_sample_name 
        if provided.

    define_sample_name : function that takes string as input
        Function mapping filename to sample name (or basename). For instance,
        you may have the basename in the path and use a regex to extract it.
        The basenames will be used as the column names. If this is not provided,
        the columns will be named as the input files.

    Returns
    -------
    out : pandas.DataFrame
        Boolean data frame indicating whether each bed file has an interval
        that overlaps each interval in the reference bed file. 

    """
    fns = []
    for i,v in enumerate(beds):
        if type(v) == str:
            fns.append(v)
            beds[i] = pbt.BedTool(v)
        else:
            fns.append(v.fn)
        beds[i] = beds[i].sort()

    names = _sample_names(fns, kwargs)
    if ref:
        if type(ref) == str:
            ref = pbt.BedTool(ref)
        ref = ref.sort()
    else:
        ref = combine(beds)
    
    ind = []
    for r in ref:
        ind.append('{}:{}-{}'.format(r.chrom, r.start, r.stop))
    bdf = pd.DataFrame(0, index=ind, columns=names)
    for i,bed in enumerate(beds):
        res = ref.intersect(bed, sorted=True, wa=True)
        ind = []
        for r in res:
            ind.append('{}:{}-{}'.format(r.chrom,
                                         r.start,
                                         r.stop))
        bdf.ix[ind, names[i]] = 1
    return bdf

def combine(beds, postmerge=True):
    """
    Combine a list of bed files or BedTool objects into a single BedTool object.

    Parameters
    ----------
    beds : list
        List of paths to bed files or BedTool objects.

    postmerge : boolean
        Whether to merge intervals after combining beds together. 

    Returns
    -------
    out : pybedtools.BedTool
        New BedTool with intervals from all input beds.

    """
    for i,v in enumerate(beds):
        if type(v) == str:
            beds[i] = pbt.BedTool(v)
        beds[i] = beds[i].sort()

    out = reduce(lambda x,y : x.cat(y, postmerge=postmerge), beds)
    return out


def write_bed_with_trackline(bed, out, trackline, add_chr=False):
    """
    Read a bed file and write a copy with a trackline. Here's a simple trackline
    example: 'track type=bed name="cool" description="A cool track."'

    Parameters
    ----------
    bed : str 
        Input bed file name.
    out : str
        Output bed file name.
    trackline : str
        UCSC trackline.
    add_chr : boolean
        Add 'chr' to the chromosomes in the input file. Necessary for
        UCSC genome browser if not present.

    """
    df = pd.read_table(bed, index_col=None, header=None)
    bt = pbt.BedTool('\n'.join(df.apply(lambda x: '\t'.join(x.astype(str)), 
                                        axis=1)) + '\n',
                     from_string=True)
    if add_chr:
        bt = add_chr_to_contig(bt)
    bt = bt.saveas(out, trackline=trackline)

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
        df = pd.read_table(bt.fn, header=None, dtype=str)
    # If the try fails, I assume that's because the file has a trackline. Note
    # that I don't preserve the trackline (I'm not sure how pybedtools keeps
    # track of it anyway).
    except pd.parser.CParserError:
        df = pd.read_table(bt.fn, header=None, skiprows=1, dtype=str)
    df[0] = df[0].apply(lambda x: x[3:])
    s = '\n'.join(df.astype(str).apply(lambda x: '\t'.join(x), axis=1)) + '\n'
    out = pbt.BedTool(s, from_string=True)
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
        df = pd.read_table(bt.fn, header=None, dtype=str)
    # If the try fails, I assume that's because the file has a trackline. Note
    # that I don't preserve the trackline (I'm not sure how pybedtools keeps
    # track of it anyway).
    except pd.parser.CParserError:
        df = pd.read_table(bt.fn, header=None, skiprows=1, dtype=str)
    df[0] = 'chr' + df[0]
    s = '\n'.join(df.astype(str).apply(lambda x: '\t'.join(x), axis=1)) + '\n'
    out = pbt.BedTool(s, from_string=True)
    return out
