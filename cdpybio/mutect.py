import pandas as pd
import pybedtools as pbt

def read_variants(fns, remove=['DBSNP'], keep_only=True ):
    """Read muTect results from the list of files fns

    Parameters
    ----------
    fns : list
        List of MuTect output files.

    Returns
    -------
    variants : pandas.DataFrame
        Pandas DataFrame summarizing variant calling results.

    remove : list
        List of site types for column "dbsnp_site" to remove.

    keep_only : boolean
        If True, only keep variants with 'KEEP' in "judgement" column.
        Otherwise, keep all variants.

    """
    import numpy.object as npobject
    variants = []
    for i, f in enumerate(fns):
        tdf = pd.read_table(f, index_col=None, header=0, skiprows=1,
                            low_memory=False, 
                            dtype={'contig':npobject})
        for t in remove:
            tdf = tdf[tdf.dbsnp_site != t]
        if keep_only:
            tdf = tdf[tdf.judgement == 'KEEP']
        variants.append(tdf)
    variants = pd.concat(variants)
    variants.index = range(variants.shape[0])
    return variants

def mutect_to_bed(df):
    """Convert MuTect results (read into dataframe) to BedTool object

    Parameters
    ----------
    df : pandas.DataFrame
        Pandas DataFrame with MuTect results.

    Returns
    -------
    bt : pybedtools.BedTool
        BedTool with variants.

    """
    s = (df.contig.astype(str) + '\t' + 
         (df.position - 1).astype(int).astype(str) + '\t' + 
         df.position.astype(int).astype(str) + '\t' + 
         df.tumor_name)
    s = '\n'.join(s.values) + '\n'
    bt = pbt.BedTool(s, from_string=True)
    return bt
