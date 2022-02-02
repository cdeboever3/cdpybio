import pandas as pd
import pybedtools as pbt

def read_variants(fns, remove=['DBSNP'], keep_only=True,
                  min_tumor_f=0.1, min_tumor_cov=14,
                  min_normal_cov=8):
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

    min_tumor_f : float between 0 and 1
        Minimum tumor allelic fraction.

    min_tumor_cov : int > 0
        Minimum coverage of the variant in the tumor.

    min_normal_cov : int > 0
        Minimum coverage of the variant in the normal.

    """
    variants = []
    for i, f in enumerate(fns):
        # If keep_only, use awk to only grab those lines for big speedup.
        if keep_only:
            from numpy import dtype
            import subprocess
            res = subprocess.check_output(
                'awk \'$35 == "KEEP"\' {}'.format(f), shell=True)
            if res.strip() != '': 
                columns = ['contig', 'position', 'context', 'ref_allele',
                           'alt_allele', 'tumor_name', 'normal_name',
                           'score', 'dbsnp_site', 'covered', 'power',
                           'tumor_power', 'normal_power', 'total_pairs',
                           'improper_pairs', 'map_Q0_reads', 't_lod_fstar',
                           'tumor_f', 'contaminant_fraction',
                           'contaminant_lod', 't_ref_count', 't_alt_count',
                           't_ref_sum', 't_alt_sum', 't_ref_max_mapq',
                           't_alt_max_mapq', 't_ins_count', 't_del_count',
                           'normal_best_gt', 'init_n_lod', 'n_ref_count',
                           'n_alt_count', 'n_ref_sum', 'n_alt_sum',
                           'judgement']
                tdf = pd.DataFrame(
                    [x.split('\t') for x in res.strip().split('\n')],
                    columns=columns)
                tdf = tdf.convert_objects(convert_numeric=True)
            else:
                tdf = pd.DataFrame(columns=columns)
            tdf['contig'] = tdf.contig.astype(object)

        else:
            tdf = pd.read_table(f, index_col=None, header=0, skiprows=1,
                                low_memory=False, 
                                dtype={'contig':object})
        for t in remove:
            tdf = tdf[tdf.dbsnp_site != t]
        tdf = tdf[tdf.tumor_f > min_tumor_f]
        tdf = tdf[tdf.t_ref_count + tdf.t_alt_count > min_tumor_cov]
        tdf = tdf[tdf.n_ref_count + tdf.n_alt_count > min_normal_cov]
        variants.append(tdf)
    variants = pd.concat(variants)
    variants.index = list(range(variants.shape[0]))
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
