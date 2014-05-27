import pandas as pd

# Column labels for SJ.out.tab.
COLUMN_NAMES = ('chrom', 'first_bp_intron', 'last_bp_intron', 'strand',
                'intron_motif', 'annotated',
                'unique_junction_reads', 'multimap_junction_reads',
                'max_overhang')
# Subset of columns from SJ.out.tab.
ANNOTATION_COLS = ('chrom', 'first_bp_intron', 'last_bp_intron', 'intron_motif',
                   'annotated')

def _sj_out_junction(row):
    return '{}:{}-{}'.format(row['chrom'], row['first_bp_intron'],
                             row['last_bp_intron'])

def _sj_out_junction_with_strand(row):
    return '{}:{}-{}:{}'.format(row['chrom'], row['first_bp_intron'],
                             row['last_bp_intron'], row['strand'])

def _sj_out_donor(row):
    if row['strand'] == '+':
        return '{}:{}:{}'.format(row['chrom'], row['first_bp_intron'], 
                                 row['strand'])
    if row['strand'] == '-':
        return '{}:{}:{}'.format(row['chrom'], row['last_bp_intron'], 
                                 row['strand'])

def _sj_out_acceptor(row):
    if row['strand'] == '+':
        return '{}:{}:{}'.format(row['chrom'], row['last_bp_intron'], 
                                 row['strand'])
    if row['strand'] == '-':
        return '{}:{}:{}'.format(row['chrom'], row['first_bp_intron'], 
                                 row['strand'])

def read_sj_out_tab(filename):
    """Read an SJ.out.tab file as produced by the RNA-STAR aligner into a
    pandas Dataframe

    Parameters
    ----------
    filename : str of filename or file handle
        Filename of the SJ.out.tab file you want to read in

    Returns
    -------
    sj : pandas.DataFrame
        Dataframe of splice junctions

    """
    def int_to_intron_motif(n):
        if n == 0:
            return 'non-canonical'
        if n == 1:
            return 'GT/AG'
        if n == 2:
            return 'CT/AC'
        if n == 3:
            return 'GC/AG'
        if n == 4:
            return 'CT/GC'
        if n == 5:
            return 'AT/AC'
        if n == 6:
            return 'GT/AT'

    sj = pd.read_table(filename, header=None, names=COLUMN_NAMES)
    sj.intron_motif = sj.intron_motif.map(int_to_intron_motif)
    sj.annotated = sj.annotated.map(bool)
    return sj

def make_sj_out_dict(fns, define_sample_name=None):
    """Read multiple sj_outs, return dict with keys as sample names and values
    as sj_out dataframes

    Parameters
    ----------
    fns : list of strs of filenames or file handles
        List of filename of the SJ.out.tab files to read in

    define_sample_name : function that takes string as input
        Function mapping filename to sample name. For instance, you may have the
        sample name in the path and use a regex to extract it.  The sample names
        will be used as the column names. If this is not provided, the columns
        will be named as the input files.

    Returns
    -------
    sj_outD : dict
        Dict whose keys are sample names and values are sj_out dataframes
    
    """
    if define_sample_name == None:
        define_sample_name = lambda x: x
    else:
        assert len(set([ define_sample_name(x) for x in fns ])) == len(fns)
    sj_outD = dict()
    for k in sj_outD.keys():
        df = sj_outD[k]
        sj_outD[k] = df[df.unique_junction_reads > 0]

    for fn in fns:
        sample = define_sample_name(fn)
        df = read_sj_out_tab(fn)
        index = df.apply(lambda x: _sj_out_junction(x), axis=1)
        assert len(index) == len(set(index))
        df.index = index
        sj_outD[sample] = df
    return sj_outD

def make_sj_out_panel(sj_outD, total_jxn_cov_cutoff=20, statsfile=None):
    """Filter junctions from many sj_out files and make panel

    Parameters
    ----------
    sj_outD : dict
        Dict whose keys are sample names and values are sj_out dataframes

    total_jxn_cov_cutoff : int
        If the unique read coverage of a junction summed over all samples is not
        greater than or equal to this value, the junction will not be included
        in the final output.

    statsfile : filename str
        If provided, some stats will be printed to this file

    Returns
    -------
    sj_outP : pandas.Panel
        Panel where each dataframe(?) corresponds to an sj_out file filtered to
        remove low coverage junctions.

    annotDF : pandas.DataFrame
        Dataframe with values ANNOTATION_COLS = ['chrom', 'first_bp_intron', 
        'last_bp_intron', 'intron_motif', 'annotated'] that are otherwise
        duplicated in the panel.
    
    """
    # Remove any junctions that don't have any uniquely mapped junction reads.
    # Even if a junction passes the cutoff in other samples, we are only
    # concerned with unique counts.
    # set of all junctions
    jxnS = reduce(lambda x,y: set(x) | set(y),
                  [ sj_outD[k].index for k in sj_outD.keys() ])

    jxn_keepS = set()
    for j in jxnS:
        if sum([ sj_outD[k].ix[j,'unique_junction_reads'] for k in sj_outD.keys()
                 if j in sj_outD[k].index ]) >= total_jxn_cov_cutoff:
            jxn_keepS.add(j)

    for k in sj_outD.keys():
        sj_outD[k] = sj_outD[k].ix[jxn_keepS]

    sj_outP = pd.Panel(sj_outD)
    for col in ['unique_junction_reads', 'multimap_junction_reads',
                'max_overhang']:
        sj_outP.ix[:,:,col] = sj_outP.ix[:,:,col].fillna(0)

    # Some dataframes will be missing information like intron_motif etc. for 
    # junctions that were not observed in that sample. We'll add that info.
    annotDF = reduce(pd.DataFrame.combine_first,
                     [ sj_outP.ix[item,:,ANNOTATION_COLS].dropna() for item in
                      sj_outP.items ])
    for col in ANNOTATION_COLS:
        sj_outP.ix[:,annotDF.index,col] = [ annotDF[col] for i in
                                               range(sj_outP.shape[0]) ]
    
    # All of the splice junction annotation information is duplicated in each
    # dataframe in the panel, so we'll make a single dataframe holding that
    # information.
    annotDF = sj_outP.ix[0,:,ANNOTATION_COLS]

    if statsfile:
        statsF = open(statsfile,'w')

        statsF.write('Number of junctions in sj_out file per sample\n')
        for k in sj_outD.keys():
            statsF.write('{0}\t{1:,}\n'.format(k,sj_outD[k].shape[0]))
        statsF.write('\n')
  
        statsF.write('sj_out panel size\t{0}\n\n'.format(sj_outP.shape))
        statsF.close()
    return sj_outP,annotDF

def read_external_annotation(fn, statsfile=None):
    """Read file with junctions from some database. This does not have to be the
    same splice junction database used with STAR

    Parameters
    ----------
    fn : filename str
        File with splice junctions from annotation. The file should have a
        header and contained the following columns  'gene', 'chrom', 'start',
        'end', 'strand', 'chr:start', 'chr:end', 'donor', 'acceptor', 'intron'.

    statsfile : string
        File to write statistics to.

    Returns
    -------
    extDF : pandas.DataFrame
        DataFrame indexed by splice junction
    
    """
    extDF = pd.read_table(fn, header=0)
    
    extDF['junction_no_strand'] = extDF.junction.apply(
        lambda x: juncRE.match(x).group().strip(':'))

    # In rare cases, a splice junction might be used by more than one gene. For
    # my purposes, these cases are confounding, so I will remove all such splice
    # junctions. 
    junctions_to_keepSE = extDF.junction_no_strand.value_counts() == 1
    # Drop duplicate junctions, but this leaves one copy of the duplicate
    # junction, so we will reindex and keep only the junctions we want.
    extDF = extDF.drop_duplicates(cols='junction_no_strand')
    extDF.index = extDF.junction_no_strand
    extDF = extDF.ix[junctions_to_keepSE]
    # Reindex with strand info.
    extDF.index = extDF.junction
    extDF = extDF.drop('junction', axis=1)

    if statsfile:
        f = open(statsfile, 'w')
        f.write('Read external annotation\t{}\n'.format(fn))
        f.write('Total number of junctions\t{}\n'.format(total_num))
        f.write('''Number of junctions used in only one
                gene\t{}'''.format(extDF.shape[0]))
        f.close()

    return extDF

def filter_jxns_donor_acceptor(sj_outP, annotDF, extDF, statsfile=None):
    """Remove junctions that do not use an annotated donor or acceptor according
    to external junction annotation. Add strand and gene information for
    junctions according to external annoation (STAR strand ignored).

    Parameters
    ----------
    sj_outP : pandas.Panel
        Panel where each dataframe(?) corresponds to an sj_out file filtered to
        remove low coverage junctions.

    annotDF : pandas.DataFrame
        Annotation information from STAR. Should have the following columns:
        'chrom', 'first_bp_intron', 'last_bp_intron', 'intron_motif',
        'annotated'.

    extDF : pandas.DataFrame
        Dataframe containing information about annotated splice sites. These can
        differ from those provided to STAR for alignment. The dataframe should
        have the following columns: 'gene', 'chrom', 'start', 'end', 'strand',
        'chr:start', 'chr:end', 'donor', 'acceptor', 'intron'

    Returns
    -------
    countDF :  pandas.DataFrame
        Number of unique junction spanning reads for each junction that passed
        filtering criteria.

    annotDF : pandas.DataFrame
        Annotation information for junctions that passed filtering criteria.
    
    """
    import re

    sjRE = re.compile('(.*:.*-.*):(\+|-)')
    juncRE = re.compile('(.*):(\d*)-(\d*):') 
    
    # Add column showing whether junction is in external annotation.
    annotDF['ext_annotated'] = False
    annotDF.ix[set(annotDF.index) & set(extDF.junction_no_strand),
               'ext_annotated'] = True

    # Add strand information to annotation of STAR junctions that are in
    # external database.
    strandSE = pd.Series(extDF.strand.values,index=extDF.junction_no_strand)
    strandSE = strandSE[set(strandSE.index) & set(annotDF.index)]
    annotDF['strand'] = '*'
    annotDF.ix[strandSE.index,'strand'] = strandSE.values

    # Add column for start and end location (chromosome plus position for
    # uniqueness).
    annotDF['chr:start'] = annotDF.apply(
        lambda x: '{}:{}'.format(x['chrom'],x['first_bp_intron']),axis=1)
    annotDF['chr:end'] = annotDF.apply(
        lambda x: '{}:{}'.format(x['chrom'],x['last_bp_intron']),axis=1)

    ext_startS = set(extDF['chr:start'].values)
    ext_endS = set(extDF['chr:end'].values)

    # Remove junctions that don't have a start or end shared with external
    # database.
    junctions_to_removeSE = annotDF[annotDF.ext_annotated == False].apply(
        lambda x: (x['chr:start'] in ext_startS) + 
        (x['chr:end'] in ext_endS) == 0,axis=1)
    annotDF = annotDF.drop(junctions_to_removeSE[junctions_to_removeSE].index)

    # Add column indicating which gene the junctions belong to for annotated
    # jxn's.
    geneSE = pd.Series(dict(zip(extDF.junction_no_strand.values,extDF.gene)))
    annotDF['gene_id'] = ''
    annotDF['gene_id'] = geneSE[annotDF.index]

    # Now we'll figure out the genes for the junctions that aren't in our
    # database. We can associate each start and end with a gene and use this.
    start_geneSE = pd.Series(
        dict(zip(extDF.apply(lambda x: '{}:{}'.format(
            x['chrom'],x['first_bp_intron']),axis=1),extDF.gene)))
    end_geneSE = pd.Series(
        dict(zip(extDF.apply(lambda x: '{}:{}'.format(
            x['chrom'],x['last_bp_intron']),axis=1),extDF.gene)))

    for ind in annotDF[annotDF.ext_annotated == False].index:
        cur_start = annotDF.ix[ind,'chr:start'] 
        if cur_start in start_geneSE.index:
            annotDF.ix[ind,'gene_id'] = start_geneSE[cur_start]
    for ind in annotDF[annotDF.ext_annotated == False].index:
        cur_end = annotDF.ix[ind,'chr:end'] 
        if cur_end in end_geneSE.index:
            annotDF.ix[ind,'gene_id'] = end_geneSE[cur_end]

    # We can use the genes to assign strand to the novel splice junctions.
    strandSE = pd.Series(dict(zip(extDF.gene,extDF.strand)))
    for ind in annotDF[annotDF.ext_annotated == False].index:
        annotDF.ix[ind,'strand'] = strandSE[annotDF.ix[ind,'gene_id']]

    # And re-index with the strand info.
    annotDF.index = [ x + ':' + annotDF.ix[x,'strand'] for x in annotDF.index ]

    # Now we'll add donor and acceptor info.
    annotDF['donor'] = annotDF.apply(lambda x: _sj_out_donor(x),axis=1)
    annotDF['acceptor'] = annotDF.apply(lambda x: _sj_out_acceptor(x),axis=1)

    # And whether the donor or acceptor is in the external database or not.
    extDF['donor'] = extDF.apply(lambda x: _sj_out_donor(x),axis=1)
    extDF['acceptor'] = extDF.apply(lambda x: _sj_out_acceptor (x),axis=1)
    ext_donorS = set(extDF.donor)
    ext_acceptorS = set(extDF.acceptor)
    annotDF['novel_donor'] = False
    annotDF['novel_acceptor'] = False
    for ind in annotDF[annotDF.ext_annotated == False].index:
        annotDF.ix[ind,'novel_donor'] = (annotDF.ix[ind,'donor'] not in 
                                         ext_donorS)
        annotDF.ix[ind,'novel_acceptor'] = (annotDF.ix[ind,'acceptor'] not in
                                            ext_acceptorS)

    # Sort by gene ID and start/end.
    annotDF = annotDF.sort(columns=['gene_id', 'first_bp_intron', 
                                    'last_bp_intron'])

    # Make file with counts for the junctions we are interested in.
    L = [ juncRE.match(x).group().strip(':') for x in annotDF.index ]
    countDF = sj_out_filteredP.ix[:, L, 'unique_junction_reads']
    countDF.index = annotDF.index

    if statsfile:
        f = open(statsfile, 'w')

        t = extDF.shape[0]
        f.write('''Number of junctions in external 
                annotation\t{0:,}\n'''.format(t))

        t = annotDF.annotated.sum()
        f.write('''Number observed junctions in STAR
                annotation\t{0:,}\n'''.format(t))

        t = annotDF.ext_annotated.sum()
        f.write('''Number observed junctions in external
                annotation\t{0:,}\n'''.format(t))

        t = annotDF.shape[0] - annotDF.ext_annotated.sum()
        f.write('''Number of observed junctions not in external
                annotation\t{0:,}\n'''.format(t))

        t = annotDF.ix[annotDF.ext_annotated == False,'annotated'].sum()
        f.write('''Number of observed junctions not in external annotation but
                in STAR annotation\t{0:,}\n'''.format(t))

        t = sum(annotDF.ix[annotDF.ext_annotated == False,'annotated'].values == 0)
        f.write('''Number of observed junctions not in external annotation and
                not in STAR annotation\t{0:,}\n'''.format(t))
  
        t = len(set(annotDF[annotDF.novel_donor].donor))
        f.write('''Number of novel donors\t{0:,}\n'''.format(t))

        t = annotDF.novel_donor.sum()
        f.write('''Number of novel junctions with novel
                donors\t{0:,}\n'''.format(t))

        t = len(set(annotDF[annotDF.novel_acceptor].acceptor))
        f.write('''Number of novel acceptors\t{0:,}\n'''.format(t))

        t = annotDF.novel_acceptor.sum()
        f.write('''Number of novel junctions with novel
                acceptors\t{0:,}\n'''.format(t))

        t = (annotDF[annotDF.ext_annotated].shape[0] - 
             sum(annotDF.novel_donor) - 
             sum(annotDF.novel_acceptor))
        f.write('''Number of novel junctions with new combination of donor and
                acceptor\t{0:,}\n'''.format(t))

        f.close()
   
    return countDF, annotDF

def combine_sj_out(fns, external_db, total_jxn_cov_cutoff=20, 
                   define_sample_name=None, statsfile=None):
    """Combine SJ.out.tab files from STAR by filtering based on coverage and
    comparing to an external annotation to discover novel junctions.

    Parameters
    ----------
    fns : list of strings 
        Filenames of SJ.out.tab files to combine.

    external_db : str
        Filename of splice junction information from external database. The file
        should have a header and contained the following columns  'gene',
        'chrom', 'start', 'end', 'strand', 'chr:start', 'chr:end', 'donor',
        'acceptor', 'intron'.

    total_jxn_cov_cutoff : int
        Discard junctions with less than this many reads summed over all
        samples.

    define_samle_name : function
        A function mapping the SJ.out.tab filenames to sample names.

    statsfile : string
        File to write statistics to. Statistics are not output if a filename is
        not provided.

    Returns
    -------
    countDF :  pandas.DataFrame
        Number of unique junction spanning reads for each junction that passed
        filtering criteria.

    annotDF : pandas.DataFrame
        Annotation information for junctions that passed filtering criteria.
    
    """
    sj_outD = make_sj_out_dict(fns, define_sample_name=define_sample_name)
    sj_outP, annotDF = make_sj_out_panel(sj_outD, total_jxn_cov_cutoff, 
                                        statsfile=statsfile)

    # I'll read in the statsfile and hold that information so I can combine into
    # a final statsfile.
    f = open(statsfile, 'r')
    lines = f.readlines()
    f.close()
    
    extDF = read_external_annotation(external_db)
    countsDF, annotDF = filter_jxns_donor_acceptor(sj_outP, annotDF, extDF, 
                                                   statsfile=statsfile)

    f = open(statsfile, 'r')
    lines2 = f.readlines()
    f.close()

    f = open(statsfile, 'w')
    f.write(lines + lines2)
    f.close()

    return countsDF, annotDF