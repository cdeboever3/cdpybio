import copy
import os
import pdb

import numpy as np
import pandas as pd

# Column labels for SJ.out.tab.
COLUMN_NAMES = ('chrom', 'start', 'end', 'strand',
                'intron_motif', 'annotated',
                'unique_junction_reads', 'multimap_junction_reads',
                'max_overhang')
# Subset of columns from SJ.out.tab that are quantitative and change sample to
# sample.
COUNT_COLS = ('unique_junction_reads', 'multimap_junction_reads',
              'max_overhang')
# Subset of columns from SJ.out.tab that do not change sample to sample because
ANNOTATION_COLS = ('chrom', 'start', 'end', 'strand',
                   'intron_motif', 'annotated')
# TODO: add strand to above list. I don't use the STAR strand because I didn't
# understand it (maybe buggy?) but I'm sure it's fixed in recent versions.

def _sj_out_junction(row):
    return '{}:{}-{}'.format(row['chrom'], row['start'],
                             row['end'])

def _sj_out_junction_with_strand(row):
    return '{}:{}-{}:{}'.format(row['chrom'], row['start'],
                             row['end'], row['strand'])

def _sj_out_donor(row):
    if row['strand'] == '+':
        return '{}:{}:{}'.format(row['chrom'], row['start'], 
                                 row['strand'])
    if row['strand'] == '-':
        return '{}:{}:{}'.format(row['chrom'], row['end'], 
                                 row['strand'])

def _sj_out_acceptor(row):
    if row['strand'] == '+':
        return '{}:{}:{}'.format(row['chrom'], row['end'], 
                                 row['strand'])
    if row['strand'] == '-':
        return '{}:{}:{}'.format(row['chrom'], row['start'], 
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
    sj.strand.astype('object')
    sj.strand = sj.strand.apply(lambda x: ['unk','+','-'][x])
    # See https://groups.google.com/d/msg/rna-star/B0Y4oH8ZSOY/NO4OJbbUU4cJ for
    # definition of strand in SJout files.
    sj = sj.sort(columns=['chrom', 'start', 'end'])
    return sj

def _make_sj_out_dict(fns, define_sample_name=None):
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

    for fn in fns:
        sample = define_sample_name(fn)
        df = read_sj_out_tab(fn)
        # Remove any junctions that don't have any uniquely mapped junction
        # reads.  Even if a junction passes the cutoff in other samples, we are
        # only concerned with unique counts.
        df = df[df.unique_junction_reads > 0]
        index = df.chrom + ':' + df.start.astype(str) + '-' +  df.end.astype(str)
        assert len(index) == len(set(index))
        df.index = index
        sj_outD[sample] = df

    return sj_outD

def _make_sj_out_panel(sj_outD, total_jxn_cov_cutoff=20):
    """Filter junctions from many sj_out files and make panel

    Parameters
    ----------
    sj_outD : dict
        Dict whose keys are sample names and values are sj_out dataframes

    total_jxn_cov_cutoff : int
        If the unique read coverage of a junction summed over all samples is not
        greater than or equal to this value, the junction will not be included
        in the final output.

    Returns
    -------
    sj_outP : pandas.Panel
        Panel where each dataframe corresponds to an sj_out file filtered to
        remove low coverage junctions. Each dataframe has COUNT_COLS =
        ('unique_junction_reads', 'multimap_junction_reads', 'max_overhang')

    annotDF : pandas.DataFrame
        Dataframe with values ANNOTATION_COLS = ('chrom', 'start', 
        'end', 'intron_motif', 'annotated') that are otherwise
        duplicated in the panel.
    
    """
    num_jxns = dict()
    # set of all junctions
    jxnS = reduce(lambda x,y: set(x) | set(y),
                  [ sj_outD[k].index for k in sj_outD.keys() ])

    jxn_keepS = set()
    jxn_setsD = dict()
    for k in sj_outD.keys():
        jxn_setsD[k] = frozenset(sj_outD[k].index)
    for j in jxnS:
        if sum([ sj_outD[k].ix[j,'unique_junction_reads'] for k in sj_outD.keys()
                 if j in jxn_setsD[k] ]) >= total_jxn_cov_cutoff:
            jxn_keepS.add(j)

    for k in sj_outD.keys():
        sj_outD[k] = sj_outD[k].ix[jxn_keepS]

    sj_outP = pd.Panel(sj_outD)
    for col in ['unique_junction_reads', 'multimap_junction_reads',
                'max_overhang']:
        sj_outP.ix[:,:,col] = sj_outP.ix[:,:,col].fillna(0)

    # Some dataframes will be missing information like intron_motif etc. for 
    # junctions that were not observed in that sample. The info is somewhere in
    # the panel though so we can get it.
    annotDF = reduce(pd.DataFrame.combine_first,
                     [ sj_outP.ix[item,:,ANNOTATION_COLS].dropna() for item in
                      sj_outP.items ])
    annotDF['start'] = annotDF['start'].astype(int)
    annotDF['end'] = annotDF['end'].astype(int)
    annotDF['annotated'] = annotDF['annotated'].astype(bool)
    # Sort annotation and panel
    annotDF = annotDF.sort(columns=['chrom', 'start', 'end'])
    sj_outP = sj_outP.ix[:, annotDF.index, :]

    sj_outP = sj_outP.ix[:,:,COUNT_COLS].astype(int)

    return sj_outP, annotDF

def read_external_annotation(fn):
    """Read file with junctions from some database. This does not have to be the
    same splice junction database used with STAR

    Parameters
    ----------
    fn : filename str
        File with splice junctions from annotation. The file should have a
        header and contained the following columns  'gene', 'chrom', 'start',
        'end', 'strand', 'chrom:start', 'chrom:end', 'donor', 'acceptor', 
        'intron'.

    Returns
    -------
    extDF : pandas.DataFrame
        DataFrame indexed by splice junction

    stats : list of strings
        Human readable statistics about the external database.
    
    """
    assert os.path.exists(fn)
    extDF = pd.read_table(fn, index_col=0, header=0)
    total_num = extDF.shape[0]
    
    # In rare cases, a splice junction might be used by more than one gene. For
    # my purposes, these cases are confounding, so I will remove all such splice
    # junctions. 
    intron_count = extDF.intron.value_counts()
    extDF['intron_count'] = extDF.intron.apply(lambda x: intron_count.ix[x])
    extDF = extDF[extDF.intron_count == 1]
    extDF = extDF.drop('intron_count', axis=1)

    stats = []
    stats.append('External database stats')
    stats.append('Read external annotation\t{}'.format(fn))
    stats.append('Total number of junctions\t{:,}'.format(total_num))
    stats.append(('Number of junctions used in only one '
                  'gene\t{:,}').format(extDF.shape[0]))

    return extDF, stats

def _filter_jxns_donor_acceptor(sj_outP, annotDF, extDF):
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
        'chrom', 'start', 'end', 'intron_motif',
        'annotated'.

    extDF : pandas.DataFrame
        Dataframe containing information about annotated splice sites. These can
        differ from those provided to STAR for alignment. The dataframe should
        have the following columns: 'gene', 'chrom', 'start', 'end', 'strand',
        'chrom:start', 'chrom:end', 'donor', 'acceptor', 'intron'

    Returns
    -------
    countDF :  pandas.DataFrame
        Number of unique junction spanning reads for each junction that passed
        filtering criteria.

    annotDF : pandas.DataFrame
        Annotation information for junctions that passed filtering criteria.
    
    stats : list of strings
        Human readable statistics about the combining process.
    
    """
    import re

    sjRE = re.compile('(.*:.*-.*):(\+|-)')
    juncRE = re.compile('(.*):(\d*)-(\d*):') 
    
    # Add column showing whether junction is in external annotation.
    annotDF = copy.deepcopy(annotDF)
    annotDF['ext_annotated'] = False
    annotDF.ix[set(annotDF.index) & set(extDF.intron),
               'ext_annotated'] = True
    
    # Replace strand information from STAR with strand information from external
    # database.
    strandSE = pd.Series(extDF.strand.values,index=extDF.intron)
    strandSE = strandSE[set(strandSE.index) & set(annotDF.index)]
    annotDF['strand'] = '*'
    annotDF.ix[strandSE.index,'strand'] = strandSE.values

    # Add column for start and end location (chromosome plus position for
    # uniqueness).
    annotDF['chrom:start'] = annotDF.apply(
        lambda x: '{}:{}'.format(x['chrom'],x['start']),axis=1)
    annotDF['chrom:end'] = annotDF.apply(
        lambda x: '{}:{}'.format(x['chrom'],x['end']),axis=1)

    ext_startS = set(extDF['chrom:start'].values)
    ext_endS = set(extDF['chrom:end'].values)

    # Remove junctions that don't have a start or end shared with external
    # database.
    junctions_to_removeSE = annotDF[annotDF.ext_annotated == False].apply(
        lambda x: (x['chrom:start'] in ext_startS) + 
        (x['chrom:end'] in ext_endS) == 0,axis=1)
    if junctions_to_removeSE.shape[0] > 0:
        annotDF = annotDF.drop(junctions_to_removeSE[junctions_to_removeSE].index)

    # Add column indicating which gene the junctions belong to for annotated
    # jxn's.
    geneSE = pd.Series(dict(zip(extDF.intron.values,extDF.gene)))
    annotDF['gene_id'] = ''
    annotDF['gene_id'] = geneSE[annotDF.index]

    # Now we'll figure out the genes for the junctions that aren't in our
    # database. We can associate each start and end with a gene and use this.
    start_gene = dict(zip(extDF['chrom:start'], extDF.gene))
    end_gene = dict(zip(extDF['chrom:end'], extDF.gene))
    
    genes = []
    for ind in annotDF[annotDF.ext_annotated == False].index:
        cur_start = annotDF.ix[ind,'chrom:start'] 
        cur_end = annotDF.ix[ind,'chrom:end'] 
        gene = start_gene.get(cur_start, '')
        if gene == '':
            gene = end_gene.get(cur_end, '')
        genes.append(gene)
    annotDF.ix[annotDF.ext_annotated == False, 'gene_id'] = genes

    # We can use the genes to assign strand to the novel splice junctions.
    strandSE = pd.Series(dict(zip(extDF.gene,extDF.strand)))
    ind = annotDF[annotDF.ext_annotated == False].index
    annotDF.ix[ind, 'strand'] = strandSE[annotDF.ix[ind, 'gene_id']].values

    # And re-index with the strand info.
    annotDF.index = [ x + ':' + annotDF.ix[x,'strand'] for x in annotDF.index ]

    # Now we'll add donor and acceptor info.
    pos_ind = annotDF[annotDF.strand == '+'].index
    neg_ind = annotDF[annotDF.strand == '-'].index
    annotDF.ix[pos_ind, 'donor'] = (annotDF.ix[pos_ind, 'chrom'] + ':' + 
                                    annotDF.ix[pos_ind, 'start'].astype(str) + 
                                    ':' +  annotDF.ix[pos_ind, 'strand'])
    annotDF.ix[pos_ind, 'acceptor'] = (annotDF.ix[pos_ind, 'chrom'] + ':' + 
                                       annotDF.ix[pos_ind, 'end'].astype(str) 
                                       + ':' + annotDF.ix[pos_ind, 'strand'])
    annotDF.ix[neg_ind, 'acceptor'] = (annotDF.ix[neg_ind, 'chrom'] + ':' + 
                                       annotDF.ix[neg_ind, 'start'].astype(str) 
                                       + ':' + annotDF.ix[neg_ind, 'strand'])
    annotDF.ix[neg_ind, 'donor'] = (annotDF.ix[neg_ind, 'chrom'] + ':' + 
                                    annotDF.ix[neg_ind, 'end'].astype(str) + 
                                    ':' + annotDF.ix[neg_ind, 'strand'])

    # And whether the donor or acceptor is in the external database or not.
    ext_donorS = frozenset(extDF.donor)
    ext_acceptorS = frozenset(extDF.acceptor)
    annotDF['novel_donor'] = False
    annotDF['novel_acceptor'] = False
    ind = annotDF[annotDF.ext_annotated == False].index
    if len(ind) > 0:
        novel_donor = []
        novel_acceptor = []
        for i in ind:
            novel_donor.append(annotDF.ix[i, 'donor'] not in ext_donorS)
            novel_acceptor.append(annotDF.ix[i, 'acceptor'] not in ext_acceptorS)
        annotDF.ix[ind, 'novel_donor'] = novel_donor
        annotDF.ix[ind, 'novel_acceptor'] = novel_acceptor

    # Sort by gene ID and start/end.
    annotDF = annotDF.sort(columns=['gene_id', 'start', 'end'])

    # Make file with counts for the junctions we are interested in.
    L = [ juncRE.match(x).group().strip(':') for x in annotDF.index ]
    countDF = sj_outP.ix[:, L, 'unique_junction_reads']
    countDF.index = annotDF.index

    stats = []
    t = extDF.shape[0]
    stats.append('Junction filtering stats')
    stats.append(('Number of junctions in external '
                  'annotation\t{0:,}').format(t))

    t = annotDF.annotated.sum()
    stats.append(('Number observed junctions in STAR '
                  'annotation\t{0:,}').format(t))

    t = annotDF.ext_annotated.sum()
    stats.append(('Number observed junctions in external '
                  'annotation\t{0:,}').format(t))

    t = annotDF.shape[0] - annotDF.ext_annotated.sum()
    stats.append(('Number of observed junctions not in external '
                  'annotation\t{0:,}').format(t))

    t = annotDF.ix[annotDF.ext_annotated == False,'annotated'].sum()
    stats.append(('Number of observed junctions not in external annotation but '
                  'in STAR annotation\t{0:,}').format(t))

    t = sum(annotDF.ix[annotDF.ext_annotated == False,'annotated'].values == 0)
    stats.append(('Number of observed junctions not in external annotation and '
                  'not in STAR annotation\t{0:,}').format(t))
  
    t = len(set(annotDF.ix[annotDF.novel_donor, 'donor']))
    stats.append(('Number of novel donors\t{0:,}').format(t))

    t = annotDF.novel_donor.sum()
    stats.append(('Number of novel junctions with novel '
                  'donors\t{0:,}').format(t))

    t = len(set(annotDF.ix[annotDF.novel_acceptor, 'acceptor']))
    stats.append(('Number of novel acceptors\t{0:,}').format(t))

    t = annotDF.novel_acceptor.sum()
    stats.append(('Number of novel junctions with novel '
                  'acceptors\t{0:,}').format(t))

    t = (annotDF[annotDF.ext_annotated == False].shape[0] - 
         sum(annotDF.novel_donor) - 
         sum(annotDF.novel_acceptor))
    stats.append(('Number of novel junctions with new combination of donor and '
                  'acceptor\t{0:,}').format(t))
    return countDF, annotDF, stats

def combine_sj_out(fns, external_db, total_jxn_cov_cutoff=20, 
                   define_sample_name=None):
    """Combine SJ.out.tab files from STAR by filtering based on coverage and
    comparing to an external annotation to discover novel junctions.

    Parameters
    ----------
    fns : list of strings 
        Filenames of SJ.out.tab files to combine.

    external_db : str
        Filename of splice junction information from external database. The file
        should have a header and contained the following columns  'gene',
        'chrom', 'start', 'end', 'strand', 'chrom:start', 'chrom:end', 'donor',
        'acceptor', 'intron'.

    total_jxn_cov_cutoff : int
        Discard junctions with less than this many reads summed over all
        samples.

    define_sample_name : function
        A function mapping the SJ.out.tab filenames to sample names.

    Returns
    -------
    countDF :  pandas.DataFrame
        Number of unique junction spanning reads for each junction that passed
        filtering criteria.

    annotDF : pandas.DataFrame
        Annotation information for junctions that passed filtering criteria.

    stats : list of strings
        Human readable statistics.
    
    """
    stats = []
    sj_outD = _make_sj_out_dict(fns, define_sample_name=define_sample_name)
    stats.append('Number of junctions in SJ.out file per sample')
    for k in sj_outD.keys():
        stats.append('{0}\t{1:,}'.format(k, sj_outD[k].shape[0]))
    stats.append('')

    sj_outP, annotDF = _make_sj_out_panel(sj_outD, total_jxn_cov_cutoff)
    stats.append('SJ.out panel size\t{0}'.format(sj_outP.shape))
    stats.append('')

    extDF, ext_stats = read_external_annotation(external_db)
    stats += ext_stats
    stats.append('')

    countsDF, annotDF, filter_stats = _filter_jxns_donor_acceptor(sj_outP, 
                                                                  annotDF, 
                                                                  extDF)
    stats += filter_stats
    return countsDF, annotDF, stats

def _make_splice_targets_dict(df, feature, strand):
    """Make dict mapping each donor to the location of all acceptors it splices
    to or each acceptor to all donors it splices from.

    Parameters
    ----------
    df : pandas.DataFrame 
        Dataframe with splice junction information from external database
        containing columns 'gene', 'chrom', 'start', 'end', 'strand',
        'chrom:start', 'chrom:end', 'donor', 'acceptor', 'intron'.

    feature : string
        Either 'donor' or 'acceptor'.

    strand : string
        Either '+' or '-'.

    Returns
    -------
    d : dict
        If feature='donor', dict whose keys are all distinct donors in df and
        whose values are the distinct locations (integers) of the acceptors that
        donor splices to in a numpy array. If feature='acceptor', dict whose
        keys are all distinct acceptors in df and whose values are the distinct
        locations (integers) of the donors that acceptor splices from in a numpy
        array.
    
    """
    g = df[df.strand == strand].groupby(feature)
    d = dict()
    if strand == '+':
        if feature == 'donor':
            target = 'end'
        if feature == 'acceptor':
            target = 'start'
    if strand == '-':
        if feature == 'donor':
            target = 'start'
        if feature == 'acceptor':
            target = 'end'

    for k in g.groups.keys():
        d[k] = np.array(list(set(df.ix[g.groups[k], target])))
        d[k].sort()
    return d

def find_novel_donor_acceptor_dist(annot, ext):
    """Find nearest annotated upstream/downstream donor/acceptor for novel
    donor/acceptors 

    Parameters
    ----------
    annot : pandas.DataFrame
        Dataframe with observed splice junctions (novel and known) with columns
        'chrom', 'first_bp_intron', 'last_bp_intron', 'strand', 'intron_motif',
        'annotated', 'ext_annotated', 'chrom:start', 'chrom:end', 'gene_id',
        'donor', 'acceptor', 'novel_donor', 'novel_acceptor'.

    ext : pandas.DataFrame 
        Dataframe with splice junction information from external database
        containing columns 'gene', 'chrom', 'start', 'end', 'strand',
        'chrom:start', 'chrom:end', 'donor', 'acceptor', 'intron'.

    Returns
    -------
    annot : pandas.DataFrame
        Annotation information with added columns for distance from novel
        donor/acceptor to nearest upstream/downstream donor/acceptor.
    
    """
    pos_donor_to_acceptors = _make_splice_targets_dict(ext, 'donor', '+')
    pos_acceptor_to_donors = _make_splice_targets_dict(ext, 'acceptor', '+')
    neg_donor_to_acceptors = _make_splice_targets_dict(ext, 'donor', '-')
    neg_acceptor_to_donors = _make_splice_targets_dict(ext, 'acceptor', '-')

    annot = copy.deepcopy(annot)
    annot['upstream_donor_dist'] = np.nan
    annot['downstream_donor_dist'] = np.nan
    annot['upstream_acceptor_dist'] = np.nan
    annot['downstream_acceptor_dist'] = np.nan

    juncs = annot[annot.novel_donor & (annot.strand == '+')].index
    up, down = _dist_to_annot_donor_acceptor(annot.ix[juncs], 
                                             pos_acceptor_to_donors,
                                             '+', 
                                             'donor')
    annot.ix[juncs, 'upstream_donor_dist'] = up
    annot.ix[juncs, 'downstream_donor_dist'] = down

    juncs = annot[annot.novel_donor & (annot.strand == '-')].index
    up, down = _dist_to_annot_donor_acceptor(annot.ix[juncs], 
                                             neg_acceptor_to_donors,
                                             '-', 
                                             'donor')
    annot.ix[juncs, 'upstream_donor_dist'] = up
    annot.ix[juncs, 'downstream_donor_dist'] = down

    juncs = annot[annot.novel_acceptor & (annot.strand == '+')].index
    up, down = _dist_to_annot_donor_acceptor(annot.ix[juncs], 
                                             pos_donor_to_acceptors, 
                                             '+', 
                                             'acceptor')
    annot.ix[juncs, 'upstream_acceptor_dist'] = up
    annot.ix[juncs, 'downstream_acceptor_dist'] = down

    juncs = annot[annot.novel_acceptor & (annot.strand == '-')].index
    up, down = _dist_to_annot_donor_acceptor(annot.ix[juncs], 
                                             neg_donor_to_acceptors, 
                                             '-', 
                                             'acceptor')
    annot.ix[juncs, 'upstream_acceptor_dist'] = up
    annot.ix[juncs, 'downstream_acceptor_dist'] = down
    return annot

def _dist_to_annot_donor_acceptor(df, d, strand, novel_feature):
    """Find nearest annotated upstream/downstream donor/acceptor for novel
    donor/acceptors 

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with observed splice junctions (novel and known) with columns
        'chrom', 'first_bp_intron', 'last_bp_intron', 'strand', 'intron_motif',
        'annotated', 'ext_annotated', 'chrom:start', 'chrom:end', 'gene_id',
        'donor', 'acceptor', 'novel_donor', 'novel_acceptor'.

    d : dict
        If df contains novel donors, should be a dict whose keys are acceptors
        and whose values are the locations (integers) of all associated donors.
        If df contains novel acceptors, should be a dict whose keys are donors
        and whose values are the locations (integers) of all associated
        accepators.

    strand : str ('+' or '-')
        Strand that features are on.

    novel_feature : str ('donor' or 'acceptor')
        Whether the dataframe contains novel donors or novel acceptors.

    Returns
    -------
    up : list
        List of distances from novel feature to nearest feature (of same type)
        upstream. If upstream feature does not exist, the list will have a nan.

    down : list
        List of distances from novel feature to nearest feature (of same type)
        downstream. If upstream feature does not exist, the list will have a 
        nan.
    
    """
    if df.shape[0] > 0:
        assert len(set(df.strand)) == 1
    if novel_feature == 'donor':
        assert df.novel_donor.sum() == df.shape[0]
    if novel_feature == 'acceptor':
        assert df.novel_acceptor.sum() == df.shape[0]
    # For a novel donor, we want to return the distance to the nearest upstream
    # and downstream donors that use the same acceptor. For a novel acceptor, we
    # want to return the distance to the nearest upstream and downstream
    # acceptors that use the same donor. In some cases there may not be one of
    # the upstream or downstream donors/acceptors. In that case we will just
    # return nan.
    if strand == '+':
        if novel_feature == 'donor':
            annot_feature = 'acceptor'
            novel_location = 'start'
        if novel_feature == 'acceptor':
            annot_feature = 'donor'
            novel_location = 'end'
    if strand == '-':
        if novel_feature == 'donor':
            annot_feature = 'acceptor'
            novel_location = 'end' 
        if novel_feature == 'acceptor':
            annot_feature = 'donor'
            novel_location = 'start'

    upstream_dists = []
    downstream_dists = []
    
    for i in df.index:
        a = df.ix[i, annot_feature]
        diff = df.ix[i, novel_location] - d[a]
        pos = diff[diff > 0]
        neg = diff[diff < 0]
        if strand == '+':
            if pos.shape[0] == 0:
                upstream_dists.append(np.nan)
            else:
                upstream_dists.append(pos.min())
            if neg.shape[0] == 0:
                downstream_dists.append(np.nan)
            else:
                downstream_dists.append(np.abs(neg).min())
        if strand == '-':
            if pos.shape[0] == 0:
                downstream_dists.append(np.nan)
            else:
                downstream_dists.append(pos.min())
            if neg.shape[0] == 0:
                upstream_dists.append(np.nan)
            else:
                upstream_dists.append(np.abs(neg).min())
    return upstream_dists, downstream_dists
