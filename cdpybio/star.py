import copy
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
        Panel where each dataframe corresponds to an sj_out file filtered to
        remove low coverage junctions. Each dataframe has COUNT_COLS =
        ('unique_junction_reads', 'multimap_junction_reads', 'max_overhang')

    annotDF : pandas.DataFrame
        Dataframe with values ANNOTATION_COLS = ('chrom', 'start', 
        'end', 'intron_motif', 'annotated') that are otherwise
        duplicated in the panel.
    
    """
    # Remove any junctions that don't have any uniquely mapped junction reads.
    # Even if a junction passes the cutoff in other samples, we are only
    # concerned with unique counts.
    # set of all junctions
    jxnS = reduce(lambda x,y: set(x) | set(y),
                  [ sj_outD[k].index for k in sj_outD.keys() ])

    jxn_keepS = set()
    sj_outD = copy.deepcopy(sj_outD)
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
    annotDF['start'] = annotDF['start'].astype(int)
    annotDF['end'] = annotDF['end'].astype(int)
    annotDF['annotated'] = annotDF['annotated'].astype(bool)

    sj_outP = sj_outP.ix[:,:,COUNT_COLS].astype(int)

    if statsfile:
        statsF = open(statsfile,'w')

        statsF.write('Number of junctions in sj_out file per sample\n')
        for k in sj_outD.keys():
            statsF.write('{0}\t{1:,}\n'.format(k,sj_outD[k].shape[0]))
        statsF.write('\n')
  
        statsF.write('sj_out panel size\t{0}\n\n'.format(sj_outP.shape))
        statsF.close()
    return sj_outP, annotDF

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
    extDF = pd.read_table(fn, index_col=0, header=0)
    
    # In rare cases, a splice junction might be used by more than one gene. For
    # my purposes, these cases are confounding, so I will remove all such splice
    # junctions. 
    intron_count = extDF.intron.value_counts()
    extDF['intron_count'] = extDF.intron.apply(lambda x: intron_count.ix[x])
    extDF = extDF[extDF.intron_count == 1]
    extDF = extDF.drop('intron_count', axis=1)

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
        'chrom', 'start', 'end', 'intron_motif',
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
    annotDF['chr:start'] = annotDF.apply(
        lambda x: '{}:{}'.format(x['chrom'],x['start']),axis=1)
    annotDF['chr:end'] = annotDF.apply(
        lambda x: '{}:{}'.format(x['chrom'],x['end']),axis=1)

    ext_startS = set(extDF['chr:start'].values)
    ext_endS = set(extDF['chr:end'].values)

    # Remove junctions that don't have a start or end shared with external
    # database.
    junctions_to_removeSE = annotDF[annotDF.ext_annotated == False].apply(
        lambda x: (x['chr:start'] in ext_startS) + 
        (x['chr:end'] in ext_endS) == 0,axis=1)
    if junctions_to_removeSE.shape[0] > 0:
        annotDF = annotDF.drop(junctions_to_removeSE[junctions_to_removeSE].index)

    # Add column indicating which gene the junctions belong to for annotated
    # jxn's.
    geneSE = pd.Series(dict(zip(extDF.intron.values,extDF.gene)))
    annotDF['gene_id'] = ''
    annotDF['gene_id'] = geneSE[annotDF.index]

    # Now we'll figure out the genes for the junctions that aren't in our
    # database. We can associate each start and end with a gene and use this.
    start_geneSE = pd.Series(dict(zip(extDF['chr:start'], extDF.gene)))
    end_geneSE = pd.Series(dict(zip(extDF['chr:end'], extDF.gene)))

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
    # extDF['donor'] = extDF.apply(lambda x: _sj_out_donor(x),axis=1)
    # extDF['acceptor'] = extDF.apply(lambda x: _sj_out_acceptor (x),axis=1)
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
    annotDF = annotDF.sort(columns=['gene_id', 'start', 
                                    'end'])

    # Make file with counts for the junctions we are interested in.
    L = [ juncRE.match(x).group().strip(':') for x in annotDF.index ]
    countDF = sj_outP.ix[:, L, 'unique_junction_reads']
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

def _make_splice_targets_dict(df, feature, strand):
    """Make dict mapping each donor to the location of all acceptors it splices
    to or each acceptor to all donors it splices from.

    Parameters
    ----------
    df : pandas.DataFrame 
        Dataframe with splice junction information from external database
        containing columns 'gene', 'chrom', 'start', 'end', 'strand',
        'chr:start', 'chr:end', 'donor', 'acceptor', 'intron'.

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
        'annotated', 'ext_annotated', 'chr:start', 'chr:end', 'gene_id',
        'donor', 'acceptor', 'novel_donor', 'novel_acceptor'.

    ext : pandas.DataFrame 
        Dataframe with splice junction information from external database
        containing columns 'gene', 'chrom', 'start', 'end', 'strand',
        'chr:start', 'chr:end', 'donor', 'acceptor', 'intron'.

    Returns
    -------
    annot : pandas.DataFrame
        Annotation information with added columns for distance from novel
        donor/acceptor to nearest upstream/downstream donor/acceptor.
    
    """
    pos_donor_to_acceptors = _make_splice_targets_dict(ext, '+', 'donor')
    pos_acceptor_to_donors = _make_splice_targets_dict(ext, '+', 'acceptor')
    neg_donor_to_acceptors = _make_splice_targets_dict(ext, '-', 'donor')
    neg_acceptor_to_donors = _make_splice_targets_dict(ext, '-', 'acceptor')

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
    annot.ix[juncs, 'upstream_donor_dist'] = up
    annot.ix[juncs, 'downstream_donor_dist'] = down

    juncs = annot[annot.novel_acceptor & (annot.strand == '-')].index
    up, down = _dist_to_annot_donor_acceptor(annot.ix[juncs], 
                                             neg_donor_to_acceptors, 
                                             '-', 
                                             'acceptor')
    annot.ix[juncs, 'upstream_donor_dist'] = up
    annot.ix[juncs, 'downstream_donor_dist'] = down
    return annot

def _dist_to_annot_donor_acceptor(df, d, strand, novel_feature):
    """Find nearest annotated upstream/downstream donor/acceptor for novel
    donor/acceptors 

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with observed splice junctions (novel and known) with columns
        'chrom', 'first_bp_intron', 'last_bp_intron', 'strand', 'intron_motif',
        'annotated', 'ext_annotated', 'chr:start', 'chr:end', 'gene_id',
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
