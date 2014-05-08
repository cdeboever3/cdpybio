import argparse
import pdb
import re

import pandas as pd

ANNOTATION_COLS = ['chrom', 'first_bp_intron', 'last_bp_intron', 'intron_motif',
                   'annotated']
COLUMN_NAMES = ('chrom', 'first_bp_intron', 'last_bp_intron', 'strand',
                'intron_motif', 'annotated',
                'unique_junction_reads', 'multimap_junction_reads',
                'max_overhang')

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

def make_sj_out_dict(fnL, define_sample_name=None):
    """Read multiple sj_outs, return dict with keys as sample names and values
    as sj_out dataframes

    Parameters
    ----------
    fnL : list of strs of filenames or file handles
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
        assert len(set([ define_sample_name(x) for x in fnL ])) == len(fnL)
    sj_outD = dict()
    for fn in fnL:
        sample = define_sample_name(fn)
        df = read_sj_out_tab(fn)
        index = df.apply(lambda x: _sj_out_junction(x), axis=1)
        assert len(index) == len(set(index))
        df.index = index
        sj_outD[sample] = df
    return sj_outD

def make_sj_out_panel(sj_outD, total_jxn_cov_cutoff=20, statsN=None):
    """Filter junctions from many sj_out files and make panel

    Parameters
    ----------
    sj_outD : dict
        Dict whose keys are sample names and values are sj_out dataframes

    total_jxn_cov_cutoff : int
        If the unique read coverage of a junction summed over all samples is not
        greater than or equal to this value, the junction will not be included
        in the final output.

    statsN : filename str
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
    for col in annot_cols:
        sj_outP.ix[:,annotDF.index,col] = [ annotDF[col] for i in
                                               range(sj_outP.shape[0]) ]
    
    # All of the splice junction annotation information is duplicated in each
    # dataframe in the panel, so we'll make a single dataframe holding that
    # information.
    annotDF = sj_outP.ix[0,:,ANNOTATION_COLS]

    if statsN:
        statsF = open(statsN,'w')

        statsF.write('Number of junctions in sj_out file per sample\n')
        for k in sj_outD.keys():
            statsF.write('{0}\t{1:,}\n'.format(k,sj_outD[k].shape[0]))
        statsF.write('\n')
  
        statsF.write('sj_out panel size\t{0}\n\n'.format(sj_outP.shape))
        statsF.close()
    return sj_outP,annotDF

# def combine_star_sj_out(fnL, total_jxn_cov_cutoff, statsN=None):
#     """Combine multiple sj_out files into a single file with the unique read
#     coverage of each splice junction after filtering
# 
#     Parameters
#     ----------
#     fnL : list of strs of filenames or file handles
#         List of filename of the SJ.out.tab files to read in
# 
#     total_jxn_cov_cutoff : int
#         If the unique coverage of a junction summed over all samples is not
#         greater than or equal to this value, the junction will not be included
#         in the final output.
# 
#     define_sample_name : function that takes string as input
#         Function mapping filename to sample name. For instance, you may have the
#         sample name in the path and use a regex to extract it.  The sample names
#         will be used as the column names. If this is not provided, the columns
#         will be named as the input files.
# 
#     Returns
#     -------
#     sj_outD : dict
#         Dict whose keys are sample names and values are sj_out dataframes
#     
#     """
#     sj_outD = make_sj_out_dict(fnL,define_sample_name)
#     sj_outP = make_sj_out_panel(sj_outD, total_jxn_cov_cutoff)
  
def read_external_annotation(fn):
    """Read file with junctions from some database. This does not have to be the
    same splice junction database used with STAR

    Parameters
    ----------
    fn: filename str
        File with splice junction from annotation.

    Returns
    -------
    juncDF : pandas.DataFrame
        DataFrame indexed by splice junction
    
    """
    juncDF = pd.read_table(
        fn, header=None, names=['junction','gene','chrom',
                                  'first_bp_intron','last_bp_intron', 'strand'])
    
    juncDF['junction_no_strand'] = juncDF.junction.apply(
        lambda x: juncRE.match(x).group().strip(':'))

    # In rare cases, a splice junction might be used by more than one gene. For
    # my purposes, these cases are confounding, so I will remove all such splice
    # junctions. 
    junctions_to_keepSE = juncDF.junction_no_strand.value_counts() == 1
    # Drop duplicate junctions, but this leaves one copy of the duplicate
    # junction, so we will reindex and keep only the junctions we want.
    juncDF = juncDF.drop_duplicates(cols='junction_no_strand')
    juncDF.index = juncDF.junction_no_strand
    juncDF = juncDF.ix[junctions_to_keepSE]
    # Reindex with strand info.
    juncDF.index = juncDF.junction
    juncDF = juncDF.drop('junction',axis=1)
    return juncDF

def filter_jxns_donor_acceptor(sj_outP, jxnN, statsN=None):
# jxnN, jxn_countN, 
# jxn_annotN, gencode_infoN, 
    """Remove junctions that do not use an annotated donor or acceptor according
    to external junction annotation. Add strand and gene information for
    junctions according to external annoation (STAR strand ignored).

    Parameters
    ----------
    sj_outP : pandas.Panel
        Panel where each dataframe(?) corresponds to an sj_out file filtered to
        remove low coverage junctions.

    jxnN: filename str
        File defining annotated splice sites. These can differ from those
        provided to STAR for alignment. TODO: describe file format briefly

    Returns
    -------
    sj_outD : dict
        Dict whose keys are sample names and values are sj_out dataframes
    
    """
    sjRE = re.compile('(.*:.*-.*):(\+|-)')
    juncRE = re.compile('(.*):(\d*)-(\d*):') 
    
    juncDF = read_external_annotation(fn)
   
    # All of the splice junction annotation information is duplicated in each
    # dataframe in the panel, so we'll make a single dataframe holding that
    # information.
    annotDF = sj_outP.ix[0,:,ANNOTATION_COLS]

    # Add column showing whether junction is in external annotation.
    annotDF['ext_annotated'] = False
    annotDF.ix[set(annotDF.index) & set(juncDF.junction_no_strand),
               'ext_annotated'] = True

    # add strand information to annotation of STAR junctions that are in gencode
    strandSE = pd.Series(juncDF.strand.values,index=juncDF.junction_no_strand)
    strandSE = strandSE[set(strandSE.index) & set(annotDF.index)]
    annotDF['strand'] = '*'
    annotDF.ix[strandSE.index,'strand'] = strandSE.values

    # add column for start and end location (chromosome plus position for uniqueness)
    annotDF['chr:start'] = annotDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['first_bp_intron']),axis=1)
    annotDF['chr:end'] = annotDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['last_bp_intron']),axis=1)

    # make sets of gencode starts and ends
    juncDF['chr:start'] = juncDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['first_bp_intron']),axis=1)
    juncDF['chr:end'] = juncDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['last_bp_intron']),axis=1)
    gencode_startS = set(juncDF['chr:start'].values)
    gencode_endS = set(juncDF['chr:end'].values)

    # remove junctions that don't have a start or end shared with gencode
    junctions_to_removeSE = annotDF[annotDF.ext_annotated == False].apply(
            lambda x: (x['chr:start'] in gencode_startS) + (x['chr:end'] in gencode_endS) == 0,axis=1)
    annotDF = annotDF.drop(junctions_to_removeSE[junctions_to_removeSE].index)

    # print number of junctions remaining
    statsF.write('Number of junctions that share start or end with gencode junction\t{0:,}\n\n'.format(annotDF.shape[0]))

    # add column indicating which gene the junctions belong to for gencode jxn's
    geneSE = pd.Series(dict(zip(juncDF.junction_no_strand.values,juncDF.gene)))
    annotDF['gene_id'] = ''
    annotDF['gene_id'] = geneSE[annotDF.index]

    # now we'll figure out the genes for the non-gencode jxn's
    # map starts and ends to genes
    start_geneSE = pd.Series(dict(zip(juncDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['first_bp_intron']),axis=1),juncDF.gene)))
    end_geneSE = pd.Series(dict(zip(juncDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['last_bp_intron']),axis=1),juncDF.gene)))

    for ind in annotDF[annotDF.ext_annotated == False].index:
        cur_start = annotDF.ix[ind,'chr:start'] 
        if cur_start in start_geneSE.index:
            annotDF.ix[ind,'gene_id'] = start_geneSE[cur_start]
    for ind in annotDF[annotDF.ext_annotated == False].index:
        cur_end = annotDF.ix[ind,'chr:end'] 
        if cur_end in end_geneSE.index:
            annotDF.ix[ind,'gene_id'] = end_geneSE[cur_end]

    # now that we have the genes, we can assign strands to all junctions
    strandSE = pd.Series(dict(zip(juncDF.gene,juncDF.strand)))
    for ind in annotDF[annotDF.ext_annotated == False].index:
        annotDF.ix[ind,'strand'] = strandSE[annotDF.ix[ind,'gene_id']]

    # and re-index with the strand info
    annotDF.index = [ x + ':' + annotDF.ix[x,'strand'] for x in annotDF.index ]

    # now we'll add donor and acceptor info
    annotDF['donor'] = annotDF.apply(lambda x: _sj_out_donor(x),axis=1)
    annotDF['acceptor'] = annotDF.apply(lambda x: _sj_out_acceptor(x),axis=1)

    # and whether the donor or acceptor is novel
    juncDF['donor'] = juncDF.apply(lambda x: _sj_out_donor(x),axis=1)
    juncDF['acceptor'] = juncDF.apply(lambda x: _sj_out_acceptor (x),axis=1)
    gencode_donorS = set(juncDF.donor)
    gencode_acceptorS = set(juncDF.acceptor)
    annotDF['novel_donor'] = False
    annotDF['novel_acceptor'] = False
    for ind in annotDF[annotDF.ext_annotated == False].index:
        annotDF.ix[ind,'novel_donor'] = annotDF.ix[ind,'donor'] not in gencode_donorS
        annotDF.ix[ind,'novel_acceptor'] = annotDF.ix[ind,'acceptor'] not in gencode_acceptorS

    # print novel donor and acceptor info
    statsF.write('Number of novel donors\t{0:,}\n'.format(len(set(annotDF[annotDF.novel_donor].donor))))
    statsF.write('Number of novel junctions with novel donors\t{0:,}\n'.format(sum(annotDF.novel_donor)))
    statsF.write('Number of novel acceptors\t{0:,}\n'.format(len(set(annotDF[annotDF.novel_acceptor].acceptor))))
    statsF.write('Number of novel junctions with novel acceptors\t{0:,}\n'.format(sum(annotDF.novel_acceptor)))
    statsF.write('Number of novel junctions with gencode donor and acceptor\t{0:,}\n'.format(annotDF[annotDF.ext_annotated].shape[0] - sum(annotDF.novel_donor) - sum(annotDF.novel_acceptor)))

    # sort by gene ID and start/end
    annotDF = annotDF.sort(columns=['gene_id','first_bp_intron','last_bp_intron'])

    annotDF.to_csv(jxn_annotN,sep='\t')
    juncDF.to_csv(gencode_infoN,sep='\t')

    # make file with counts for the junctions we are interested in
    countDF = sj_out_filteredP.ix[:,[ juncRE.match(x).group().strip(':') for x in annotDF.index ],'unique_junction_reads']
    countDF.index = annotDF.index
    countDF.to_csv(jxn_countN,sep='\t')

    # TODO: new stats location, everything above should come down
    statsF = open(statsN,'w')
    statsF.write('Number of annotated junctions\t{0:,}\n\n'.format(juncDF.shape[0]))
    statsF.write('Number STAR annotated: {0:,}\n\n'.format(annotDF['annotated'].sum()))

    # print number of unique gencode junctions
    statsF.write('Number of gencode junctions used only in one gene\t{0:,}\n\n'.format(juncDF.shape[0]))

    # print number of junctions in gencode
    statsF.write('Number of observed junctions in gencode\t{0:,}\n'.format(sum(annotDF.ext_annotated)))
    statsF.write('Number of observed junctions not in gencode\t{0:,}\n'.format(annotDF.shape[0] - sum(annotDF.ext_annotated)))
    statsF.write('Number of observed junctions not in gencode but in STAR sj db\t{0:,}\n'.format(sum(annotDF.ix[annotDF.ext_annotated == False,'annotated'])))
    statsF.write('Number of observed junctions not in gencode and not in STAR sj db\t{0:,}\n\n'.format(sum(annotDF.ix[annotDF.ext_annotated == False,'annotated'].values == 0)))
  
    statsF.close()

def main():
    ### magic variables ###
    total_jxn_cov_cutoff= 20
    jxnN        = '/raid/databases/hg19/gene_annotations/gencode_v14/splice_junctions.tsv'
    jxn_countN          = 'junction_counts.tsv'
    jxn_annotN          = 'junction_info.tsv'
    gencode_infoN       = 'uniq_gencode_info.tsv'
    statsN              = 'combined_sj_out_stats.txt'

    ### gather arguments from command line ###
    parser = argparse.ArgumentParser(
        description='''This script takes a list sj_out files fom STAR alignments
        and combines them into a single python data structure after some
        filtering.''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('sj_out_files', nargs='+', help='STAR sj_out files.')
    parser.add_argument('-c', metavar='coverage_cutoff', type=int,
                        default=total_jxn_cov_cutoff, help='''If a junction is
                        covered by less than this many unique reads summed over
                        all samples, it will not be output in the counts
                        file.''')
    parser.add_argument('-g', metavar='gencode_junctions', default=jxnN,
                        help='''Tsv file describing gencode splice 
                        junctions.''')
    parser.add_argument('-jc', metavar='jxn_counts', default=jxn_countN,
                        help='''File name for splice junction 
                        counts.''')
    parser.add_argument('-jg', metavar='jxn_genes', default=jxn_annotN,
                        help='''File name for splice junction annotations. This
                        file is especially useful for novel junctions and adds
                        strand information.''')
    parser.add_argument('-gi', metavar='gencode_info', default=gencode_infoN,
                        help='''Output file for gencode splice junction info
                        filtered to include only unique junctions. Also includes
                        several extra columns''')
    parser.add_argument('-f', metavar='stats_file', default=statsN, help='''File
                        to print some statistics to.''')
    parser.add_argument('--debug', action='store_true', help='''Enable python
                        debugger.''')
    
    args = parser.parse_args()
   
    temp_fnL            = args.sj_out_files
    total_jxn_cov_cutoff= args.c
    jxnN        = args.g
    jxn_countN          = args.jc
    jxn_annotN          = args.jg
    gencode_infoN       = args.gi
    statsN              = args.f
    debug               = args.debug

    ### start main ###
    combine_star_sj_out(fnL,total_jxn_cov_cutoff,jxnN,jxn_countN,jxn_annotN,gencode_infoN,statsN)

if __name__ == '__main__':
    main()
