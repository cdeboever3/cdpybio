import argparse
import pdb
import re

import pandas as pd

COLUMN_NAMES = ('chrom', 'first_bp_intron', 'last_bp_intron', 'strand',
                'intron_motif', 'annotated',
                'unique_junction_reads', 'multimap_junction_reads',
                'max_overhang')

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

def make_sjout_dict(fnL,define_sample_name=None):
    """Read multiple sjouts, return dict with keys as sample names and values
    as sjout dataframes

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
    sjoutD : dict
        Dict whose keys are sample names and values are sjout dataframes
    
    """
    if define_sample_name == None:
        define_sample_name = lambda x: x
    else:
        assert len(set([ define_sample_name(x) for x in fnL ])) == len(fnL)
    sjoutD = dict()
    for fn in fnL:
        sample = define_sample_name(fn)
        sjoutD[sample] = parse_sjout(fn)
    return sjoutD

def define_donor(row):
    if row['strand'] == '+':
        return '{}:{}:{}'.format(row['chrom'], row['first_bp_intron'], 
                                 row['strand'])
    if row['strand'] == '-':
        return '{}:{}:{}'.format(row['chrom'], row['last_bp_intron'], 
                                 row['strand'])

def define_acceptor(row):
    if row['strand'] == '+':
        return '{}:{}:{}'.format(row['chrom'], row['last_bp_intron'], 
                                 row['strand'])
    if row['strand'] == '-':
        return '{}:{}:{}'.format(row['chrom'], row['first_bp_intron'], 
                                 row['strand'])

def name_splice_site(L):
    return '{0}:{1}-{2}'.format(L[0],L[1],L[2])

def parse_sjout(fn):
    df = pd.read_table(fn, header=None, names=COLUMN_NAMES)
    df.index = [ '{0}:{1}-{2}'.format(df.ix[x,'chrom'],df.ix[x,'first_bp_intron'],df.ix[x,'last_bp_intron']) for x in df.index ]
    return df

def make_sjout_panel(sjoutD, total_jxn_cov_cutoff, statsN=None):
    """Filter junctions from many sjout files and make panel

    Parameters
    ----------
    sjoutD : dict
        Dict whose keys are sample names and values are sjout dataframes

    total_jxn_cov_cutoff : int
        If the unique read coverage of a junction summed over all samples is not
        greater than or equal to this value, the junction will not be included
        in the final output.

    statsN : filename str
        If provided, some stats will be printed to this file

    Returns
    -------
    sjoutP : pandas.Panel
        Panel where each dataframe(?) corresponds to an sjout file filtered to
        remove low coverage junctions.
    
    """
    # set of all junctions
    jxnS = reduce(lambda x,y: set(x) | set(y),
                  [ sjoutD[k].index for k in sjoutD.keys() ])

    jxn_keepS = set()
    for jxn in jxnS:
        if sum([ sjoutD[k].ix[jxn,'unique_junction_reads'] for k in sjoutD.keys() 
                 if jxn in sjoutD[k].index ]) >= total_jxn_cov_cutoff:
            jxn_keepS.add(jxn)

    for k in sjoutD.keys():
        sjoutD[k] = sjoutD[k].ix[jxn_keepS]

    sjoutP = pd.Panel(sjoutD)
    for col in ['unique_junction_reads', 'multimap_junction_reads',
                'max_splice']:
        sjoutP.ix[:,:,col] = sjoutP.ix[:,:,col].fillna(0)
    
    if statsN:
        statsF = open(statsN,'w')

        statsF.write('Number of junctions in sjout file per sample\n')
        for k in sjoutD.keys():
            statsF.write('{0}\t{1:,}\n'.format(k,sjoutD[k].shape[0]))
        statsF.write('\n')
  
        statsF.write('sjout panel size\t{0}\n\n'.format(sjoutP.shape))
        statsF.close()
    return sjoutP

def combine_star_sjout(fnL, total_jxn_cov_cutoff, statsN=None):
    """Combine multiple sjout files into a single file with the unique read
    coverage of each splice junction after filtering

    Parameters
    ----------
    fnL : list of strs of filenames or file handles
        List of filename of the SJ.out.tab files to read in

    total_jxn_cov_cutoff : int
        If the unique coverage of a junction summed over all samples is not
        greater than or equal to this value, the junction will not be included
        in the final output.

    define_sample_name : function that takes string as input
        Function mapping filename to sample name. For instance, you may have the
        sample name in the path and use a regex to extract it.  The sample names
        will be used as the column names. If this is not provided, the columns
        will be named as the input files.

    Returns
    -------
    sjoutD : dict
        Dict whose keys are sample names and values are sjout dataframes
    
    """
    sjoutD = make_sjout_dict(fnL,define_sample_name)
    sjoutP = make_sjout_panel(sjoutD, total_jxn_cov_cutoff)
   
def filter_jxns_donor_acceptor(sjoutP, gencode_jxnN, statsN=None):
gencode_jxnN, jxn_countN, 
jxn_annotN, gencode_infoN, 
    """Remove junctions that do not use a known donor or acceptor

    Parameters
    ----------
    fnL : list of strs of filenames or file handles
        List of filename of the SJ.out.tab files to read in

    total_jxn_cov_cutoff: int
        If the unique coverage of a junction summed over all samples is not
        greater than or equal to this value, the junction will not be included
        in the final output.

    define_sample_name: function that takes string as input
        Function mapping filename to sample name. For instance, you may have the
        sample name in the path and use a regex to extract it.  The sample names
        will be used as the column names. If this is not provided, the columns
        will be named as the input files.

    Returns
    -------
    sjoutD : dict
        Dict whose keys are sample names and values are sjout dataframes
    
    """
    # regular expression for parsing splice site definitions
    sjRE = re.compile('(.*:.*-.*):(\+|-)')
    # regular expression to get splice jxn without strand
    juncRE = re.compile('(.*):(\d*)-(\d*):') 
    
    # read gencode splice junction annotation 
    gencode_juncDF = pd.read_table(gencode_jxnN, header=None, 
                                   names=['junction','gene','chrom',
                                          'first_bp_intron','last_bp_intron',
                                          'strand'])
    statsF.write('Number of gencode junctions\t{0:,}\n\n'.format(gencode_juncDF.shape[0]))
  
    ### Filter data ###
    # we already filtered according to total_jxn_cov_cutoff, so we'll name the panel
    sjout_filteredP = sjoutP

    # make dataframe for junction annotation information. These are all of the jxn's from the star output after filtering (either in our splice jxn definitions or not) annotated with the information from the sjout file. I'm removing strand because I don't understand how star is assigning strand
    annotDF = reduce(pd.DataFrame.combine_first,[ sjout_filteredP.ix[item,:,['chrom','first_bp_intron','last_bp_intron','motif','annotated']].dropna() for item in sjout_filteredP.items ])

    # rename annotated column
    annotDF['star_annotated'] = annotDF.annotated == 1
    annotDF = annotDF.drop('annotated',axis=1)

    # make sure start and stop are ints
    annotDF.start = [ int(x) for x in annotDF.start]
    annotDF.end = [ int(x) for x in annotDF.end]
    
    # remove junction annotation information (for memory and speed)
    # sjout_filteredP = sjout_filteredP.drop(['chrom','first_bp_intron','last_bp_intron','strand','motif','annotated'],axis=2)
    
    statsF.write('Number of splice junctions after coverage filtering: {0:,}\n'.format(sjout_filteredP.shape[1]))
    statsF.write('Number STAR annotated: {0:,}\n\n'.format(annotDF['star_annotated'].sum()))
    
    ### make combined sjout file ###
    # add column for junction without strand. This is how the STAR output is currently indexed
    gencode_juncDF['junction_no_strand'] = gencode_juncDF.junction.apply(lambda x: juncRE.match(x).group().strip(':'))

    # find unique gencode junctions and keep them. Remove junctions used in more than one gene
    junctions_to_keepSE = gencode_juncDF.junction_no_strand.value_counts() == 1
    # drop duplicate junctions, but this leaves one copy of the duplicate junction
    uniq_gencode_juncDF = gencode_juncDF.drop_duplicates(cols='junction_no_strand')
    # so we will reindex and keep only the junctions we want
    uniq_gencode_juncDF.index = uniq_gencode_juncDF.junction_no_strand
    uniq_gencode_juncDF = uniq_gencode_juncDF.ix[junctions_to_keepSE]
    # reindex with strand info
    uniq_gencode_juncDF.index = uniq_gencode_juncDF.junction
    uniq_gencode_juncDF = uniq_gencode_juncDF.drop('junction',axis=1)

    # print number of unique gencode junctions
    statsF.write('Number of gencode junctions used only in one gene\t{0:,}\n\n'.format(uniq_gencode_juncDF.shape[0]))

    # make column showing whether junction is in gencode
    annotDF['gencode_annotated'] = False
    annotDF.ix[set(annotDF.index) & set(uniq_gencode_juncDF.junction_no_strand),'gencode_annotated'] = True

    # print number of junctions in gencode
    statsF.write('Number of observed junctions in gencode\t{0:,}\n'.format(sum(annotDF.gencode_annotated)))
    statsF.write('Number of observed junctions not in gencode\t{0:,}\n'.format(annotDF.shape[0] - sum(annotDF.gencode_annotated)))
    statsF.write('Number of observed junctions not in gencode but in STAR sj db\t{0:,}\n'.format(sum(annotDF.ix[annotDF.gencode_annotated == False,'star_annotated'])))
    statsF.write('Number of observed junctions not in gencode and not in STAR sj db\t{0:,}\n\n'.format(sum(annotDF.ix[annotDF.gencode_annotated == False,'star_annotated'].values == 0)))
  
    # add strand information to annotation of STAR junctions that are in gencode
    strandSE = pd.Series(uniq_gencode_juncDF.strand.values,index=uniq_gencode_juncDF.junction_no_strand)
    strandSE = strandSE[set(strandSE.index) & set(annotDF.index)]
    annotDF['strand'] = '*'
    annotDF.ix[strandSE.index,'strand'] = strandSE.values

    # add column for start and end location (chromosome plus position for uniqueness)
    annotDF['chr:start'] = annotDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['first_bp_intron']),axis=1)
    annotDF['chr:end'] = annotDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['last_bp_intron']),axis=1)

    # make sets of gencode starts and ends
    uniq_gencode_juncDF['chr:start'] = uniq_gencode_juncDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['first_bp_intron']),axis=1)
    uniq_gencode_juncDF['chr:end'] = uniq_gencode_juncDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['last_bp_intron']),axis=1)
    gencode_startS = set(uniq_gencode_juncDF['chr:start'].values)
    gencode_endS = set(uniq_gencode_juncDF['chr:end'].values)

    # remove junctions that don't have a start or end shared with gencode
    junctions_to_removeSE = annotDF[annotDF.gencode_annotated == False].apply(
            lambda x: (x['chr:start'] in gencode_startS) + (x['chr:end'] in gencode_endS) == 0,axis=1)
    annotDF = annotDF.drop(junctions_to_removeSE[junctions_to_removeSE].index)

    # print number of junctions remaining
    statsF.write('Number of junctions that share start or end with gencode junction\t{0:,}\n\n'.format(annotDF.shape[0]))

    # add column indicating which gene the junctions belong to for gencode jxn's
    geneSE = pd.Series(dict(zip(uniq_gencode_juncDF.junction_no_strand.values,uniq_gencode_juncDF.gene)))
    annotDF['gene_id'] = ''
    annotDF['gene_id'] = geneSE[annotDF.index]

    # now we'll figure out the genes for the non-gencode jxn's
    # map starts and ends to genes
    start_geneSE = pd.Series(dict(zip(uniq_gencode_juncDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['first_bp_intron']),axis=1),uniq_gencode_juncDF.gene)))
    end_geneSE = pd.Series(dict(zip(uniq_gencode_juncDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['last_bp_intron']),axis=1),uniq_gencode_juncDF.gene)))

    for ind in annotDF[annotDF.gencode_annotated == False].index:
        cur_start = annotDF.ix[ind,'chr:start'] 
        if cur_start in start_geneSE.index:
            annotDF.ix[ind,'gene_id'] = start_geneSE[cur_start]
    for ind in annotDF[annotDF.gencode_annotated == False].index:
        cur_end = annotDF.ix[ind,'chr:end'] 
        if cur_end in end_geneSE.index:
            annotDF.ix[ind,'gene_id'] = end_geneSE[cur_end]

    # now that we have the genes, we can assign strands to all junctions
    strandSE = pd.Series(dict(zip(uniq_gencode_juncDF.gene,uniq_gencode_juncDF.strand)))
    for ind in annotDF[annotDF.gencode_annotated == False].index:
        annotDF.ix[ind,'strand'] = strandSE[annotDF.ix[ind,'gene_id']]

    # and re-index with the strand info
    annotDF.index = [ x + ':' + annotDF.ix[x,'strand'] for x in annotDF.index ]

    # now we'll add donor and acceptor info
    annotDF['donor'] = annotDF.apply(lambda x: define_donor(x),axis=1)
    annotDF['acceptor'] = annotDF.apply(lambda x: define_acceptor(x),axis=1)

    # and whether the donor or acceptor is novel
    uniq_gencode_juncDF['donor'] = uniq_gencode_juncDF.apply(lambda x: define_donor(x),axis=1)
    uniq_gencode_juncDF['acceptor'] = uniq_gencode_juncDF.apply(lambda x: define_acceptor (x),axis=1)
    gencode_donorS = set(uniq_gencode_juncDF.donor)
    gencode_acceptorS = set(uniq_gencode_juncDF.acceptor)
    annotDF['novel_donor'] = False
    annotDF['novel_acceptor'] = False
    for ind in annotDF[annotDF.gencode_annotated == False].index:
        annotDF.ix[ind,'novel_donor'] = annotDF.ix[ind,'donor'] not in gencode_donorS
        annotDF.ix[ind,'novel_acceptor'] = annotDF.ix[ind,'acceptor'] not in gencode_acceptorS

    # print novel donor and acceptor info
    statsF.write('Number of novel donors\t{0:,}\n'.format(len(set(annotDF[annotDF.novel_donor].donor))))
    statsF.write('Number of novel junctions with novel donors\t{0:,}\n'.format(sum(annotDF.novel_donor)))
    statsF.write('Number of novel acceptors\t{0:,}\n'.format(len(set(annotDF[annotDF.novel_acceptor].acceptor))))
    statsF.write('Number of novel junctions with novel acceptors\t{0:,}\n'.format(sum(annotDF.novel_acceptor)))
    statsF.write('Number of novel junctions with gencode donor and acceptor\t{0:,}\n'.format(annotDF[annotDF.gencode_annotated].shape[0] - sum(annotDF.novel_donor) - sum(annotDF.novel_acceptor)))

    # sort by gene ID and start/end
    annotDF = annotDF.sort(columns=['gene_id','first_bp_intron','last_bp_intron'])

    annotDF.to_csv(jxn_annotN,sep='\t')
    uniq_gencode_juncDF.to_csv(gencode_infoN,sep='\t')
    statsF.close()

    # make file with counts for the junctions we are interested in
    countDF = sjout_filteredP.ix[:,[ juncRE.match(x).group().strip(':') for x in annotDF.index ],'unique_junction_reads']
    countDF.index = annotDF.index
    countDF.to_csv(jxn_countN,sep='\t')

def main():
    ### magic variables ###
    total_jxn_cov_cutoff= 20
    gencode_jxnN        = '/raid/databases/hg19/gene_annotations/gencode_v14/splice_junctions.tsv'
    jxn_countN          = 'junction_counts.tsv'
    jxn_annotN          = 'junction_info.tsv'
    gencode_infoN       = 'uniq_gencode_info.tsv'
    statsN              = 'combined_sjout_stats.txt'

    ### gather arguments from command line ###
    parser = argparse.ArgumentParser(
        description='''This script takes a list sjout files fom STAR alignments
        and combines them into a single python data structure after some
        filtering.''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('sjout_files', nargs='+', help='STAR sjout files.')
    parser.add_argument('-c', metavar='coverage_cutoff', type=int,
                        default=total_jxn_cov_cutoff, help='''If a junction is
                        covered by less than this many unique reads summed over
                        all samples, it will not be output in the counts
                        file.''')
    parser.add_argument('-g', metavar='gencode_junctions', default=gencode_jxnN,
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
   
    temp_fnL            = args.sjout_files
    total_jxn_cov_cutoff= args.c
    gencode_jxnN        = args.g
    jxn_countN          = args.jc
    jxn_annotN          = args.jg
    gencode_infoN       = args.gi
    statsN              = args.f
    debug               = args.debug

    ### start main ###
    combine_star_sjout(fnL,total_jxn_cov_cutoff,gencode_jxnN,jxn_countN,jxn_annotN,gencode_infoN,statsN)

if __name__ == '__main__':
    main()
