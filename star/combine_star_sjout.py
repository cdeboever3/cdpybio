import cPickle,re,sys,argparse,glob,pdb
import pandas as pd
import numpy as np
from project_tools import *

def define_donor(row):
    if row['strand'] == '+':
        return row['chrom'] + ':' + str(int(row['start'])) + ':' + row['strand']
    elif row['strand'] == '-':
        return row['chrom'] + ':' + str(int(row['end'])) + ':' + row['strand']

def define_acceptor(row):
    if row['strand'] == '+':
        return row['chrom'] + ':' + str(int(row['end'])) + ':' + row['strand']
    elif row['strand'] == '-':
        return row['chrom'] + ':' + str(int(row['start'])) + ':' + row['strand']

def name_splice_site(L):
    return '{0}:{1}-{2}'.format(L[0],L[1],L[2])

def parse_SJout(fn):
    # 1: chrom, 2: first base of intron (1-based), 3: last base of intron (1-based), 4: strand, 
    # 5: intron motif (0: non-canonical, 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT)
    # 6: (0 unannotated, 1 annotated)
    # 7: number of unique reads crossing junction
    # 8: number of multi reads crossing junction
    # 9: maximum spliced alignment overhang
    df = pd.read_table(fn,header=None,names=['chrom','start','end','strand','motif','annotated','num_uniq','num_mult','max_splice'])
    df.index = [ '{0}:{1}-{2}'.format(df.ix[x,'chrom'],df.ix[x,'start'],df.ix[x,'end']) for x in df.index ]
    return df

def combine_star_sjout(fnL,total_jxn_cov_cutoff,gencode_jxnN,jxn_countN,jxn_annotN,gencode_infoN,statsN):
    # open file for writing stats
    statsF = open(statsN,'w')

    # regular expression for parsing splice site definitions
    sjRE = re.compile('(.*:.*-.*):(\+|-)')
    # regular expression to get splice jxn without strand
    juncRE = re.compile('(.*):(\d*)-(\d*):') 
    
    ### Read Data and Make Data Structures ###
    # read gencode splice junction annotation 
    gencode_juncDF = pd.read_table(gencode_jxnN,header=None,names=['junction','gene','chrom','start','end','strand'])
    statsF.write('Number of gencode junctions\t{0:,}\n\n'.format(gencode_juncDF.shape[0]))
   
    # make dict to hold splice junction coverage info for each sample
    SJoutD = dict()
    for fn in fnL:
        sample = define_basename('default',fn,basenameL)
        SJoutD[sample] = parse_SJout(fn)

    # write the size of each SJout file to the stats file
    statsF.write('Number of junctions in SJout file per sample\n')
    for k in SJoutD.keys():
        statsF.write('{0}\t{1:,}\n'.format(k,SJoutD[k].shape[0]))
    statsF.write('\n')
  
    # filter junctions before making panel. This is necessary when we have many
    # samples
    # first, make set of all junctions
    jxnS = reduce(lambda x,y: set(x) | set(y),[ SJoutD[k].index for k in SJoutD.keys() ])

    # now, count how many reads span each junction and only keep those that pass filter
    jxn_keepS = set()
    for jxn in jxnS:
        if sum([ SJoutD[k].ix[jxn,'num_uniq'] for k in SJoutD.keys() 
                 if jxn in SJoutD[k].index ]) >= total_jxn_cov_cutoff:
            jxn_keepS.add(jxn)

    # keep only the junctions that pass the filter for each of the dataframes
    for k in SJoutD.keys():
        SJoutD[k] = SJoutD[k].ix[jxn_keepS]

    # make panel out of dict of dataframes
    SJoutP = pd.Panel(SJoutD)
   
    # Write the size of the SJout panel to the stats file
    statsF.write('SJout panel size\t{0}\n\n'.format(SJoutP.shape))
    
    # replace NaN's with zeros for counting columns
    SJoutP.ix[:,:,['num_uniq','num_mult','max_splice']] = SJoutP.ix[:,:,['num_uniq','num_mult','max_splice']].fillna(0)
    
    ### Filter data ###
    # we already filtered according to total_jxn_cov_cutoff, so we'll name the panel
    SJout_filteredP = SJoutP

    # make dataframe for junction annotation information. These are all of the jxn's from the star output after filtering (either in our splice jxn definitions or not) annotated with the information from the SJout file. I'm removing strand because I don't understand how star is assigning strand
    annotDF = reduce(pd.DataFrame.combine_first,[ SJout_filteredP.ix[item,:,['chrom','start','end','motif','annotated']].dropna() for item in SJout_filteredP.items ])

    # rename annotated column
    annotDF['star_annotated'] = annotDF.annotated == 1
    annotDF = annotDF.drop('annotated',axis=1)

    # make sure start and stop are ints
    annotDF.start = [ int(x) for x in annotDF.start]
    annotDF.end = [ int(x) for x in annotDF.end]
    
    # remove junction annotation information (for memory and speed)
    # SJout_filteredP = SJout_filteredP.drop(['chrom','start','end','strand','motif','annotated'],axis=2)
    
    statsF.write('Number of splice junctions after coverage filtering: {0:,}\n'.format(SJout_filteredP.shape[1]))
    statsF.write('Number STAR annotated: {0:,}\n\n'.format(annotDF['star_annotated'].sum()))
    
    ### make combined SJout file ###
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
    annotDF['chr:start'] = annotDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['start']),axis=1)
    annotDF['chr:end'] = annotDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['end']),axis=1)

    # make sets of gencode starts and ends
    uniq_gencode_juncDF['chr:start'] = uniq_gencode_juncDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['start']),axis=1)
    uniq_gencode_juncDF['chr:end'] = uniq_gencode_juncDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['end']),axis=1)
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
    start_geneSE = pd.Series(dict(zip(uniq_gencode_juncDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['start']),axis=1),uniq_gencode_juncDF.gene)))
    end_geneSE = pd.Series(dict(zip(uniq_gencode_juncDF.apply(lambda x: '{0}:{1}'.format(x['chrom'],x['end']),axis=1),uniq_gencode_juncDF.gene)))

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
    annotDF = annotDF.sort(columns=['gene_id','start','end'])

    annotDF.to_csv(jxn_annotN,sep='\t')
    uniq_gencode_juncDF.to_csv(gencode_infoN,sep='\t')
    statsF.close()

    # make file with counts for the junctions we are interested in
    countDF = SJout_filteredP.ix[:,[ juncRE.match(x).group().strip(':') for x in annotDF.index ],'num_uniq']
    countDF.index = annotDF.index
    countDF.to_csv(jxn_countN,sep='\t')

def main():
    ### magic variables ###
    total_jxn_cov_cutoff= 20
    gencode_jxnN        = '/raid/databases/hg19/gene_annotations/gencode_v14/splice_junctions.tsv'
    jxn_countN          = 'junction_counts.tsv'
    jxn_annotN          = 'junction_info.tsv'
    gencode_infoN       = 'uniq_gencode_info.tsv'
    statsN              = 'combined_SJout_stats.txt'

    ### gather arguments from command line ###
    parser = argparse.ArgumentParser(description='This script takes a list SJout files fom STAR alignments and combines them into a single python data structure after some filtering.')
    parser.add_argument('SJout_files',nargs='+',help='STAR SJout files.')
    parser.add_argument('-c',metavar='coverage_cutoff',type=int,default=total_jxn_cov_cutoff,help='If a junction is covered by less than this many reads in all samples, it will not be output in the counts file.')
    parser.add_argument('-g',metavar='gencode_junctions',default=gencode_jxnN,help='Tsv file describing gencode splice junctions. Default: {0}'.format(gencode_jxnN))
    parser.add_argument('-jc',metavar='jxn_counts',default=jxn_countN,help='File name for splice junction counts.'.format(jxn_countN))
    parser.add_argument('-jg',metavar='jxn_genes',default=jxn_annotN,help='File name for splice junction annotations. This file is especially useful for novel junctions and adds strand information.'.format(jxn_annotN))
    parser.add_argument('-gi',metavar='gencode_info',default=gencode_infoN,help='Output file for gencode splice junction info filtered to include only unique junctions. Also includes several extra columns'.format(gencode_infoN))
    parser.add_argument('-f',metavar='stats_file',default=statsN,help='File to print some statistics to. Default: {0}'.format(statsN))
    parser.add_argument('--debug', action='store_true', help='Enable python debugger.')
    
    args = parser.parse_args()
   
    temp_fnL            = args.SJout_files
    total_jxn_cov_cutoff= args.c
    gencode_jxnN        = args.g
    jxn_countN          = args.jc
    jxn_annotN          = args.jg
    gencode_infoN       = args.gi
    statsN              = args.f
    debug               = args.debug

    ### start main ###
    fnL = []
    for fn in temp_fnL:
        fnL += glob.glob(fn)
    
    combine_star_sjout(fnL,total_jxn_cov_cutoff,gencode_jxnN,jxn_countN,jxn_annotN,gencode_infoN,statsN)

if __name__ == '__main__':
    main()
