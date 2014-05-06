import argparse
import copy
import glob
import os
import pdb
import sys

import pandas as pd

def combine_express_output(fnL,
                           column,
                           namesL=None,
                           transcript_to_geneN=None,
                           define_basename=None,
                           debug=False):
    """
    Combine eXpress output files

    Parameters:
    -----------

    fnL: list of strings
        List of paths to results.xprs files.

    column: string
        Column name of eXpress output to combine.

    namesL: list of strings
        Names to use for columns of output files. Overrides define_basename if
        provided.

    transcript_to_geneN: string
        File with transcript-to-gene mapping. Transcripts should be in first
        column and genes in second column. 

    define_basename: function that takes string as input
        Function mapping filename to sample name (or basename). For instance,
        you may have the basename in the path and use a regex to extract it.
        The basenames will be used as the column names. If this is not provided,
        the columns will be named as the input files.

    debug: boolean
        Passing True will trigger any debugging statements.

    """
    if namesL is not None:
        assert len(namesL) == len(fnL)
    if define_basename is None:
        define_basename = lambda x: x
    
    transcriptL = []
    for i,fn in enumerate(fnL):
        if namesL is not None:
            bn = namesL[i]
        else:
            bn = define_basename(fn)
        tDF = pd.read_table(fn,index_col=1,header=0)
        se = tDF[column]
        se.name = bn
        transcriptL.append(se)
    transcriptDF = pd.DataFrame(transcriptL).T
    transcriptDF.index.name = 'transcript'

    if transcript_to_geneN is not None:
        tgDF = pd.read_table(transcript_to_geneN,
                             index_col=0,
                             header=None,
                             names=['gene_id'])
        geneDF = copy.deepcopy(transcriptDF)
        geneDF['gene'] = tgDF.ix[geneDF.index]
        geneDF = geneDF.groupby('gene').sum()

    if geneDF is not None:
        return transcriptDF,geneDF
    else:
        return transcriptDF,None

def main():
    column = 'eff_counts'
    transcript_outN = 'express_transcript_[column].tsv'
    gene_outN = 'express_gene_[column].tsv'

    parser = argparse.ArgumentParser(
        description='''This script combines the eXpress output files for 
        multiple samples into a single file. If you provide a file defining 
        the genes for each transcript, the values of the transcripts 
        for that gene will be summed together to provide gene-level 
        estimates.''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('files', nargs='+', help='''eXpress output files 
                        (results.xprs). Wildcards may be used.''')
    parser.add_argument('--column',
                        choices=['tot_counts','uniq_counts','est_counts',
                                 'eff_counts','fpkm'],
                        help='Column to combine from eXpress output.',
                        default=column)
    parser.add_argument('--names',
                        nargs='+',
                        help='Names to be used for columns of output. Must '
                        'match the number of input files and be in the same '
                        'order.')
    parser.add_argument('-t',
                        metavar='transcript output file',
                        help='Output file for transcripts.',
                        default=transcript_outN)
    parser.add_argument('-tg', metavar='transcript to gene mapping', 
                        help='''File with transcripts in first column and genes 
                        in second column. Used to aggregate 
                        transcript values into gene-level values. If not 
                        provided, only transcript results are returned.''')
    parser.add_argument('-g', metavar='gene output file', help='''Output file 
                        for genes if a conversion file is provided.''',
                        default='express_gene_[column].tsv')
    parser.add_argument('--debug', action='store_true', help='''Enable python 
                        debugger.''')

    args = parser.parse_args()

    temp_fnL = args.files
    namesL = args.names
    column = args.column
    transcript_outN = args.t
    transcript_to_geneN = args.tg
    gene_outN = args.g
    debug = args.debug

    if transcript_outN == 'express_transcript_[column].tsv':
        transcript_outN = 'express_transcript_{0}.tsv'.format(column)
    if gene_outN == 'express_gene_[column].tsv':
        gene_outN = 'express_gene_{0}.tsv'.format(column)

    fnL = []
    for fn in temp_fnL:
        fnL += glob.glob(fn)
    if namesL is not None:
        assert len(namesL) == len(fnL)

    transcriptDF,geneDF = combine_express_output(fnL,
                                                 column,
                                                 namesL,
                                                 transcript_to_geneN,
                                                 None,
                                                 debug)
    transcriptDF.to_csv(transcript_outN,sep='\t')
    if geneDF is not None:
        geneDF.to_csv(gene_outN,sep='\t')

if __name__ == '__main__':
    main()
