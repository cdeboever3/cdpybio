from copy import deepcopy
import pandas as pd
from pandas.util.testing import assert_frame_equal
from pandas.util.testing import assert_series_equal
import pytest

import cdpybio as cpb

class TestMakeTranscriptGeneSe:
    def test_two_genes(self):
        transcripts = [u'ENST00000456328.2', u'ENST00000515242.2',
                       u'ENST00000518655.2', u'ENST00000450305.2',
                       u'ENST00000466430.1', u'ENST00000477740.1',
                       u'ENST00000471248.1', u'ENST00000453576.2']
        genes = ['ENSG00000223972.4', 'ENSG00000223972.4', 'ENSG00000223972.4',
                 'ENSG00000223972.4', 'ENSG00000238009.2', 'ENSG00000238009.2',
                 'ENSG00000238009.2', 'ENSG00000238009.2']
        se = pd.Series(genes, index=transcripts)
        se2 = cpb.gencode.make_transcript_gene_se('annot.gtf')
        assert_series_equal(se, se2)

class TestMakeGeneInfoDf:
    def test_make_gene_info_df(self):

        df = pd.DataFrame([['pseudogene', 'KNOWN', 'DDX11L1'],
                           ['antisense', 'NOVEL', 'RP11-34P13.7']], 
                          index=[u'ENSG00000223972.4', u'ENSG00000238009.2'],
                          columns=[u'gene_type', u'gene_status', u'gene_name'])
        df.index.name = 'gene_id'
        df2 = cpb.gencode.make_gene_info_df('annot.gtf')
        assert_frame_equal(df, df2)

class TestMakeSpliceJunctionDf:
    def test_pos(self):
        df = pd.DataFrame(index=[u'chr1:12228-12612:+', u'chr1:12722-13220:+',
                                 u'chr1:12722-13224:+', u'chr1:12228-12594:+',
                                 u'chr1:12722-13402:+', u'chr1:13656-13660:+',
                                 u'chr1:12058-12178:+', u'chr1:12698-12974:+',
                                 u'chr1:13053-13220:+', u'chr1:13375-13452:+'])

        df['gene'] = 'ENSG00000223972.4'
        df['chrom'] = 'chr1'
        df['start'] = [12228, 12722, 12722, 12228, 12722, 13656, 12058, 12698,
                       13053, 13375]
        df['end'] = [12612, 13220, 13224, 12594, 13402, 13660, 12178, 12974,
                     13220, 13452]
        df['strand'] = '+'
        df['chrom:start'] = df.chrom + ':' + df.start.astype(str)
        df['chrom:end'] = df.chrom + ':' + df.end.astype(str)
        df['donor'] = df['chrom:start'] + ':' + df.strand
        df['acceptor'] = df['chrom:end'] + ':' + df.strand
        df['intron'] = (df.chrom + ':' + df.start.astype(str) + '-' + 
                        df.end.astype(str))


        df2 = cpb.gencode.make_splice_junction_df('DDX11L1.gtf')
        assert_frame_equal(df, df2)
    
    def test_neg(self):
        df = pd.DataFrame(index=[u'chr1:91630-92090:-', 
                                 u'chr1:92241-112699:-',
                                 u'chr1:112805-120774:-',
                                 u'chr1:112805-120720:-',
                                 u'chr1:120933-129054:-',
                                 u'chr1:111358-112699:-',
                                 u'chr1:112805-129054:-',
                                 u'chr1:129224-133373:-'])

        df['gene'] = 'ENSG00000238009.2'
        df['chrom'] = 'chr1'
        df['start'] = [91630, 92241, 112805, 112805, 120933, 111358, 112805,
                       129224]
        df['end'] = [92090, 112699, 120774, 120720, 129054, 112699, 129054,
                     133373]
        df['strand'] = '-'
        df['chrom:start'] = df.chrom + ':' + df.start.astype(str)
        df['chrom:end'] = df.chrom + ':' + df.end.astype(str)
        df['donor'] = df['chrom:end'] + ':' + df.strand
        df['acceptor'] = df['chrom:start'] + ':' + df.strand
        df['intron'] = (df.chrom + ':' + df.start.astype(str) + '-' + 
                        df.end.astype(str))


        df2 = cpb.gencode.make_splice_junction_df('RP11-34P13.7.gtf')
        assert_frame_equal(df, df2)
