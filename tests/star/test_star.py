from copy import deepcopy
import os

from numpy import array
from numpy import nan
import numpy as np
import pandas as pd
from pandas.util.testing import assert_frame_equal
from pandas.util.testing import assert_panel_equal
import pytest

import cdpybio as cpb

# Note: I use pos and neg in this file to refer to the plus and minus strands
# respectively.

# TODO: I might want to include some more tests. I haven't tested whether the
# stats in the statsfiles are correct. I might want to check to make sure the
# results aren't sensitive to strand. I could also test the define_sample_name
# functionality.

def add_root(fn):
    return os.path.join(cpb._root, 'tests', 'star', fn)

SJ_NEG_NEW_A = add_root('SJ.out.tab.neg_new_a')
SJ_NEG_NEW_B = add_root('SJ.out.tab.neg_new_b')
SJ_NEG_NONEW_A = add_root('SJ.out.tab.neg_nonew_a')
SJ_NEG_NONEW_B = add_root('SJ.out.tab.neg_nonew_b')
SJ_NEW = add_root('SJ.out.tab.new')
SJ_NEW_A = add_root('SJ.out.tab.new_a')
SJ_NEW_B = add_root('SJ.out.tab.new_b')
SJ_NONEW_A = add_root('SJ.out.tab.nonew_a')
SJ_NONEW_B = add_root('SJ.out.tab.nonew_b')
SJ_UNK_NONEW_A = add_root('SJ.out.tab.unk_nonew_a')

EXT = add_root('ext.tsv')

class TestMisc:
    def test_read_ext(self):
        vals = [['gene1', 'chr1', 10, 20, '+', 'chr1:10', 'chr1:20',
                 'chr1:10:+', 'chr1:20:+', 'chr1:10-20'],
                ['gene1', 'chr1', 5, 25, '+', 'chr1:5', 'chr1:25', 'chr1:5:+',
                 'chr1:25:+', 'chr1:5-25'],
                ['gene1', 'chr1', 2, 20, '+', 'chr1:2', 'chr1:20', 'chr1:2:+',
                 'chr1:20:+', 'chr1:2-20'],
                ['gene1', 'chr1', 5, 20, '+', 'chr1:5', 'chr1:20', 'chr1:5:+',
                 'chr1:20:+', 'chr1:5-20'],
                ['gene2', 'chr2', 10, 20, '-', 'chr2:10', 'chr2:20',
                 'chr2:20:-', 'chr2:10:-', 'chr2:10-20'],
                ['gene2', 'chr2', 5, 25, '-', 'chr2:5', 'chr2:25', 'chr2:25:-',
                 'chr2:5:-', 'chr2:5-25'],
                ['gene2', 'chr2', 2, 20, '-', 'chr2:2', 'chr2:20', 'chr2:20:-',
                 'chr2:2:-', 'chr2:2-20'],
                ['gene2', 'chr2', 5, 20, '-', 'chr2:5', 'chr2:20', 'chr2:20:-',
                 'chr2:5:-', 'chr2:5-20']]
        ind = [u'chr1:10-20:+', u'chr1:5-25:+', u'chr1:2-20:+', u'chr1:5-20:+',
               u'chr2:10-20:-', u'chr2:5-25:-', u'chr2:2-20:-', u'chr2:5-20:-']
        cols=['gene', 'chrom', 'start', 'end', 'strand', 'chrom:start',
              'chrom:end', 'donor', 'acceptor', 'intron']
        df = pd.DataFrame(vals, index=ind, columns=cols)
        df2, stats = cpb.star.read_external_annotation(EXT)
        assert_frame_equal(df, df2)

    def test_read_sj_out_pos(self):
        vals = [['chr1', 2, 20, '+', 'GT/AG', True, 5, 1, 10],
                ['chr1', 5, 20, '+', 'GT/AG', True, 20, 1, 14],
                ['chr1', 5, 25, '+', 'CT/AC', True, 10, 1, 7],
                ['chr1', 10, 20, '+', 'CT/AC', True, 20, 1, 7]]
        cols = [u'chrom', u'start', u'end', u'strand', u'intron_motif',
                u'annotated', u'unique_junction_reads',
                u'multimap_junction_reads', u'max_overhang']
        df = pd.DataFrame(vals, columns=cols)
        df2 = cpb.star.read_sj_out_tab(SJ_NONEW_A)
        assert_frame_equal(df, df2)
    
    def test_read_sj_out_neg(self):
        vals = [['chr2', 2, 20, '-', 'GT/AG', True, 5, 1, 10],
                ['chr2', 5, 20, '-', 'GT/AG', True, 20, 1, 14],
                ['chr2', 5, 25, '-', 'CT/AC', True, 10, 1, 7],
                ['chr2', 10, 20, '-', 'CT/AC', True, 20, 1, 7]]
        cols = [u'chrom', u'start', u'end', u'strand', u'intron_motif',
                u'annotated', u'unique_junction_reads',
                u'multimap_junction_reads', u'max_overhang']
        df = pd.DataFrame(vals, columns=cols)
        df2 = cpb.star.read_sj_out_tab(SJ_NEG_NONEW_A)
        assert_frame_equal(df, df2)

    def test_read_sj_out_unk(self):
        df = pd.DataFrame([['chr3', 2, 20, 'unk', 'GT/AG', True, 5, 1, 10],
                           ['chr3', 5, 20, 'unk', 'GT/AG', True, 20, 1, 14],
                           ['chr3', 5, 25, 'unk', 'CT/AC', True, 10, 1, 7],
                           ['chr3', 10, 20, 'unk', 'CT/AC', True, 20, 1, 7]],
                          columns=[u'chrom', u'start',
                                   u'end', u'strand',
                                   u'intron_motif', u'annotated',
                                   u'unique_junction_reads',
                                   u'multimap_junction_reads', u'max_overhang'])
        df2 = cpb.star.read_sj_out_tab(SJ_UNK_NONEW_A)
        assert_frame_equal(df, df2)

    # TODO: I'm running into some kind of error when I compare the dataframes. I
    # see some rumbling that there may be some numpy/pandas difficulties so it
    # might not be my problem.
    # def test_read_log(self):
    #     ind = [u'Started job on', u'Started mapping on', u'Finished on', 
    #            u'Mapping speed, Million of reads per hour', 
    #            u'Number of input reads', u'Average input read length', 
    #            u'Uniquely mapped reads number', u'Uniquely mapped reads %', 
    #            u'Average mapped length', u'Number of splices: Total', 
    #            u'Number of splices: Annotated (sjdb)', 
    #            u'Number of splices: GT/AG', u'Number of splices: GC/AG', 
    #            u'Number of splices: AT/AC', u'Number of splices: Non-canonical',
    #            u'Mismatch rate per base, %', u'Deletion rate per base', 
    #            u'Deletion average length', u'Insertion rate per base', 
    #            u'Insertion average length', 
    #            u'Number of reads mapped to multiple loci', 
    #            u'% of reads mapped to multiple loci', 
    #            u'Number of reads mapped to too many loci', 
    #            u'% of reads mapped to too many loci', 
    #            u'% of reads unmapped: too many mismatches', 
    #            u'% of reads unmapped: too short', u'% of reads unmapped: other']
    #     cols = [add_root('Log.final.out.a')]
    #     vals= [['Mar 06 17:38:15'], ['Mar 06 17:53:05'], ['Mar 06 20:13:16'],
    #            ['62.51'], ['146042756'], ['135'], ['103778365'], ['71.06%'],
    #            ['119.74'], ['37420413'], ['35853326'], ['36980144'], ['351650'],
    #            ['17910'], ['70709'], ['1.13%'], ['0.01%'], ['1.51'], ['0.01%'],
    #            ['1.29'], ['42173939'], ['28.88%'], ['536'], ['0.00%'],
    #            ['0.00%'], ['0.00%'], ['0.06%']]
    #     df = pd.DataFrame(vals, index=ind, columns=cols)
    #     df2 = cpb.star._read_log(add_root('Log.final.out.a'))
    #     assert_frame_equal(df, df2)

    # TODO: I'm running into some kind of error when I compare the dataframes. I
    # see some rumbling that there may be some numpy/pandas difficulties so it
    # might not be my problem.
    # def test_make_logs_df(self):
    #     cols = [u'Started job on', u'Started mapping on', u'Finished on', 
    #             u'Mapping speed, Million of reads per hour', 
    #             u'Number of input reads', u'Average input read length', 
    #             u'Uniquely mapped reads number', u'Uniquely mapped reads %', 
    #             u'Average mapped length', u'Number of splices: Total', 
    #             u'Number of splices: Annotated (sjdb)', 
    #             u'Number of splices: GT/AG', u'Number of splices: GC/AG', 
    #             u'Number of splices: AT/AC', u'Number of splices: Non-canonical',
    #             u'Mismatch rate per base, %', u'Deletion rate per base', 
    #             u'Deletion average length', u'Insertion rate per base', 
    #             u'Insertion average length', 
    #             u'Number of reads mapped to multiple loci', 
    #             u'% of reads mapped to multiple loci', 
    #             u'Number of reads mapped to too many loci', 
    #             u'% of reads mapped to too many loci', 
    #             u'% of reads unmapped: too many mismatches', 
    #             u'% of reads unmapped: too short', u'% of reads unmapped: other']
    #     ind = [add_root(x) for x in ['Log.final.out.a', u'Log.final.out.b']]
    #     vals = [['Mar 06 17:38:15', 'Mar 06 17:53:05', 'Mar 06 20:13:16', 62.51,
    #              146042756.0, 135.0, 103778365.0, 71.06, 119.74, 37420413.0,
    #              '35853326', 36980144.0, 351650.0, 17910.0, 70709.0, 1.13, 0.01,
    #              '1.51', 0.01, '1.29', 42173939.0, 28.88, 536.0, 0.0, 0.0, 0.0,
    #              0.06], 
    #             ['Mar 04 19:39:13', 'Mar 04 19:49:11', 'Mar 04 21:13:01', 84.92,
    #              118648978.0, 136.0, 105411961.0, 88.84, 132.3, 30047584.0,
    #              '29100214', 29616122.0, 351932.0, 21726.0, 57804.0, 0.69, 0.01,
    #              '1.51', 0.01, '1.25', 13141675.0, 11.08, 951.0, 0.0, 0.0, 0.0,
    #              0.08]]
    #     df = pd.DataFrame(vals, index=ind, columns=cols)
    #     df2 = cpb.star.make_logs_df(
    #         [add_root(x) for x in ['Log.final.out.a', 'Log.final.out.b']])
    #     assert_frame_equal(df, df2)

class TestMakeSJOutDict:
    def test_make_sj_out_dict_pos(self):
        d = cpb.star._make_sj_out_dict([SJ_NONEW_A,
                                       SJ_NONEW_B])
        a = cpb.star.read_sj_out_tab(SJ_NONEW_A)
        a.index = a.apply(lambda x: cpb.star._sj_out_junction(x), axis=1)
        b = cpb.star.read_sj_out_tab(SJ_NONEW_B)
        b.index = b.apply(lambda x: cpb.star._sj_out_junction(x), axis=1)
        assert_frame_equal(a, d[SJ_NONEW_A])
        assert_frame_equal(b, d[SJ_NONEW_B])

    def test_make_sj_out_dict_neg(self):
        d = cpb.star._make_sj_out_dict([SJ_NEG_NONEW_A,
                                       SJ_NEG_NONEW_B])
        a = cpb.star.read_sj_out_tab(SJ_NEG_NONEW_A)
        a.index = a.apply(lambda x: cpb.star._sj_out_junction(x), axis=1)
        b = cpb.star.read_sj_out_tab(SJ_NEG_NONEW_B)
        b.index = b.apply(lambda x: cpb.star._sj_out_junction(x), axis=1)
        assert_frame_equal(a, d[SJ_NEG_NONEW_A])
        assert_frame_equal(b, d[SJ_NEG_NONEW_B])

class TestMakeSJOutPanel:
    def test_make_sj_out_panel_pos(self):
        ind = [u'chr1:5-20', u'chr1:5-25', u'chr1:10-20']
        d = cpb.star._make_sj_out_dict([SJ_NONEW_A,
                                       SJ_NONEW_B])
        df = d[SJ_NONEW_A].ix[ind, cpb.star.COUNT_COLS]
        df2 = d[SJ_NONEW_B].ix[ind, cpb.star.COUNT_COLS]
        df2 = df2.fillna(0)

        p = pd.Panel({SJ_NONEW_A:df,
                      SJ_NONEW_B:df2})
        p = p.astype(int)
        a = pd.DataFrame([['chr1', 5, 20, '+', 'GT/AG', True],
                          ['chr1', 5, 25, '+', 'CT/AC', True],
                          ['chr1', 10, 20, '+', 'CT/AC', True]],
                         index=ind,
                         columns=[u'chrom', u'start',
                                  u'end', u'strand', u'intron_motif',
                                  u'annotated'])
        p2, a2 = cpb.star._make_sj_out_panel(d)
        assert_frame_equal(a, a2)
        assert_panel_equal(p, p2)

    def test_make_sj_out_panel_neg(self):
        ind = [u'chr2:5-20', u'chr2:5-25', u'chr2:10-20']
        d = cpb.star._make_sj_out_dict([SJ_NEG_NONEW_A,
                                       SJ_NEG_NONEW_B])
        df = d[SJ_NEG_NONEW_A].ix[ind, cpb.star.COUNT_COLS]
        df2 = d[SJ_NEG_NONEW_B].ix[ind, cpb.star.COUNT_COLS]
        df2 = df2.fillna(0)

        p = pd.Panel({SJ_NEG_NONEW_A:df,
                      SJ_NEG_NONEW_B:df2})
        p = p.astype(int)
        a = pd.DataFrame([['chr2', 5, 20, '-', 'GT/AG', True],
                          ['chr2', 5, 25, '-', 'CT/AC', True],
                          ['chr2', 10, 20, '-', 'CT/AC', True]],
                         index=ind,
                         columns=[u'chrom', u'start',
                                  u'end', u'strand', u'intron_motif',
                                  u'annotated'])
        p2, a2 = cpb.star._make_sj_out_panel(d)
        assert_frame_equal(a, a2)
        assert_panel_equal(p, p2)

    def test_new_junctions_pos(self):
        ind = [u'chr1:2-25', u'chr1:3-25', u'chr1:5-20', u'chr1:5-30',
               u'chr1:10-20', u'chr1:30-40']
        d = cpb.star._make_sj_out_dict([SJ_NONEW_A,
                                       SJ_NEW_A])
        df = d[SJ_NONEW_A].ix[ind, cpb.star.COUNT_COLS]
        df = df.fillna(0)
        df2 = d[SJ_NEW_A].ix[ind, cpb.star.COUNT_COLS]
        df2 = df2.fillna(0)

        p = pd.Panel({SJ_NONEW_A:df,
                      SJ_NEW_A:df2})
        p = p.astype(int)
        a = pd.DataFrame(
            [['chr1', 2, 25, '+', 'GT/AG', False],
             ['chr1', 3, 25, '+', 'CT/AC', False],
             ['chr1', 5, 20, '+', 'GT/AG', True],
             ['chr1', 5, 30, '+', 'GT/AG', False],
             ['chr1', 10, 20, '+', 'CT/AC', True],
             ['chr1', 30, 40, '+', 'CT/AC', False]],
            index=ind,
            columns=[u'chrom', u'start', u'end', u'strand',
                     u'intron_motif', u'annotated']
        )
        p2, a2 = cpb.star._make_sj_out_panel(d)
        assert_frame_equal(a, a2)
        assert_panel_equal(p, p2)

    def test_new_junctions_neg(self):
        ind = [u'chr2:2-25', u'chr2:3-25', u'chr2:5-20', u'chr2:5-30',
               u'chr2:10-20', u'chr2:30-40']
        d = cpb.star._make_sj_out_dict([SJ_NEG_NONEW_A,
                                       SJ_NEG_NEW_A])
        df = d[SJ_NEG_NONEW_A].ix[ind, cpb.star.COUNT_COLS]
        df = df.fillna(0)
        df2 = d[SJ_NEG_NEW_A].ix[ind, cpb.star.COUNT_COLS]
        df2 = df2.fillna(0)

        p = pd.Panel({SJ_NEG_NONEW_A:df,
                      SJ_NEG_NEW_A:df2})
        p = p.astype(int)
        a = pd.DataFrame(
            [['chr2', 2, 25, '-', 'GT/AG', False],
             ['chr2', 3, 25, '-', 'CT/AC', False],
             ['chr2', 5, 20, '-', 'GT/AG', True],
             ['chr2', 5, 30, '-', 'GT/AG', False],
             ['chr2', 10, 20, '-', 'CT/AC', True],
             ['chr2', 30, 40, '-', 'CT/AC', False]],
            index=ind,
            columns=[u'chrom', u'start', u'end', u'strand',
                     u'intron_motif', u'annotated']
        )
        p2, a2 = cpb.star._make_sj_out_panel(d)
        assert_frame_equal(a, a2)
        assert_panel_equal(p, p2)

class TestFilterJxnsDonorAcceptor:
    def test_filter_jxns_donor_acceptor_pos(self):
        d = cpb.star._make_sj_out_dict([SJ_NONEW_A,
                                       SJ_NONEW_B])
        p, a = cpb.star._make_sj_out_panel(d)
        ext, stats = cpb.star.read_external_annotation(EXT)
        c2, a2, stats = cpb.star._filter_jxns_donor_acceptor(p, a, ext)
        a = pd.DataFrame(
            [['chr1', 5, 20, '+', 'GT/AG', True, True, 'chr1:5', 
              'chr1:20', 'gene1', 'chr1:5:+', 'chr1:20:+', False, False], 
             ['chr1', 5, 25, '+', 'CT/AC', True, True, 'chr1:5', 
              'chr1:25', 'gene1', 'chr1:5:+', 'chr1:25:+', False, False], 
             ['chr1', 10, 20, '+', 'CT/AC', True, True, 'chr1:10', 
              'chr1:20', 'gene1', 'chr1:10:+', 'chr1:20:+', False, False]],
            index=[u'chr1:5-20:+', u'chr1:5-25:+', u'chr1:10-20:+'],
            columns=[u'chrom', u'start', u'end', 
                     u'strand', u'intron_motif', u'annotated', 
                     u'ext_annotated', u'chrom:start', u'chrom:end', 
                     u'gene_id', u'donor', u'acceptor', u'novel_donor', 
                     u'novel_acceptor'])
        c = pd.DataFrame(array([[20,  0],[10, 10],[20, 20]]),
                         index=[u'chr1:5-20:+', u'chr1:5-25:+',
                                u'chr1:10-20:+'],
                         columns=[SJ_NONEW_A, SJ_NONEW_B])

        assert_frame_equal(a, a2)
        assert_frame_equal(c, c2)
    
    def test_filter_jxns_donor_acceptor_neg(self):
        d = cpb.star._make_sj_out_dict([SJ_NEG_NONEW_A,
                                       SJ_NEG_NONEW_B])
        p, a = cpb.star._make_sj_out_panel(d)
        ext, stats = cpb.star.read_external_annotation(EXT)
        c2, a2, stats = cpb.star._filter_jxns_donor_acceptor(p, a, ext)
        a = pd.DataFrame(
            [['chr2', 5, 20, '-', 'GT/AG', True, True, 'chr2:5', 
              'chr2:20', 'gene2', 'chr2:20:-', 'chr2:5:-', False, False], 
             ['chr2', 5, 25, '-', 'CT/AC', True, True, 'chr2:5', 
              'chr2:25', 'gene2', 'chr2:25:-', 'chr2:5:-', False, False], 
             ['chr2', 10, 20, '-', 'CT/AC', True, True, 'chr2:10', 
              'chr2:20', 'gene2', 'chr2:20:-', 'chr2:10:-', False, False]],
            index=[u'chr2:5-20:-', u'chr2:5-25:-', u'chr2:10-20:-'],
            columns=[u'chrom', u'start', u'end', 
                     u'strand', u'intron_motif', u'annotated', 
                     u'ext_annotated', u'chrom:start', u'chrom:end', 
                     u'gene_id', u'donor', u'acceptor', u'novel_donor', 
                     u'novel_acceptor'])
        c = pd.DataFrame(array([[20,  0],[10, 10],[20, 20]]),
                         index=[u'chr2:5-20:-', u'chr2:5-25:-',
                                u'chr2:10-20:-'],
                         columns=[SJ_NEG_NONEW_A, 
                                  SJ_NEG_NONEW_B])

        assert_frame_equal(a, a2)
        assert_frame_equal(c, c2)
    
    def test_filter_new_jxns(self):
        d = cpb.star._make_sj_out_dict([SJ_NEW_A,
                                       SJ_NONEW_A])
        p, a = cpb.star._make_sj_out_panel(d)
        ext, stats = cpb.star.read_external_annotation(EXT)
        c2, a2, stats = cpb.star._filter_jxns_donor_acceptor(p, a, ext)
        a = pd.DataFrame(
                [['chr1', 2, 25, '+', 'GT/AG', False, False, 'chr1:2',
                  'chr1:25', 'gene1', 'chr1:2:+', 'chr1:25:+', False, False],
                 ['chr1', 3, 25, '+', 'CT/AC', False, False, 'chr1:3', 
                  'chr1:25', 'gene1', 'chr1:3:+', 'chr1:25:+', True, False],
                 ['chr1', 5, 20, '+', 'GT/AG', True, True, 'chr1:5', 'chr1:20',
                  'gene1', 'chr1:5:+', 'chr1:20:+', False, False],
                 ['chr1', 5, 30, '+', 'GT/AG', False, False, 'chr1:5', 
                  'chr1:30', 'gene1', 'chr1:5:+', 'chr1:30:+', False, True],
                 ['chr1', 10, 20, '+', 'CT/AC', True, True, 'chr1:10', 
                  'chr1:20', 'gene1', 'chr1:10:+', 'chr1:20:+', False, False]],
                index=[u'chr1:2-25:+', u'chr1:3-25:+', u'chr1:5-20:+',
                       u'chr1:5-30:+', u'chr1:10-20:+'],
                columns=[u'chrom', u'start', u'end', u'strand', u'intron_motif',
                         u'annotated', u'ext_annotated', u'chrom:start',
                         u'chrom:end', u'gene_id', u'donor', u'acceptor',
                         u'novel_donor', u'novel_acceptor'])

        c = pd.DataFrame([[25,  0],
                          [30,  0],
                          [ 0, 20],
                          [20,  0],
                          [ 0, 20]],
                         index=[u'chr1:2-25:+', u'chr1:3-25:+', u'chr1:5-20:+',
                                u'chr1:5-30:+', u'chr1:10-20:+'],
                         columns=[SJ_NEW_A, SJ_NONEW_A])

        assert_frame_equal(a, a2)
        assert_frame_equal(c, c2)

    def test_filter_new_jxns_neg(self):
        d = cpb.star._make_sj_out_dict([SJ_NEG_NEW_A,
                                       SJ_NEG_NONEW_A])
        p, a = cpb.star._make_sj_out_panel(d)
        ext, es = cpb.star.read_external_annotation(EXT)
        c2, a2, s2 = cpb.star._filter_jxns_donor_acceptor(p, a, ext)
        a = pd.DataFrame(
                [['chr2', 2, 25, '-', 'GT/AG', False, False, 'chr2:2',
                  'chr2:25', 'gene2', 'chr2:25:-', 'chr2:2:-', False, False],
                 ['chr2', 3, 25, '-', 'CT/AC', False, False, 'chr2:3', 
                  'chr2:25', 'gene2', 'chr2:25:-', 'chr2:3:-', False, True],
                 ['chr2', 5, 20, '-', 'GT/AG', True, True, 'chr2:5', 'chr2:20',
                  'gene2', 'chr2:20:-', 'chr2:5:-', False, False],
                 ['chr2', 5, 30, '-', 'GT/AG', False, False, 'chr2:5', 
                  'chr2:30', 'gene2', 'chr2:30:-', 'chr2:5:-', True, False],
                 ['chr2', 10, 20, '-', 'CT/AC', True, True, 'chr2:10', 
                  'chr2:20', 'gene2', 'chr2:20:-', 'chr2:10:-', False, False]],
                index=[u'chr2:2-25:-', u'chr2:3-25:-', u'chr2:5-20:-',
                       u'chr2:5-30:-', u'chr2:10-20:-'],
                columns=[u'chrom', u'start', u'end', u'strand', u'intron_motif',
                         u'annotated', u'ext_annotated', u'chrom:start',
                         u'chrom:end', u'gene_id', u'donor', u'acceptor',
                         u'novel_donor', u'novel_acceptor'])

        c = pd.DataFrame([[25,  0],
                          [30,  0],
                          [ 0, 20],
                          [20,  0],
                          [ 0, 20]],
                         index=[u'chr2:2-25:-', u'chr2:3-25:-', u'chr2:5-20:-',
                                u'chr2:5-30:-', u'chr2:10-20:-'],
                         columns=[SJ_NEG_NEW_A, 
                                  SJ_NEG_NONEW_A])

        assert_frame_equal(a, a2)
        assert_frame_equal(c, c2)

class TestFindNovelDonorAcceptorDist:
    def test_make_splice_targets_dict_donor_pos(self):
        df, stats = cpb.star.read_external_annotation(EXT)
        strand = '+'
        feature = 'donor'
        d = cpb.star._make_splice_targets_dict(df, feature, strand)
        d2 = {'chr1:10:+': array([20]),
              'chr1:2:+': array([20]),
              'chr1:5:+': array([20, 25])}
        assert d.keys() == d2.keys()
        for k in d.keys():
            assert (d[k] == d2[k]).all()

    def test_make_splice_targets_dict_donor_neg(self):
        df, stats = cpb.star.read_external_annotation(EXT)
        strand = '-'
        feature = 'donor'
        d = cpb.star._make_splice_targets_dict(df, feature, strand)
        d2 = {'chr2:20:-': array([2, 5, 10]),
              'chr2:25:-': array([5])}
        assert d.keys() == d2.keys()
        for k in d.keys():
            assert (d[k] == d2[k]).all()

    def test_make_splice_targets_dict_acceptor_pos(self):
        df, stats = cpb.star.read_external_annotation(EXT)
        strand = '+'
        feature = 'acceptor'
        d = cpb.star._make_splice_targets_dict(df, feature, strand)
        d2 = {'chr1:20:+': array([2, 5, 10]), 
              'chr1:25:+': array([5])}
        assert d.keys() == d2.keys()
        for k in d.keys():
            assert (d[k] == d2[k]).all()

    def test_make_splice_targets_dict_acceptor_neg(self):
        df, stats = cpb.star.read_external_annotation(EXT)
        strand = '-'
        feature = 'acceptor'
        d = cpb.star._make_splice_targets_dict(df, feature, strand)
        d2 = {'chr2:2:-': array([20]), 
              'chr2:5:-': array([20, 25]),
              'chr2:10:-': array([20])}
        assert d.keys() == d2.keys()
        for k in d.keys():
            assert (d[k] == d2[k]).all()

    def test_dist_to_annot_donor_acceptor_pos(self):
        ext, stats = cpb.star.read_external_annotation(EXT)
        strand = '+'
        feature = 'donor'
        # d is a dict whose keys are donors and whose values are sets that
        # contain the positions of all acceptors associated with this donor.
        d = cpb.star._make_splice_targets_dict(ext, feature, strand)

        sjd = cpb.star._make_sj_out_dict([SJ_NEW_A,
                                         SJ_NONEW_A])
        p, a = cpb.star._make_sj_out_panel(sjd)
        c, a, s = cpb.star._filter_jxns_donor_acceptor(p, a, ext)
        novel_feature = 'acceptor'
        a = a[(a.strand == strand) & (a.novel_acceptor)]
        up, down = cpb.star._dist_to_annot_donor_acceptor(a, d, strand, 
                                                          novel_feature)
        assert up == [5]
        assert down == [nan]

    def test_dist_to_annot_donor_acceptor_neg(self):
        ext, stats = cpb.star.read_external_annotation(EXT)
        strand = '-'
        feature = 'donor'
        # d is a dict whose keys are donors and whose values are sets that
        # contain the positions of all acceptors associated with this donor.
        d = cpb.star._make_splice_targets_dict(ext, feature, strand)

        sjd = cpb.star._make_sj_out_dict([SJ_NEG_NEW_A,
                                         SJ_NEG_NONEW_A])
        p, a = cpb.star._make_sj_out_panel(sjd)
        c, a, s = cpb.star._filter_jxns_donor_acceptor(p, a, ext)
        novel_feature = 'acceptor'
        a = a[(a.strand == strand) & (a.novel_acceptor)]
        up, down = cpb.star._dist_to_annot_donor_acceptor(a, d, strand, 
                                                          novel_feature)
        assert up == [2]
        assert down == [nan]

    def test_dist_to_annot_donor_acceptor_pos(self):
        ext, stats = cpb.star.read_external_annotation(EXT)
        strand = '+'
        feature = 'acceptor'
        # d is a dict whose keys are acceptors and whose values are sets that
        # contain the positions of all donors associated with this acceptor.
        d = cpb.star._make_splice_targets_dict(ext, feature, strand)

        sjd = cpb.star._make_sj_out_dict([SJ_NEW_A,
                                         SJ_NONEW_A])
        p, a = cpb.star._make_sj_out_panel(sjd)
        c, a, s = cpb.star._filter_jxns_donor_acceptor(p, a, ext)
        novel_feature = 'donor'
        a = a[(a.strand == strand) & (a.novel_donor)]
        up, down = cpb.star._dist_to_annot_donor_acceptor(a, d, strand, 
                                                          novel_feature)
        assert up == [np.nan]
        assert down == [2]

    def test_dist_to_annot_donor_acceptor_neg(self):
        ext, stats = cpb.star.read_external_annotation(EXT)
        strand = '-'
        feature = 'acceptor'
        # d is a dict whose keys are acceptors and whose values are sets that
        # contain the positions of all donors associated with this acceptor.
        d = cpb.star._make_splice_targets_dict(ext, feature, strand)

        sjd = cpb.star._make_sj_out_dict([SJ_NEG_NEW_A,
                                         SJ_NEG_NONEW_A])
        p, a = cpb.star._make_sj_out_panel(sjd)
        c, a, s = cpb.star._filter_jxns_donor_acceptor(p, a, ext)
        novel_feature = 'donor' 
        a = a[(a.strand == strand) & (a.novel_donor)]
        up, down = cpb.star._dist_to_annot_donor_acceptor(a, d, strand, 
                                                          novel_feature)
        assert up == [np.nan]
        assert down == [5]

    def test_find_novel_donor_acceptor_dist_pos_a(self):
        ext, stats = cpb.star.read_external_annotation(EXT)
        sjd = cpb.star._make_sj_out_dict([SJ_NEW_A,
                                         SJ_NONEW_A])
        p, a = cpb.star._make_sj_out_panel(sjd)
        c, a, s = cpb.star._filter_jxns_donor_acceptor(p, a, ext)
        df = cpb.star._find_novel_donor_acceptor_dist(a, ext)

        df2 = pd.DataFrame([['chr1', 2, 25, '+', 'GT/AG', False, False,
                             'chr1:2', 'chr1:25', 'gene1', 'chr1:2:+',
                             'chr1:25:+', False, False, nan, nan, nan, nan],
                            ['chr1', 3, 25, '+', 'CT/AC', False, False,
                             'chr1:3', 'chr1:25', 'gene1', 'chr1:3:+',
                             'chr1:25:+', True, False, nan, 2.0, nan, nan],
                            ['chr1', 5, 20, '+', 'GT/AG', True, True, 'chr1:5',
                             'chr1:20', 'gene1', 'chr1:5:+', 'chr1:20:+', False,
                             False, nan, nan, nan, nan], 
                            ['chr1', 5, 30, '+', 'GT/AG', False, False,
                             'chr1:5', 'chr1:30', 'gene1', 'chr1:5:+',
                             'chr1:30:+', False, True, nan, nan, 5.0, nan],
                            ['chr1', 10, 20, '+', 'CT/AC', True, True,
                             'chr1:10', 'chr1:20', 'gene1', 'chr1:10:+',
                             'chr1:20:+', False, False, nan, nan, nan, nan]],
                           index=[u'chr1:2-25:+', u'chr1:3-25:+',
                                  u'chr1:5-20:+', u'chr1:5-30:+',
                                  u'chr1:10-20:+'],
                           columns=[u'chrom', u'start', u'end', u'strand',
                                    u'intron_motif', u'annotated',
                                    u'ext_annotated', u'chrom:start',
                                    u'chrom:end', u'gene_id', u'donor',
                                    u'acceptor', u'novel_donor',
                                    u'novel_acceptor', u'upstream_donor_dist',
                                    u'downstream_donor_dist',
                                    u'upstream_acceptor_dist',
                                    u'downstream_acceptor_dist'])
        assert_frame_equal(df, df2)
    
    def test_find_novel_donor_acceptor_dist_neg_a(self):
        ext, stats = cpb.star.read_external_annotation(EXT)
        sjd = cpb.star._make_sj_out_dict([SJ_NEG_NEW_A,
                                         SJ_NEG_NONEW_A])
        p, a = cpb.star._make_sj_out_panel(sjd)
        c, a, s = cpb.star._filter_jxns_donor_acceptor(p, a, ext)
        df = cpb.star._find_novel_donor_acceptor_dist(a, ext)

        df2 = pd.DataFrame([['chr2', 2, 25, '-', 'GT/AG', False, False,
                             'chr2:2', 'chr2:25', 'gene2', 'chr2:25:-',
                             'chr2:2:-', False, False, nan, nan, nan, nan],
                            ['chr2', 3, 25, '-', 'CT/AC', False, False,
                             'chr2:3', 'chr2:25', 'gene2', 'chr2:25:-',
                             'chr2:3:-', False, True, nan, nan, 2.0, nan],
                            ['chr2', 5, 20, '-', 'GT/AG', True, True, 'chr2:5',
                             'chr2:20', 'gene2', 'chr2:20:-', 'chr2:5:-', False,
                             False, nan, nan, nan, nan], 
                            ['chr2', 5, 30, '-', 'GT/AG', False, False,
                             'chr2:5', 'chr2:30', 'gene2', 'chr2:30:-',
                             'chr2:5:-', True, False, nan, 5.0, nan, nan],
                            ['chr2', 10, 20, '-', 'CT/AC', True, True,
                             'chr2:10', 'chr2:20', 'gene2', 'chr2:20:-',
                             'chr2:10:-', False, False, nan, nan, nan, nan]],
                           index=[u'chr2:2-25:-', u'chr2:3-25:-',
                                  u'chr2:5-20:-', u'chr2:5-30:-',
                                  u'chr2:10-20:-'],
                           columns=[u'chrom', u'start', u'end', u'strand',
                                    u'intron_motif', u'annotated',
                                    u'ext_annotated', u'chrom:start',
                                    u'chrom:end', u'gene_id', u'donor',
                                    u'acceptor', u'novel_donor',
                                    u'novel_acceptor', u'upstream_donor_dist',
                                    u'downstream_donor_dist',
                                    u'upstream_acceptor_dist',
                                    u'downstream_acceptor_dist'])
        assert_frame_equal(df, df2)
    
    def test_find_novel_donor_acceptor_dist_pos_b(self):
        ext, stats = cpb.star.read_external_annotation(EXT)
        sjd = cpb.star._make_sj_out_dict([SJ_NEW_B,
                                         SJ_NONEW_A])
        p, a = cpb.star._make_sj_out_panel(sjd)
        c, a, s = cpb.star._filter_jxns_donor_acceptor(p, a, ext)
        df = cpb.star._find_novel_donor_acceptor_dist(a, ext)

        vals = \
            [['chr1', 2, 25, '+', 'GT/AG', False, False, 'chr1:2', 'chr1:25',
              'gene1', 'chr1:2:+', 'chr1:25:+', False, False, nan, nan, nan,
              nan],
             ['chr1', 3, 20, '+', 'CT/AC', False, False, 'chr1:3', 'chr1:20',
              'gene1', 'chr1:3:+', 'chr1:20:+', True, False, 1.0, 2.0, nan,
              nan],
             ['chr1', 5, 20, '+', 'GT/AG', True, True, 'chr1:5', 'chr1:20',
              'gene1', 'chr1:5:+', 'chr1:20:+', False, False, nan, nan, nan,
              nan],
             ['chr1', 5, 22, '+', 'GT/AG', False, False, 'chr1:5', 'chr1:22',
              'gene1', 'chr1:5:+', 'chr1:22:+', False, True, nan, nan, 2.0,
              3.0],
             ['chr1', 5, 30, '+', 'GT/AG', False, False, 'chr1:5', 'chr1:30',
              'gene1', 'chr1:5:+', 'chr1:30:+', False, True, nan, nan, 5.0,
              nan],
             ['chr1', 10, 20, '+', 'CT/AC', True, True, 'chr1:10', 'chr1:20',
              'gene1', 'chr1:10:+', 'chr1:20:+', False, False, nan, nan, nan,
              nan]]
        cols = [u'chrom', u'start', u'end', u'strand', u'intron_motif',
                u'annotated', u'ext_annotated', u'chrom:start', u'chrom:end',
                u'gene_id', u'donor', u'acceptor', u'novel_donor',
                u'novel_acceptor', u'upstream_donor_dist',
                u'downstream_donor_dist', u'upstream_acceptor_dist',
                u'downstream_acceptor_dist']
        ind = [u'chr1:2-25:+', u'chr1:3-20:+', u'chr1:5-20:+', u'chr1:5-22:+',
               u'chr1:5-30:+', u'chr1:10-20:+']

        df2 = pd.DataFrame(vals, index=ind, columns=cols)
        assert_frame_equal(df, df2)
    
    def test_find_novel_donor_acceptor_dist_neg_b(self):
        ext, stats = cpb.star.read_external_annotation(EXT)
        sjd = cpb.star._make_sj_out_dict([SJ_NEG_NEW_B,
                                         SJ_NEG_NONEW_A])
        p, a = cpb.star._make_sj_out_panel(sjd)
        c, a, s = cpb.star._filter_jxns_donor_acceptor(p, a, ext)
        df = cpb.star._find_novel_donor_acceptor_dist(a, ext)

        vals = [['chr2', 2, 25, '-', 'GT/AG', False, False, 'chr2:2', 'chr2:25',
                 'gene2', 'chr2:25:-', 'chr2:2:-', False, False, nan, nan, nan,
                 nan],
                ['chr2', 3, 20, '-', 'CT/AC', False, False, 'chr2:3', 'chr2:20',
                 'gene2', 'chr2:20:-', 'chr2:3:-', False, True, nan, nan, 2.0,
                 1.0],
                ['chr2', 5, 20, '-', 'GT/AG', True, True, 'chr2:5', 'chr2:20',
                 'gene2', 'chr2:20:-', 'chr2:5:-', False, False, nan, nan, nan,
                 nan],
                ['chr2', 5, 22, '-', 'GT/AG', False, False, 'chr2:5', 'chr2:22',
                 'gene2', 'chr2:22:-', 'chr2:5:-', True, False, 3.0, 2.0, nan,
                 nan],
                ['chr2', 5, 30, '-', 'GT/AG', False, False, 'chr2:5', 'chr2:30',
                 'gene2', 'chr2:30:-', 'chr2:5:-', True, False, nan, 5.0, nan,
                 nan],
                ['chr2', 10, 20, '-', 'CT/AC', True, True, 'chr2:10', 'chr2:20',
                 'gene2', 'chr2:20:-', 'chr2:10:-', False, False, nan, nan, nan,
                 nan]]
        cols = [u'chrom', u'start', u'end', u'strand', u'intron_motif',
                u'annotated', u'ext_annotated', u'chrom:start', u'chrom:end',
                u'gene_id', u'donor', u'acceptor', u'novel_donor',
                u'novel_acceptor', u'upstream_donor_dist',
                u'downstream_donor_dist', u'upstream_acceptor_dist',
                u'downstream_acceptor_dist']
        ind = [u'chr2:2-25:-', u'chr2:3-20:-', u'chr2:5-20:-', u'chr2:5-22:-',
               u'chr2:5-30:-', u'chr2:10-20:-']
        df2 = pd.DataFrame(vals, index=ind, columns=cols)

        assert_frame_equal(df, df2)
