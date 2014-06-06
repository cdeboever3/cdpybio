from copy import deepcopy
from numpy import array
import pandas as pd
from pandas.util.testing import assert_frame_equal
from pandas.util.testing import assert_panel_equal
import pytest

import cdpybio as cpb

# TODO: I might want to include some more tests. I haven't tested whether the
# stats in the statsfiles are correct. I might want to check to make sure the
# results aren't sensitive to strand. I could also test the define_sample_name
# functionality.

EXTDF = pd.DataFrame([['gene1', 'chr1', 10, 20, '+', 'chr1:10', 'chr1:20',
                       'chr1:10:+', 'chr1:20:+', 'chr1:10-20'],
                      ['gene1', 'chr1', 5, 25, '+', 'chr1:5', 'chr1:25',
                       'chr1:5:+', 'chr1:25:+', 'chr1:5-25'],
                      ['gene1', 'chr1', 2, 20, '+', 'chr1:2', 'chr1:20',
                       'chr1:2:+', 'chr1:20:+', 'chr1:2-20'],
                      ['gene1', 'chr1', 5, 20, '+', 'chr1:5', 'chr1:20',
                       'chr1:5:+', 'chr1:20:+', 'chr1:5-20']],
                     index=['chr1:10-20:+', 'chr1:5-25:+', 'chr1:2-20:+',
                            'chr1:5-20:+'],
                     columns=['gene', 'chrom', 'start', 'end', 'strand',
                              'chr:start', 'chr:end', 'donor', 'acceptor',
                              'intron'])

class TestMisc:
    def test_read_ext(self):
        df = cpb.star.read_external_annotation('ext.tsv')
        assert_frame_equal(df, EXTDF)

    def test_read_sj_out(self):
        df = pd.DataFrame([['chr1', 2, 20, '+', 'GT/AG', True, 5, 1, 10],
                           ['chr1', 5, 20, '+', 'GT/AG', True, 20, 1, 14],
                           ['chr1', 5, 25, '+', 'CT/AC', True, 10, 1, 7],
                           ['chr1', 10, 20, '+', 'CT/AC', True, 20, 1, 7]],
                          columns=[u'chrom', u'first_bp_intron',
                                   u'last_bp_intron', u'strand',
                                   u'intron_motif', u'annotated',
                                   u'unique_junction_reads',
                                   u'multimap_junction_reads', u'max_overhang'])
        df2 = cpb.star.read_sj_out_tab('SJ.out.tab.nonew_a')
        assert_frame_equal(df, df2)

class TestMakeSJOutDict:
    def test_make_sj_out_dict(self):
        d = cpb.star.make_sj_out_dict(['SJ.out.tab.nonew_a',
                                       'SJ.out.tab.nonew_b'])
        a = cpb.star.read_sj_out_tab('SJ.out.tab.nonew_a')
        a.index = a.apply(lambda x: cpb.star._sj_out_junction(x), axis=1)
        b = cpb.star.read_sj_out_tab('SJ.out.tab.nonew_b')
        b.index = b.apply(lambda x: cpb.star._sj_out_junction(x), axis=1)
        assert_frame_equal(a, d['SJ.out.tab.nonew_a'])
        assert_frame_equal(b, d['SJ.out.tab.nonew_b'])

class TestMakeSJOutPanel:
    def test_make_sj_out_panel(self):
        ind = [u'chr1:5-20', u'chr1:5-25', u'chr1:10-20']
        d = cpb.star.make_sj_out_dict(['SJ.out.tab.nonew_a',
                                       'SJ.out.tab.nonew_b'])
        df = d['SJ.out.tab.nonew_a'].ix[ind, cpb.star.COUNT_COLS]
        df2 = d['SJ.out.tab.nonew_b'].ix[ind, cpb.star.COUNT_COLS]
        df2 = df2.fillna(0)

        p = pd.Panel({'SJ.out.tab.nonew_a':df,
                      'SJ.out.tab.nonew_b':df2})
        p = p.astype(int)
        a = pd.DataFrame([['chr1', 5, 20, '+', 'GT/AG', True],
                          ['chr1', 5, 25, '+', 'CT/AC', True],
                          ['chr1', 10, 20, '+', 'CT/AC', True]],
                         index=[u'chr1:5-20', u'chr1:5-25', u'chr1:10-20'],
                         columns=[u'chrom', u'first_bp_intron',
                                  u'last_bp_intron', u'strand', u'intron_motif',
                                  u'annotated'])
        p2, a2 = cpb.star.make_sj_out_panel(d)
        assert_frame_equal(a, a2)
        assert_panel_equal(p, p2)

    def test_new_junctions(self):
        ind = [u'chr1:2-25', u'chr1:30-40', u'chr1:5-30', u'chr1:5-20',
               u'chr1:3-25', u'chr1:10-20']
        d = cpb.star.make_sj_out_dict(['SJ.out.tab.nonew_a',
                                       'SJ.out.tab.new'])
        df = d['SJ.out.tab.nonew_a'].ix[ind, cpb.star.COUNT_COLS]
        df = df.fillna(0)
        df2 = d['SJ.out.tab.new'].ix[ind, cpb.star.COUNT_COLS]
        df2 = df2.fillna(0)

        p = pd.Panel({'SJ.out.tab.nonew_a':df,
                      'SJ.out.tab.new':df2})
        p = p.astype(int)
        a = pd.DataFrame(
            [['chr1', 2, 25, '+', 'GT/AG', False],
             ['chr1', 30, 40, '+', 'CT/AC', False],
             ['chr1', 5, 30, '+', 'GT/AG', False],
             ['chr1', 5, 20, '+', 'GT/AG', True],
             ['chr1', 3, 25, '+', 'CT/AC', False],
             ['chr1', 10, 20, '+', 'CT/AC', True]],
            index=[u'chr1:2-25', u'chr1:30-40', u'chr1:5-30', u'chr1:5-20',
                   u'chr1:3-25', u'chr1:10-20'],
            columns=[u'chrom', u'first_bp_intron', u'last_bp_intron', u'strand',
                     u'intron_motif', u'annotated']
        )
        p2, a2 = cpb.star.make_sj_out_panel(d)
        assert_frame_equal(a, a2)
        assert_panel_equal(p, p2)

class TestFilterJxnsDonorAcceptor:
    def test_filter_jxns_donor_acceptor(self):
        d = cpb.star.make_sj_out_dict(['SJ.out.tab.nonew_a',
                                       'SJ.out.tab.nonew_b'])
        p, a = cpb.star.make_sj_out_panel(d)
        ext = cpb.star.read_external_annotation('ext.tsv')
        c2, a2 = cpb.star.filter_jxns_donor_acceptor(p, a, ext)
        a = pd.DataFrame(
            [['chr1', 5, 20, '+', 'GT/AG', True, True, 'chr1:5', 
              'chr1:20', 'gene1', 'chr1:5:+', 'chr1:20:+', False, False], 
             ['chr1', 5, 25, '+', 'CT/AC', True, True, 'chr1:5', 
              'chr1:25', 'gene1', 'chr1:5:+', 'chr1:25:+', False, False], 
             ['chr1', 10, 20, '+', 'CT/AC', True, True, 'chr1:10', 
              'chr1:20', 'gene1', 'chr1:10:+', 'chr1:20:+', False, False]],
            index=[u'chr1:5-20:+', u'chr1:5-25:+', u'chr1:10-20:+'],
            columns=[u'chrom', u'first_bp_intron', u'last_bp_intron', 
                     u'strand', u'intron_motif', u'annotated', 
                     u'ext_annotated', u'chr:start', u'chr:end', 
                     u'gene_id', u'donor', u'acceptor', u'novel_donor', 
                     u'novel_acceptor'])
        c = pd.DataFrame(array([[20,  0],[10, 10],[20, 20]]),
                         index=[u'chr1:5-20:+', u'chr1:5-25:+',
                                u'chr1:10-20:+'],
                         columns=[u'SJ.out.tab.nonew_a', u'SJ.out.tab.nonew_b'])

        assert_frame_equal(a, a2)
        assert_frame_equal(c, c2)
