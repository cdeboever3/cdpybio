from copy import deepcopy
from numpy import array
import pandas as pd
from pandas.util.testing import assert_frame_equal
from pandas.util.testing import assert_panel_equal
import pytest

import cdpybio as cpb

EXTDF = pd.DataFrame([['gene1', 'chr1', 10, 20, '-', 'chr1:10', 'chr1:20',
                       'chr1:20:-', 'chr1:10:-', 'chr1:10:20'],
                      ['gene1', 'chr1', 5, 25, '-', 'chr1:5', 'chr1:25',
                       'chr1:25:-', 'chr1:5:-', 'chr1:5:25'],
                      ['gene1', 'chr1', 2, 20, '-', 'chr1:2', 'chr1:20',
                       'chr1:20:-', 'chr1:2:-', 'chr1:2:20'],
                      ['gene1', 'chr1', 5, 20, '-', 'chr1:5', 'chr1:20',
                       'chr1:20:-', 'chr1:5:-', 'chr1:5:20']],
                     index=['chr1:10-20:-', 'chr1:5-25:-', 'chr1:2-20:-',
                            'chr1:5-20:-'],
                     columns=['gene', 'chrom', 'start', 'end', 'strand',
                              'chr:start', 'chr:end', 'donor', 'acceptor',
                              'intron'])

# TODO: Things to test: Filtering. Make sure junctions not observed are given
# zero values. Inclusion of three types of novel junctions (donor, acceptor,
# combination). Statsfiles? Make sure filtering works correctly. Everything
# works for both strands? Test define_sample_name?

class TestMisc:
    def test_read_ext(self):
        df = cpb.star.read_external_annotation('ext.tsv')
        assert_frame_equal(df, EXTDF)

    def test_read_sj_out(self):
        df = pd.DataFrame([['chr1', 2, 20, 1, 'GT/AG', True, 5, 1, 10],
                           ['chr1', 5, 20, 1, 'GT/AG', True, 20, 1, 14],
                           ['chr1', 5, 25, 1, 'CT/AC', True, 10, 1, 7],
                           ['chr1', 10, 20, 1, 'CT/AC', True, 20, 1, 7]],
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
        col = [u'unique_junction_reads', u'multimap_junction_reads',
               u'max_overhang']
        df = pd.DataFrame(
            array([[20.0, 1.0, 14.0],
                   [10.0, 1.0, 7.0],
                   [20.0, 1.0, 7.0]], dtype=object),
            index=ind,
            columns=col
        )
        df2 = pd.DataFrame(
            array([[0.0, 0.0, 0.0],
                   [10.0, 1.0, 7.0],
                   [20.0, 1.0, 7.0]], dtype=object),
            index=ind,
            columns=col
        )

        p = pd.Panel({'SJ.out.tab.nonew_a':df,
                      'SJ.out.tab.nonew_b':df2})
        a = pd.DataFrame([['chr1', 5, 20, 'GT/AG', True],
                          ['chr1', 5, 25, 'CT/AC', True],
                          ['chr1', 10, 20, 'CT/AC', True]],
                         index=[u'chr1:5-20', u'chr1:5-25', u'chr1:10-20'],
                         columns=[u'chrom', u'first_bp_intron',
                                  u'last_bp_intron', u'intron_motif',
                                  u'annotated'])
        d = cpb.star.make_sj_out_dict(['SJ.out.tab.nonew_a',
                                       'SJ.out.tab.nonew_b'])
        p2, a2 = cpb.star.make_sj_out_panel(d)
        assert_frame_equal(a, a2)
        assert_panel_equal(p, p2)

