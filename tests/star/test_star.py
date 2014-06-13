from copy import deepcopy
from numpy import array
from numpy import nan
import numpy as np
import pandas as pd
from pandas.util.testing import assert_frame_equal
from pandas.util.testing import assert_panel_equal
import pytest

import cdpybio as cpb

# TODO: I might want to include some more tests. I haven't tested whether the
# stats in the statsfiles are correct. I might want to check to make sure the
# results aren't sensitive to strand. I could also test the define_sample_name
# functionality.

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
        df2, stats = cpb.star.read_external_annotation('ext.tsv')
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
        df2 = cpb.star.read_sj_out_tab('SJ.out.tab.nonew_a')
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
        df2 = cpb.star.read_sj_out_tab('SJ.out.tab.neg_nonew_a')
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
        df2 = cpb.star.read_sj_out_tab('SJ.out.tab.unk_nonew_a')
        assert_frame_equal(df, df2)

class TestMakeSJOutDict:
    def test_make_sj_out_dict_pos(self):
        d = cpb.star._make_sj_out_dict(['SJ.out.tab.nonew_a',
                                       'SJ.out.tab.nonew_b'])
        a = cpb.star.read_sj_out_tab('SJ.out.tab.nonew_a')
        a.index = a.apply(lambda x: cpb.star._sj_out_junction(x), axis=1)
        b = cpb.star.read_sj_out_tab('SJ.out.tab.nonew_b')
        b.index = b.apply(lambda x: cpb.star._sj_out_junction(x), axis=1)
        assert_frame_equal(a, d['SJ.out.tab.nonew_a'])
        assert_frame_equal(b, d['SJ.out.tab.nonew_b'])

    def test_make_sj_out_dict_neg(self):
        d = cpb.star._make_sj_out_dict(['SJ.out.tab.neg_nonew_a',
                                       'SJ.out.tab.neg_nonew_b'])
        a = cpb.star.read_sj_out_tab('SJ.out.tab.neg_nonew_a')
        a.index = a.apply(lambda x: cpb.star._sj_out_junction(x), axis=1)
        b = cpb.star.read_sj_out_tab('SJ.out.tab.neg_nonew_b')
        b.index = b.apply(lambda x: cpb.star._sj_out_junction(x), axis=1)
        assert_frame_equal(a, d['SJ.out.tab.neg_nonew_a'])
        assert_frame_equal(b, d['SJ.out.tab.neg_nonew_b'])

class TestMakeSJOutPanel:
    def test_make_sj_out_panel_pos(self):
        ind = [u'chr1:5-20', u'chr1:5-25', u'chr1:10-20']
        d = cpb.star._make_sj_out_dict(['SJ.out.tab.nonew_a',
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
                         index=ind,
                         columns=[u'chrom', u'start',
                                  u'end', u'strand', u'intron_motif',
                                  u'annotated'])
        p2, a2 = cpb.star._make_sj_out_panel(d)
        assert_frame_equal(a, a2)
        assert_panel_equal(p, p2)

    def test_make_sj_out_panel_neg(self):
        ind = [u'chr2:5-20', u'chr2:5-25', u'chr2:10-20']
        d = cpb.star._make_sj_out_dict(['SJ.out.tab.neg_nonew_a',
                                       'SJ.out.tab.neg_nonew_b'])
        df = d['SJ.out.tab.neg_nonew_a'].ix[ind, cpb.star.COUNT_COLS]
        df2 = d['SJ.out.tab.neg_nonew_b'].ix[ind, cpb.star.COUNT_COLS]
        df2 = df2.fillna(0)

        p = pd.Panel({'SJ.out.tab.neg_nonew_a':df,
                      'SJ.out.tab.neg_nonew_b':df2})
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
        d = cpb.star._make_sj_out_dict(['SJ.out.tab.nonew_a',
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
        d = cpb.star._make_sj_out_dict(['SJ.out.tab.neg_nonew_a',
                                       'SJ.out.tab.neg_new'])
        df = d['SJ.out.tab.neg_nonew_a'].ix[ind, cpb.star.COUNT_COLS]
        df = df.fillna(0)
        df2 = d['SJ.out.tab.neg_new'].ix[ind, cpb.star.COUNT_COLS]
        df2 = df2.fillna(0)

        p = pd.Panel({'SJ.out.tab.neg_nonew_a':df,
                      'SJ.out.tab.neg_new':df2})
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
        d = cpb.star._make_sj_out_dict(['SJ.out.tab.nonew_a',
                                       'SJ.out.tab.nonew_b'])
        p, a = cpb.star._make_sj_out_panel(d)
        ext, stats = cpb.star.read_external_annotation('ext.tsv')
        c2, a2, stats = cpb.star.filter_jxns_donor_acceptor(p, a, ext)
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
                         columns=[u'SJ.out.tab.nonew_a', u'SJ.out.tab.nonew_b'])

        assert_frame_equal(a, a2)
        assert_frame_equal(c, c2)
    
    def test_filter_jxns_donor_acceptor_neg(self):
        d = cpb.star._make_sj_out_dict(['SJ.out.tab.neg_nonew_a',
                                       'SJ.out.tab.neg_nonew_b'])
        p, a = cpb.star._make_sj_out_panel(d)
        ext, stats = cpb.star.read_external_annotation('ext.tsv')
        c2, a2, stats = cpb.star.filter_jxns_donor_acceptor(p, a, ext)
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
                         columns=[u'SJ.out.tab.neg_nonew_a', 
                                  u'SJ.out.tab.neg_nonew_b'])

        assert_frame_equal(a, a2)
        assert_frame_equal(c, c2)
    
    def test_filter_new_jxns(self):
        d = cpb.star._make_sj_out_dict(['SJ.out.tab.new',
                                       'SJ.out.tab.nonew_a'])
        p, a = cpb.star._make_sj_out_panel(d)
        ext, stats = cpb.star.read_external_annotation('ext.tsv')
        c2, a2, stats = cpb.star.filter_jxns_donor_acceptor(p, a, ext)
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
                         columns=[u'SJ.out.tab.new', u'SJ.out.tab.nonew_a'])

        assert_frame_equal(a, a2)
        assert_frame_equal(c, c2)

    def test_filter_new_jxns_neg(self):
        d = cpb.star._make_sj_out_dict(['SJ.out.tab.neg_new',
                                       'SJ.out.tab.neg_nonew_a'])
        p, a = cpb.star._make_sj_out_panel(d)
        ext, es = cpb.star.read_external_annotation('ext.tsv')
        c2, a2, s2 = cpb.star.filter_jxns_donor_acceptor(p, a, ext)
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
                         columns=[u'SJ.out.tab.neg_new', 
                                  u'SJ.out.tab.neg_nonew_a'])

        assert_frame_equal(a, a2)
        assert_frame_equal(c, c2)

class TestFindNovelDonorAcceptorDist:
    def test_make_splice_targets_dict_donor(self):
        df, stats = cpb.star.read_external_annotation('ext.tsv')
        strand = '+'
        feature = 'donor'
        d = cpb.star._make_splice_targets_dict(df, feature, strand)
        d2 = {'chr1:10:+': array([20]),
              'chr1:2:+': array([20]),
              'chr1:5:+': array([20, 25])}
        assert d.keys() == d2.keys()
        for k in d.keys():
            assert (d[k] == d2[k]).all()

    def test_make_splice_targets_dict_acceptor(self):
        df, stats = cpb.star.read_external_annotation('ext.tsv')
        strand = '+'
        feature = 'acceptor'
        d = cpb.star._make_splice_targets_dict(df, feature, strand)
        d2 = {'chr1:20:+': array([2, 5, 10]), 
              'chr1:25:+': array([5])}
        assert d.keys() == d2.keys()
        for k in d.keys():
            assert (d[k] == d2[k]).all()

    def test_dist_to_annot_donor_acceptor(self):
        ext, stats = cpb.star.read_external_annotation('ext.tsv')
        strand = '+'
        feature = 'donor'
        # d is a dict whose keys are donors and whose values are sets that
        # contain the positions of all acceptors associated with this donor.
        d = cpb.star._make_splice_targets_dict(ext, feature, strand)

        sjd = cpb.star._make_sj_out_dict(['SJ.out.tab.new',
                                         'SJ.out.tab.nonew_a'])
        p, a = cpb.star._make_sj_out_panel(sjd)
        c, a, s = cpb.star.filter_jxns_donor_acceptor(p, a, ext)
        novel_feature = 'acceptor'
        a = a[(a.strand == strand) & (a.novel_acceptor)]
        up, down = cpb.star._dist_to_annot_donor_acceptor(a, d, strand, 
                                                          novel_feature)
        assert up == [5]
        assert down == [np.nan]

    def test_dist_to_annot_donor_acceptor(self):
        ext, stats = cpb.star.read_external_annotation('ext.tsv')
        strand = '+'
        feature = 'acceptor'
        # d is a dict whose keys are acceptors and whose values are sets that
        # contain the positions of all donors associated with this acceptor.
        d = cpb.star._make_splice_targets_dict(ext, feature, strand)

        sjd = cpb.star._make_sj_out_dict(['SJ.out.tab.new',
                                         'SJ.out.tab.nonew_a'])
        p, a = cpb.star._make_sj_out_panel(sjd)
        c, a, s = cpb.star.filter_jxns_donor_acceptor(p, a, ext)
        novel_feature = 'donor'
        a = a[(a.strand == strand) & (a.novel_donor)]
        up, down = cpb.star._dist_to_annot_donor_acceptor(a, d, strand, 
                                                          novel_feature)
        assert up == [np.nan]
        assert down == [2]

    def test_find_novel_donor_acceptor_dist(self):
        ext, stats = cpb.star.read_external_annotation('ext.tsv')
        sjd = cpb.star._make_sj_out_dict(['SJ.out.tab.new',
                                         'SJ.out.tab.nonew_a'])
        p, a = cpb.star._make_sj_out_panel(sjd)
        c, a, s = cpb.star.filter_jxns_donor_acceptor(p, a, ext)
        df = cpb.star.find_novel_donor_acceptor_dist(a, ext)

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
                             'chr1:30:+', False, True, 5.0, nan, nan, nan],
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
