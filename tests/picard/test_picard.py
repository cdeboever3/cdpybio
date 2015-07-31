import os
import pandas as pd
from pandas.util.testing import assert_series_equal
from pandas.util.testing import assert_frame_equal
import pytest

import cdpybio as cpb

def add_root(fn):
    return os.path.join(cpb._root, 'tests', 'picard', fn)

class TestParseAlignmentSummaryMetrics:
    def test_parse_alignment_summary_metrics(self):
        df = cpb.picard.parse_alignment_summary_metrics(
            add_root('alignment_summary_metrics.txt'))
        df2 = pd.read_pickle(
            add_root('alignment_summary_metrics_output.pickle'))
        assert_frame_equal(df, df2)

class TestParseMarkDuplicateMetrics:
    def test_parse_mark_duplicate_metrics(self):
        metrics, hist = cpb.picard.parse_mark_duplicate_metrics(
            add_root('duplicate_metrics.txt'))
        metrics2 = pd.read_pickle(add_root('duplicate_metrics_metrics.pickle'))
        hist2 = pd.read_pickle(add_root('duplicate_metrics_hist.pickle'))
        assert_series_equal(metrics, metrics2)
        assert_series_equal(hist, hist2)

class TestParseInsertMetrics:
    def test_parse_insert_metrics(self):
        metrics, hist = cpb.picard.parse_insert_metrics(
            add_root('insert_size_metrics.txt'))
        metrics2 = pd.read_pickle(
            add_root('insert_size_metrics_metrics.pickle'))
        hist2 = pd.read_pickle(add_root('insert_size_metrics_hist.pickle'))
        assert_series_equal(metrics, metrics2)
        assert_series_equal(hist, hist2)

class TestParseRnaSeqMetrics:
    def test_parse_rna_seq_metrics(self):
        metrics, hist = cpb.picard.parse_rna_seq_metrics(
            add_root('rna_seq_metrics.txt'))
        metrics2 = pd.read_pickle(add_root('rna_seq_metrics_metrics.pickle'))
        hist2 = pd.read_pickle(add_root('rna_seq_metrics_hist.pickle'))
        assert_series_equal(metrics, metrics2)
        assert_series_equal(hist, hist2)
