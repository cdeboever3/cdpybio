import pandas as pd
from pandas.util.testing import assert_series_equal
from pandas.util.testing import assert_frame_equal
import pytest

import cdpybio as cpb

class TestParseAlignmentSummaryMetrics:
    def test_parse_alignment_summary_metrics(self):
        df = cpb.picard.parse_alignment_summary_metrics(
            'alignment_summary_metrics.txt')
        df2 = pd.read_pickle('alignment_summary_metrics_output.pickle')
        assert_frame_equal(df, df2)

class TestParseMarkDuplicateMetrics:
    def test_parse_mark_duplicate_metrics(self):
        metrics, hist = cpb.picard.parse_mark_duplicate_metrics(
            'duplicate_metrics.txt')
        metrics2 = pd.read_pickle('duplicate_metrics_metrics.pickle')
        hist2 = pd.read_pickle('duplicate_metrics_hist.pickle')
        assert_series_equal(metrics, metrics2)
        assert_series_equal(hist, hist2)

class TestParseInsertMetrics:
    def test_parse_insert_metrics(self):
        metrics, hist = cpb.picard.parse_insert_metrics(
            'insert_size_metrics.txt')
        metrics2 = pd.read_pickle('insert_size_metrics_metrics.pickle')
        hist2 = pd.read_pickle('insert_size_metrics_hist.pickle')
        assert_series_equal(metrics, metrics2)
        assert_series_equal(hist, hist2)

class TestParseRnaSeqMetrics:
    def test_parse_rna_seq_metrics(self):
        metrics, hist = cpb.picard.parse_rna_seq_metrics('rna_seq_metrics.txt')
        metrics2 = pd.read_pickle('rna_seq_metrics_metrics.pickle')
        hist2 = pd.read_pickle('rna_seq_metrics_hist.pickle')
        assert_series_equal(metrics, metrics2)
        assert_series_equal(hist, hist2)
