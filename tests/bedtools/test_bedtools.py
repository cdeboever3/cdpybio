from copy import deepcopy
import os
from numpy import array
import pandas as pd
from pandas.util.testing import assert_frame_equal
import pytest

import cdpybio as cpb

class TestBedsToBoolean:
    def test_normal(self):
        beds = [os.path.join(cpb._root, 'tests', 'bedtools', x) for x in
                ['a.bed', 'b.bed', 'c.bed']]
        df = cpb.bedtools.beds_to_boolean(beds)
        tdf = pd.DataFrame(
            array([[1, 1, 1],
                   [1, 1, 1],
                   [1, 1, 1]]),
            index=[u'chr1:5-60', u'chr1:90-160', u'chr1:190-310'],
            columns=beds)
        assert_frame_equal(df, tdf, check_names = True)
    
    def test_not_all_one(self):
        beds = [os.path.join(cpb._root, 'tests', 'bedtools', x) for x in
                ['a.bed', 'b.bed', 'c.bed', 'd.bed']]
        df = cpb.bedtools.beds_to_boolean(beds)
        tdf = pd.DataFrame(
            array([[1, 1, 1, 0],
                   [1, 1, 1, 0],
                   [1, 1, 1, 0],
                   [0, 0, 0, 1]]),
            index=[u'chr1:5-60', u'chr1:90-160', u'chr1:190-310',
                   u'chr2:200-300'],
            columns=beds)
        assert_frame_equal(df, tdf, check_names = True)
    
class TestCombine:
    def test_normal(self):
        beds = [os.path.join(cpb._root, 'tests', 'bedtools', x) for x in
                ['a.bed', 'b.bed', 'c.bed']]
        bt = cpb.bedtools.combine(beds)
        df = bt.to_dataframe()
        tdf = pd.DataFrame(
            [['chr1', 5, 60],
             ['chr1', 90, 160],
             ['chr1', 190, 310]],
            index=array([0, 1, 2], dtype='object'),
            columns=[u'chrom', u'start', u'end'])
        assert_frame_equal(df, tdf, check_names=True)

class TestIntervalsToBed:
    def test_intervals_to_bed(self):
        intervals = ['chr1:100-200', 'chr2:200-300']
        cpb.bedtools.intervals_to_bed(intervals)

    def test_intervals_to_bed_strand(self):
        intervals = ['chr1:100-200:+', 'chr2:200-300:-']
        cpb.bedtools.intervals_to_bed(intervals)

    def test_intervals_to_bed_mixed(self):
        intervals = ['chr1:100-200:+', 'chr2:200-300']
        cpb.bedtools.intervals_to_bed(intervals)
