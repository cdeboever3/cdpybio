from copy import deepcopy
import pandas as pd
from pandas.util.testing import assert_frame_equal
import pytest

import cdpybio as cpb

FL = [ 'results.{}.xprs'.format(x) for x in ['a','b','c'] ]
    
TDF = pd.DataFrame([[ 68.,  44.,  50.],
                    [ 98.,  33.,  25.],
                    [ 82.,  27.,  24.]],
                   index=['TA', 'TB', 'TC'],
                   columns=['results.a.xprs', 'results.b.xprs', 
                           'results.c.xprs'])
TDF.index.name = 'transcript'

GDF = pd.DataFrame([[ 166.,   77.,   75.],
                    [  82.,   27.,   24.]],
                   index=['GA', 'GB'],
                   columns=['results.a.xprs', 'results.b.xprs',
                            'results.c.xprs'])
GDF.index.name = 'gene'

class TestCombineExpressOutput:
    def test_normal(self):
        df = cpb.express.combine_express_output(FL)[0]
        assert_frame_equal(df, TDF, check_names = True)
    
    def test_none(self):
        assert cpb.express.combine_express_output(FL)[1] is None
    
    def test_est_counts(self):
        df = pd.DataFrame([[ 25.,  54.,  81.],
                           [ 34.,  78.,  69.],
                           [ 71.,  88.,  35.]],
                          index=TDF.index,
                          columns=TDF.columns)
        df.index.name = 'transcript'
        df2 = cpb.express.combine_express_output(FL, column='est_counts')[0]
        assert_frame_equal(df, df2)
        df = pd.DataFrame([[  59.,  132.,  150.],
                           [  71.,   88.,   35.]],
                          index=GDF.index,
                          columns=GDF.columns)
        df.index.name = 'gene'
        df2 = cpb.express.combine_express_output(FL, column='est_counts',
                                                 tg='tg.tsv')[1]
        assert_frame_equal(df, df2)
    
    def test_gene(self):
        df = cpb.express.combine_express_output(FL, tg='tg.tsv')[1] 
        assert_frame_equal(df, GDF)
    
    def test_missing_values(self):
        with pytest.raises(SystemExit):
            cpb.express.combine_express_output(['results.a.xprs',
                                                'results.missing.xprs'])
    
    def test_names(self):
        df = deepcopy(TDF)
        df.columns = ['a','b','c']
        df2 = cpb.express.combine_express_output(FL, names=['a','b','c'])[0] 
        assert_frame_equal(df, df2)
        df = deepcopy(GDF)
        df.columns = ['a','b','c']
        df2 = cpb.express.combine_express_output(FL, names=['a','b','c'], 
                                                 tg='tg.tsv')[1]
        assert_frame_equal(df, df2)
    
    def test_define_sample_names(self):
        df = deepcopy(TDF)
        df.columns = ['a','b','c']
        fnc = lambda x: x.split('.')[1]
        df2 = cpb.express.combine_express_output(FL, define_sample_name=fnc)[0]
        assert_frame_equal(df, df2)
        df = deepcopy(GDF)
        df.columns = ['a','b','c']
        df2 = cpb.express.combine_express_output(FL, names=['a','b','c'], 
                                                 tg='tg.tsv')[1]
        assert_frame_equal(df, df2)
