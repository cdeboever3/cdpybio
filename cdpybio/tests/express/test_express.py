from copy import deepcopy
import pandas as pd
import pytest

import cdpybio as cpb

FL = [ 'results.{}.xprs'.format(x) for x in ['a','b','c'] ]
    
TDF = pd.DataFrame([[ 68.,  44.,  50.],
                    [ 98.,  33.,  25.],
                    [ 82.,  27.,  24.]],
                   index=['TA', 'TB', 'TC'],
                   columns=['results.a.xprs', 'results.b.xprs', 
                           'results.c.xprs'])

GDF = pd.DataFrame([[ 166.,   77.,   75.],
                    [  82.,   27.,   24.]],
                   index=['GA', 'GB'],
                   columns=['results.a.xprs', 'results.b.xprs',
                            'results.c.xprs'])

class TestCombineExpressOutput:
    def test_normal(self):
        assert cpb.express.combine_express_output(FL)[0] == TDF
    
    def test_none(self):
        assert cpb.express.combine_express_output(FL)[1] is None
    
    def test_est_counts(self):
        df = TODO
        assert cpb.express.combine_express_output(FL, column='est_counts')[0] == df
    
    def test_gene(self):
        assert cpb.express.combine_express_output(FL, tgN='tg.tsv')[1] == GDF
    
    def test_missing_values(self):
        with pytest.raises(SystemExit):
            cpb.express.combine_express_output(['results.a.xprs','results.missing.xprs'])
    
    def test_names(self):
        df = deepcopy(TDF)
        df.columns = ['a','b','c']
        assert cpb.express.combine_express_output(FL, names=['a','b','c'])[0] == df
        df = deepcopy(GDF)
        df.columns = ['a','b','c']
        assert cpb.express.combine_express_output(FL, 
                                          names=['a','b','c'], 
                                          tg='tg.tsv')[1] == df
    
    def test_define_sample_names(self):
        df = deepcopy(TDF)
        df.columns = ['a','b','c']
        fnc = lambda x: x.split('.')[1]
        assert cpb.express.combine_express_output(FL, define_sample_name=fnc)
        df = deepcopy(GDF)
        df.columns = ['a','b','c']
        assert cpb.express.combine_express_output(FL, 
                                          names=['a','b','c'], 
                                          tg='tg.tsv')[1] == df
