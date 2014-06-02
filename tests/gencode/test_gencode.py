from copy import deepcopy
import pandas as pd
from pandas.util.testing import assert_frame_equal
import pytest

import cdpybio as cpb

# class TestMakeTranscriptGeneSe:
#     2
# 
# class TestMakeGeneInfoDf:
#     2
# 
class TestMakeSpliceJunctionDf:
    def test_plus(self):
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
