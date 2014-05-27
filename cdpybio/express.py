import sys

import pandas as pd

def combine_express_output(fnL,
                           column='eff_counts',
                           namesL=None,
                           tgN=None,
                           define_sample_name=None,
                           debug=False):
    """
    Combine eXpress output files

    Parameters:
    -----------

    fnL : list of strs of filenames
        List of paths to results.xprs files.

    column : string
        Column name of eXpress output to combine.

    namesL : list of strings
        Names to use for columns of output files. Overrides define_sample_name 
        if provided.

    tgN : string
        File with transcript-to-gene mapping. Transcripts should be in first
        column and genes in second column. 

    define_sample_name : function that takes string as input
        Function mapping filename to sample name (or basename). For instance,
        you may have the basename in the path and use a regex to extract it.
        The basenames will be used as the column names. If this is not provided,
        the columns will be named as the input files.

    debug : boolean
        Passing True will trigger any debugging statements.

    """
    if namesL is not None:
        assert len(namesL) == len(fnL)
    if define_sample_name is None:
        define_sample_name = lambda x: x
    
    transcriptL = []
    for i,fn in enumerate(fnL):
        if namesL is not None:
            bn = namesL[i]
        else:
            bn = define_sample_name(fn)
        tDF = pd.read_table(fn, index_col=1, header=0)
        se = tDF[column]
        se.name = bn
        transcriptL.append(se)
    transcriptDF = pd.DataFrame(transcriptL).T
    transcriptDF.index.name = 'transcript'
    # There should not be any missing values.
    if transcriptDF.shape != transcriptDF.dropna().shape:
        sys.stderr.write('''Missing values in eXpress output. Check that the
                         same reference was used for all output files.\n''')
        sys.exit(1)

    if tgN is not None:
        tgDF = pd.read_table(tgN,
                             index_col=0,
                             header=None,
                             names=['gene_id'])
        import copy
        geneDF = copy.deepcopy(transcriptDF)
        geneDF['gene'] = tgDF.ix[geneDF.index]
        geneDF = geneDF.groupby('gene').sum()
        return transcriptDF, geneDF
    else:
        return transcriptDF, None
