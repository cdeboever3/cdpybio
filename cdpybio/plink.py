import numpy as np
import sys

import pandas as pd

def read_linear2(fn, header=True):
    """Read a plink 2 output file of type glm.linear into a pandas DataFrame.

    Parameters
    ----------
    fn : str 
        Path to the plink file. The file can be gzipped or not.

    header : str
        True if the file has a header (this is generally the case unless the
        file has been processed after it was created). False if no header. Pass
        None if it's unknown whether the file has a header.

    Returns
    -------
    res : pandas.DataFrame
        Dataframe with results.

    """
    dtypes = {'#CHROM':str, 'POS':int, 'ID':str, 'REF':str, 'ALT1':str,
              'TEST':str, 'OBS_CT':int, 'BETA':float, 'SE':float,
              'T_STAT':float, 'P':float}
    if header is None:
        if fn[-3:] == '.gz':
            from gzip import open
            with gzip.open(fn, 'r') as f:
                line = f.readline()
        else:
            with open(fn, 'r') as f:
                line = f.readline()
        header = line [0] == '#'
    if header:
        res = pd.read_table(fn, index_col=2, dtype=dtypes, low_memory=False)
    else:
        cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT1', 'TEST', 'OBS_CT',
                'BETA', 'SE', 'T_STAT', 'P']
        res = pd.read_table(fn, index_col=2, dtype=dtypes, names=cols,
                            low_memory=False)
    res.columns = [x.replace('#', '') for x in res.columns]
    return(res)
            
def read_logistic2(fn, header=True):
    """Read a plink 2 output file of type glm.logistic into a pandas DataFrame.

    Parameters
    ----------
    fn : str 
        Path to the plink file. The file can be gzipped or not.

    header : str
        True if the file has a header (this is generally the case unless the
        file has been processed after it was created).

    Returns
    -------
    res : pandas.DataFrame
        Dataframe with results.

    """
    dtypes = {'#CHROM':str, 'POS':int, 'ID':str, 'REF':str, 'ALT1':str,
              'FIRTH?':str, 'TEST':str, 'OBS_CT':int, 'OR':float, 'SE':float,
              'T_STAT':float, 'P':float}
    if header is None:
        if fn[-3:] == '.gz':
            from gzip import open
        with open(fn, 'r') as f:
            line = f.readline()
        header = line [0] == '#'
    if header:
        res = pd.read_table(fn, index_col=2, dtype=dtypes, low_memory=False)
    else:
        cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT1', 'FIRTH?', 'TEST',
                'OBS_CT', 'OR', 'SE', 'T_STAT', 'P']
        res = pd.read_table(fn, index_col=2, dtype=dtypes, names=cols,
                            low_memory=False)
    res.columns = [x.replace('#', '') for x in res.columns]
    return(res)

def parse_log2(fn):
    """Parse out some information from a plink 2 log. This function currently
    only supports log files from linear or logistic regression.

    Parameters
    ----------
    fn : str 
        Path to the plink log file.

    Returns
    -------
    res : pandas.DataFrame
        Dataframe with log file information.

    """
    with open(fn) as f:
        lines = f.readlines()
    if len(lines) == 0:
        sys.stderr.write('Empty log file: {}.\n'.format(fn))
        return(None)
    logtype = None
    # TODO: Eventually, I will look for other arguments that indicate which
    # plink analysis was run.
    if len([x for x in lines if '--glm standard-beta' in x]):
        logtype = 'linear'
    elif len([x for x in lines if '--glm firth-fallback' in x]):
        logtype = 'logistic'
    if logtype is None:
        return(None)
        sys.stderr.write('Log file not supported: {}.\n'.format(fn))
    try:
        lines = [x for x in lines if 'remaining after' in x]
        i = 0
        x = lines[i].split()
        samples = int(x[0])
        females = int(x[2][1:])
        males = int(x[4])
        i += 1
        cases = np.nan
        controls = np.nan
        if logtype == 'logistic':
            x = lines[i].split()
            cases = int(x[0])
            controls = int(x[3])
            i += 1
        variants = int(lines[i].split()[0])
    except:
        sys.stderr.write('Error parsing log file: {}.\n'.format(fn))
        return(None)
    se = pd.Series([samples, females, males, cases, controls, variants],
                     index=['samples', 'females', 'males', 'cases',
                            'controls', 'variants']).dropna()
    return(se)
