import re

import numpy as np
import pandas as pd
import scipy.stats as stats

R_REGEX = re.compile('(.*):(.*)-(.*)')
R_REGEX_STRAND = re.compile('(.*):(.*)-(.*):(.*)')

def estimate_allele_frequency(ac, an, a=1, b=100):
    """
    Make sample (or other) names.

    Parameters:
    -----------

    ac : array-like
        Array-like object with the observed allele counts for each variant. If
        ac is a pandas Series, the output dataframe will have the same index as
        ac.

    an : array-like
        Array-like object with the number of haplotypes that were genotyped.

    a : float
        Parameter for prior distribution beta(a, b).

    b : float
        Parameter for prior distribution beta(a, b).

    Returns
    -------
    out : pandas.DataFrame
        Pandas dataframe with allele frequency estimate

    """
    # Credible interval is 95% highest posterior density
    td = dict(zip(['ci_lower', 'ci_upper'], 
                  stats.beta(a + ac, b + an - ac).interval(0.95)))
    td['af'] = (a + ac) / (a + b + an)
    td['af_mle'] = np.array(ac).astype(float) / np.array(an)
    out = pd.DataFrame(td)[['af_mle', 'af', 'ci_lower', 'ci_upper']]
    if type(ac) == pd.Series:
        out.index = ac.index
    return(out)

def transform_standard_normal(df):
    """Transform a series or the rows of a dataframe to the values of a standard
    normal based on rank."""
    import pandas as pd
    import scipy.stats as stats
    if type(df) == pd.core.frame.DataFrame:
        gc_ranks = df.rank(axis=1)
        gc_ranks = gc_ranks / (gc_ranks.shape[1] + 1)
        std_norm = stats.norm.ppf(gc_ranks)
        std_norm = pd.DataFrame(std_norm, index=gc_ranks.index, 
                                columns=gc_ranks.columns)
    elif type(df) == pd.core.series.Series:
        gc_ranks = df.rank()
        gc_ranks = gc_ranks / (gc_ranks.shape[0] + 1)
        std_norm = stats.norm.ppf(gc_ranks)
        std_norm = pd.Series(std_norm, index=df.index)
    return std_norm

def read_gzipped_text_url(url):
    """Read a gzipped text file from a URL and return 
    contents as a string."""
    import urllib2
    import zlib
    from StringIO import StringIO

    opener = urllib2.build_opener() 
    request = urllib2.Request(url)
    request.add_header('Accept-encoding', 'gzip')
    respond = opener.open(request)
    compressedData = respond.read()
    respond.close()
    opener.close()
    compressedDataBuf = StringIO(compressedData)
    d = zlib.decompressobj(16+zlib.MAX_WBITS)
    buffer = compressedDataBuf.read(1024)
    #saveFile = open('/tmp/test.txt', "wb")
    s = []
    while buffer:
        s.append(d.decompress(buffer))
        buffer = compressedDataBuf.read(1024)
    s = ''.join(s)
    return s

def parse_region(region):
    """
    Parse region of type chr1:10-20 or chr1:10-20:+

    Parameters:
    -----------

    region : str
        Region of type chr1:10-20 or chr1:10-20:+.

    Returns
    -------
    groups : tuple
        Tuple of groups from regex e.g. (chr1, 10, 20) or (chr1, 10, 20, +).

    """
    m = R_REGEX_STRAND.search(region)
    if not m:
        m = R_REGEX.search(region)
    if m:
        groups = m.groups()
        return groups
    else:
        return None

def _sample_names(files, kwargs):
    """
    Make sample (or other) names.

    Parameters:
    -----------

    files : list of string
        Typically a list of file paths although could be any list of strings
        that you want to make names for. If neither names nor define_sample_name
        are provided, then files is returned as is.

    kwargs : dict
        kwargs from another function. Can include the following keys with
        appropriate arguments.

    names : list of strings
        Names to use. Overrides define_sample_name if provided.

    define_sample_name : function that takes string as input
        Function mapping string to name. For instance, you may have a sample
        name in a file path and use a regex to extract it.

    """
    if 'define_sample_name' not in kwargs.keys():
        define_sample_name = lambda x: x
    else:
        define_sample_name = kwargs['define_sample_name']
    
    if 'names' in kwargs.keys():
        names = kwargs['names']
    else:
        names = [define_sample_name(f) for f in files]
    
    assert len(names) == len(files)
    return names
