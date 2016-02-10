import copy
import pandas as pd
import pybedtools as pbt

from general import _sample_names

class AnnotatedBed:
    def __init__(
        self,   
        bed,
        annot_beds,
        completely_contains=None,
    ):  
        """
        Initialize AnnotatedBed object.
            
        Parameters
        ----------
        bed : str or pybedtools.Bedtool
            Bed file to annotate.
        
        annot_beds : dict
            Dict whose keys are names (like 'gene', 'promoter', etc.) and whose
            values are bed files to annotate the input bed file with.
            
        completely_contains : list
            List of keys from annot_beds. For these beds, we will check whether
            the features is entirely contained by features "bed."
        
        """            
        self._initialize_bt(bed)
        self._num_cols = len(self.bt[0].fields)
        self._has_name_col = self._num_cols > 3
        self._initialize_dfs()
        self._initialize_annot_beds(annot_beds)
        for k in annot_beds.keys():
            self.annotate_bed(self.bt, k, k)
        for k in completely_contains:
            self.annotate_bed(self.bt, k, k, complete=True)
        # Remove bt object because pybedtools relies on the file on the disk and
        # we don't know if the file will always be around.
        self.bt = None
        
    def _initialize_annot_beds(
        self, 
        annot_beds,
    ):
        import pybedtools as pbt
        self.annot_beds = dict()
        for k in annot_beds.keys():
            if type(annot_beds[k]) == str:
                self.annot_beds[k] = pbt.BedTool(annot_beds[k])
            else:
                self.annot_beds[k] = annot_beds[k]

    def _initialize_dfs(
        self,
    ):  
        self.df = self.bt.to_dataframe()
        if self._has_name_col:
            if len(set(self.df.name)) != self.df.shape[0]:
                self._has_name_col = False
        if self._has_name_col:
            self.df.index = self.df.name
        else:
            self.df.index = (self.df.chrom.astype + ':' +
                             self.df.start.astype(str) + '-' +
                             self.df.end.astype(str))
        self.feature_to_df = pd.DataFrame(index=self.df.index)

    def _initialize_bt(
        self,
        bed,           
    ):                 
        import pybedtools as pbt
        if type(bed) == str:
            self.bt = pbt.BedTool(bed)
        else:
            self.bt = bed
        self.bt = self.bt.sort()
        
    def bt_from_df(self):
        """Make a BedTool object for the input bed file."""
        import pybedtools as pbt
        s = ('\n'.join(df.astype(str).apply(lambda x: '\t'.join(x), axis=1)) +
             '\n')
        df.bt = pbt.BedTool(s, from_string=True)
        
    def annotate_bed(
        self,
        bt,
        name,
        col_name,
        complete=None,
    ):
        """name is the key for the bed file in self.annot_beds. col_name is used
        to name the columns created in the dataframe."""
        import numpy as np
        import pandas as pd
        has_name_col = len(self.annot_beds[name][0].fields) > 3
        if complete:
            res = bt.intersect(self.annot_beds[name], sorted=True, wo=True, F=1)
            col_name = col_name + '_complete'
        else:
            res = bt.intersect(self.annot_beds[name], sorted=True, wo=True)
        df = res.to_dataframe(names=range(len(res[0].fields)))
        if self._has_name_col:
            ind = df[3].values
        else:
            ind = list(df[0].astype(str) + ':' +
                       df[1].astype(str) + '-' +
                       df[2].astype(str))
        if has_name_col:
            vals = df[self._num_cols + 3].values
        else:
            vals = list(df[self._num_cols + 0].astype(str) + ':' +
                        df[self._num_cols + 1].astype(str) + '-' +
                        df[self._num_cols + 2].astype(str))
        self.df[col_name] = False
        self.df.ix[set(ind), col_name] = True
        se = pd.Series(vals, index=ind)
        vc = pd.Series(se.index).value_counts()
        self.feature_to_df[col_name] = np.nan
        self.feature_to_df.ix[list(vc[vc == 1].index), col_name] \
                = se[list(vc[vc == 1].index)].apply(lambda x: set([x]))
        m = list(set(vc[vc > 1].index))
        v = []
        for i in m:
            v.append(set(se[i].values))
        self.feature_to_df.ix[m, col_name] = v
        
class AnnotatedInteractions(AnnotatedBed):
    def __init__(
        self,   
        df,
        annot_beds,
        completely_contains=None,
    ):  
        """
        Initialize AnnotatedBed object.
            
        Parameters
        ----------
        df : pandas.DataFrame
            Dataframe with peaks. Must contain columns chrom1, start1, end1,
            chrom2, start2, end2 and have unique index.
        
        annot_beds : dict
            Dict whose keys are names (like 'gene', 'promoter', etc.) and whose
            values are bed files to annotate the input bed file with.
        
        """            
        self.df = df.copy(deep=True)
        self.df['name'] = self.df.index
        self.feature_to_df = pd.DataFrame(index=self.df.index)
        self.annotate_interactions()
        self.bt_from_df()
        self._initialize_annot_beds(annot_beds)
        suffixes = ['1', '2', '_loop', '_loop_inner']
        for i,bt in enumerate([self.bt1, self.bt2, self.bt_loop,
                               self.bt_loop_inner]):
            for k in annot_beds.keys():
                self.annotate_bed(bt, k, '{}{}'.format(k, suffixes[i]))
        for i,bt in enumerate([self.bt1, self.bt2, self.bt_loop,
                               self.bt_loop_inner]):
            for k in completely_contains:
                self.annotate_bed(bt, k, '{}{}'.format(k, suffixes[i]),
                                  complete=True)
        # Remove BedTools objects because pybedtools relies on the file on the
        # disk and we don't know if the file will always be around.
        self.bt1 = None
        self.bt2 = None
        self.bt_loop = None
        self.bt_loop_inner = None
        
    def annotate_bed(
        self,
        bt,
        name,
        col_name,
        complete=None,
    ):
        """name is the key for the bed file in self.annot_beds. col_name is used
        to name the columns created in the dataframe."""
        import numpy as np
        import pandas as pd
        has_name_col = len(self.annot_beds[name][0].fields) > 3
        if complete:
            res = bt.intersect(self.annot_beds[name], sorted=True, wo=True, F=1)
            col_name = col_name + '_complete'
        else:
            res = bt.intersect(self.annot_beds[name], sorted=True, wo=True)
        try:
            df = res.to_dataframe(names=range(len(res[0].fields)))
            ind = df[3].values
            if has_name_col:
                vals = df[7].values
            else:
                vals = list(df[4].astype(str) + ':' +
                            df[5].astype(str) + '-' +
                            df[6].astype(str))
            self.df[col_name] = False
            self.df.ix[set(ind), col_name] = True
            se = pd.Series(vals, index=ind)
            vc = pd.Series(se.index).value_counts()
            self.feature_to_df[col_name] = np.nan
            self.feature_to_df.ix[list(vc[vc == 1].index), col_name] \
                    = se[list(vc[vc == 1].index)].apply(lambda x: set([x]))
            m = list(set(vc[vc > 1].index))
            v = []
            for i in m:
                v.append(set(se[i].values))
            self.feature_to_df.ix[m, col_name] = v
        except IndexError:
            pass
        
    def annotate_interactions(self):
        import numpy as np
        self.df['anchor1'] = (self.df.chrom1.astype(str) + ':' +
                              self.df.start1.astype(str) + '-' +
                              self.df.end1.astype(str))
        self.df['anchor2'] = (self.df.chrom2.astype(str) + ':' +
                              self.df.start2.astype(str) + '-' +
                              self.df.end2.astype(str))
        self.df['intra'] = True
        self.df.ix[self.df.chrom1 != self.df.chrom2, 'intra'] = False
        ind = self.df[self.df.intra].index
        self.df['loop'] = np.nan
        self.df.ix[ind, 'loop'] = (
            self.df.chrom1 + ':' + 
            self.df[['start1', 'start2']].min(axis=1).astype(str) + '-' + 
            self.df[['end1', 'end2']].max(axis=1).astype(str))
        self.df['loop_inner'] = np.nan
        self.df.ix[ind, 'loop_inner'] = (
            self.df.chrom1 + ':' + 
            self.df[['end1', 'end2']].min(axis=1).astype(str) + '-' + 
            self.df[['start1', 'start2']].max(axis=1).astype(str))
        self.df['loop_length'] = (self.df[['end1', 'end2']].max(axis=1) - 
                                  self.df[['start1', 'start2']].min(axis=1))
        self.df['loop_inner_length'] = (
            self.df[['start1', 'start2']].max(axis=1) - 
            self.df[['end1', 'end2']].min(axis=1))
        
    def bt_from_df(self):         
        import pybedtools as pbt
        s = '\n'.join(self.df.chrom1.astype(str) + '\t' +
                      self.df.start1.astype(str) + '\t' +
                      self.df.end1.astype(str) + '\t' + self.df.name) + '\n'
        self.bt1 = pbt.BedTool(s, from_string=True).sort()
        s = '\n'.join(self.df.chrom2.astype(str) + '\t' +
                      self.df.start2.astype(str) + '\t' +
                      self.df.end2.astype(str) + '\t' + self.df.name) + '\n'
        self.bt2 = pbt.BedTool(s, from_string=True).sort()
        s = '\n'.join(self.df.chrom1.astype(str) + '\t' + 
                      self.df[['start1', 'start2']].min(axis=1).astype(str) +
                      '\t' + self.df[['end1', 'end2']].max(axis=1).astype(str) +
                      '\t' + self.df.name) + '\n'
        self.bt_loop = pbt.BedTool(s, from_string=True).sort()
        s = '\n'.join(self.df.chrom1.astype(str) + '\t' + 
                      self.df[['end1', 'end2']].min(axis=1).astype(str) + '\t' +
                      self.df[['start1', 'start2']].max(axis=1).astype(str)  +
                      '\t' + self.df.name) + '\n'
        self.bt_loop_inner = pbt.BedTool(s, from_string=True).sort()

def beds_to_boolean(beds, ref=None, beds_sorted=False, ref_sorted=False,
                    **kwargs):
    """
    Compare a list of bed files or BedTool objects to a reference bed file and
    create a boolean matrix where each row is an interval and each column is a 1
    if that file has an interval that overlaps the row interval and a 0
    otherwise. If no reference bed is provided, the provided bed files will be
    merged into a single bed and compared to that.

    Parameters
    ----------
    beds : list
        List of paths to bed files or BedTool objects.

    ref : str or BedTool
        Reference bed file to compare against. If no reference bed is provided,
        the provided bed files will be merged into a single bed and compared to
        that.

    beds_sorted : boolean
        Whether the bed files in beds are already sorted. If False, all bed
        files in beds will be sorted.

    ref_sorted : boolean
        Whether the reference bed file is sorted. If False, ref will be sorted.

    names : list of strings
        Names to use for columns of output files. Overrides define_sample_name 
        if provided.

    define_sample_name : function that takes string as input
        Function mapping filename to sample name (or basename). For instance,
        you may have the basename in the path and use a regex to extract it.
        The basenames will be used as the column names. If this is not provided,
        the columns will be named as the input files.

    Returns
    -------
    out : pandas.DataFrame
        Boolean data frame indicating whether each bed file has an interval
        that overlaps each interval in the reference bed file. 

    """
    beds = copy.deepcopy(beds)
    fns = []
    for i,v in enumerate(beds):
        if type(v) == str:
            fns.append(v)
            beds[i] = pbt.BedTool(v)
        else:
            fns.append(v.fn)
        if not beds_sorted:
            beds[i] = beds[i].sort()

    names = _sample_names(fns, kwargs)
    if ref:
        if type(ref) == str:
            ref = pbt.BedTool(ref)
        if not ref_sorted:
            ref = ref.sort()
    else:
        ref = combine(beds)
    
    ind = []
    for r in ref:
        ind.append('{}:{}-{}'.format(r.chrom, r.start, r.stop))
    bdf = pd.DataFrame(0, index=ind, columns=names)
    for i,bed in enumerate(beds):
        res = ref.intersect(bed, sorted=True, wa=True)
        ind = []
        for r in res:
            ind.append('{}:{}-{}'.format(r.chrom,
                                         r.start,
                                         r.stop))
        bdf.ix[ind, names[i]] = 1
    return bdf

def combine(beds, beds_sorted=False, postmerge=True):
    """
    Combine a list of bed files or BedTool objects into a single BedTool object.

    Parameters
    ----------
    beds : list
        List of paths to bed files or BedTool objects.

    beds_sorted : boolean
        Whether the bed files in beds are already sorted. If False, all bed
        files in beds will be sorted.

    postmerge : boolean
        Whether to merge intervals after combining beds together. 

    Returns
    -------
    out : pybedtools.BedTool
        New sorted BedTool with intervals from all input beds.

    """
    beds = copy.deepcopy(beds)
    for i,v in enumerate(beds):
        if type(v) == str:
            beds[i] = pbt.BedTool(v)
        if not beds_sorted:
            beds[i] = beds[i].sort()

    # For some reason, doing the merging in the reduce statement doesn't work. I
    # think this might be a pybedtools bug. In any fashion, I can merge
    # afterward although I think it makes a performance hit because the combined
    # bed file grows larger than it needs to.
    out = reduce(lambda x,y : x.cat(y, postmerge=False), beds)
    out = out.sort()
    if postmerge:
        out = out.merge()
    return out


def write_bed_with_trackline(bed, out, trackline, add_chr=False):
    """
    Read a bed file and write a copy with a trackline. Here's a simple trackline
    example: 'track type=bed name="cool" description="A cool track."'

    Parameters
    ----------
    bed : str 
        Input bed file name.
    out : str
        Output bed file name.
    trackline : str
        UCSC trackline.
    add_chr : boolean
        Add 'chr' to the chromosomes in the input file. Necessary for
        UCSC genome browser if not present.

    """
    df = pd.read_table(bed, index_col=None, header=None)
    bt = pbt.BedTool('\n'.join(df.apply(lambda x: '\t'.join(x.astype(str)), 
                                        axis=1)) + '\n',
                     from_string=True)
    if add_chr:
        bt = add_chr_to_contig(bt)
    bt = bt.saveas(out, trackline=trackline)

def strip_chr(bt):
    """Strip 'chr' from chromosomes for BedTool object

    Parameters
    ----------
    bt : pybedtools.BedTool
        BedTool to strip 'chr' from.

    Returns
    -------
    out : pybedtools.BedTool
        New BedTool with 'chr' stripped from chromosome names.

    """
    try:
        df = pd.read_table(bt.fn, header=None, dtype=str)
    # If the try fails, I assume that's because the file has a trackline. Note
    # that I don't preserve the trackline (I'm not sure how pybedtools keeps
    # track of it anyway).
    except pd.parser.CParserError:
        df = pd.read_table(bt.fn, header=None, skiprows=1, dtype=str)
    df[0] = df[0].apply(lambda x: x[3:])
    s = '\n'.join(df.astype(str).apply(lambda x: '\t'.join(x), axis=1)) + '\n'
    out = pbt.BedTool(s, from_string=True)
    return out

def add_chr(bt):
    """Add 'chr' to chromosomes for BedTool object

    Parameters
    ----------
    bt : pybedtools.BedTool
        BedTool to add 'chr' to.

    Returns
    -------
    out : pybedtools.BedTool
        New BedTool with 'chr' added to chromosome names.

    """
    try:
        df = pd.read_table(bt.fn, header=None, dtype=str)
    # If the try fails, I assume that's because the file has a trackline. Note
    # that I don't preserve the trackline (I'm not sure how pybedtools keeps
    # track of it anyway).
    except pd.parser.CParserError:
        df = pd.read_table(bt.fn, header=None, skiprows=1, dtype=str)
    df[0] = 'chr' + df[0]
    s = '\n'.join(df.astype(str).apply(lambda x: '\t'.join(x), axis=1)) + '\n'
    out = pbt.BedTool(s, from_string=True)
    return out

def intervals_to_bed(intervals):
    """
    Convert list of intervals of format chr1:100-200 or chr1:100-200:+ to 
    BedTool object.

    Parameters
    ----------
    intervals : array-like
        List of intervals.

    Returns
    -------
    bt : pybedtools.BedTool
        BedTool with one line for each interval.

    """
    import re
    strand = re.compile('(.*):(.*)-(.*):(\+|-)')
    no_strand = re.compile('(.*):(.*)-(.*)')
    bed_lines = []
    s = False
    for i in intervals:
        m = strand.match(i)
        if m:
            bed_lines.append('\t'.join([m.group(x) for x in range(1, 5)]))
        else:
            m = no_strand.match(i)
            if m:
                bed_lines.append('\t'.join([m.group(x) for x in range(1, 4)]))
    bt = pbt.BedTool('\n'.join(bed_lines) + '\n', from_string=True)
    return bt
