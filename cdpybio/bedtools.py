import copy
import pandas as pd
import pybedtools as pbt

from general import _sample_names

class AnnotatedInteractions:
    def __init__(
        self,   
        df,
        annot_beds,
        completely_contains=None,
    ):  
        """
        Initialize AnnotatedInteractions object.
            
        Parameters
        ----------
        df : pandas.DataFrame
            Dataframe with peaks. Must contain columns chrom1, start1, end1,
            chrom2, start2, and end2. Other columns will not be removed but 
            may be overwritten if they clash with column names created here.
            Interactions must be unique.
        
        annot_beds : dict
            Dict whose keys are names (like 'gene', 'promoter', etc.) and whose
            values are bed files to annotate the input bed file with.
        
        """            
        self.df = df.copy(deep=True)
        self.df.index = (self.df.chrom1.astype(str) + ':' +
                         self.df.start1.astype(str) + '-' +
                         self.df.end1.astype(str) + '==' +
                         self.df.chrom2.astype(str) + ':' +
                         self.df.start1.astype(str) + '-' +
                         self.df.end2.astype(str))
        assert len(set(self.df.index)) == self.df.shape[0]
        self.df['name'] = self.df.index
        self.feature_to_df = pd.DataFrame(index=self.df.index)
        self.annotate_interactions()
        self.bts_from_df()
        self._initialize_annot_beds(annot_beds)
        for k in annot_beds.keys():
            self.annotate_bed(bt=self.bt1, name=k, col_name='{}1'.format(k),
                              df_col='anchor1')
            if k in completely_contains:
                self.annotate_bed(bt=self.bt1, name=k,
                                  col_name='{}1_complete'.format(k),
                                  df_col='anchor1', complete=True)
        for k in annot_beds.keys():
            self.annotate_bed(bt=self.bt2, name=k, col_name='{}2'.format(k),
                              df_col='anchor2')
            if k in completely_contains:
                self.annotate_bed(bt=self.bt2, name=k,
                                  col_name='{}2_complete'.format(k),
                                  df_col='anchor2', complete=True)
        for k in annot_beds.keys():
            self.annotate_bed(bt=self.bt1, name=k, col_name='{}_loop'.format(k),
                              df_col='loop')
            if k in completely_contains:
                self.annotate_bed(bt=self.bt1, name=k,
                                  col_name='{}_loop_complete'.format(k),
                                  df_col='loop', complete=True)
        for k in annot_beds.keys():
            self.annotate_bed(bt=self.bt1, name=k,
                              col_name='{}_loop_inner'.format(k),
                              df_col='loop_inner')
            if k in completely_contains:
                self.annotate_bed(bt=self.bt1, name=k,
                                  col_name='{}_loop_inner_complete'.format(k), 
                                  df_col='loop_inner', complete=True)
        self._bt1_path = None
        self._bt2_path = None
        self._bt_loop_path = None
        self._bt_loop_inner_path = None

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

    def load_saved_bts(self):
        """If the AnnotatedInteractions object was saved to a pickle and
        reloaded, this method remakes the BedTool objects."""
        if self._bt1_path:
            self.bt1 = pbt.BedTool(self._bt1_path)
        if self._bt2_path:
            self.bt2 = pbt.BedTool(self._bt2_path)
        if self._bt_loop_path:
            self.bt_loop = pbt.BedTool(self._bt_loop_path)
        if self._bt_loop_inner_path:
            self.bt_loop_inner = pbt.BedTool(self._bt_loop_inner_path)
        
    def save(
        self,
        path,
        name,
    ):
        """
        Save AnnotatedInteractions object and bed files. The object is stored in
        a pickle and the bed files are saved as separate bed files. The object
        can be reloaded by reading the pickle using cPickle and the BedTool
        objects can be recreated using .load_saved_bts().
            
        Parameters
        ----------
        path : str
            Path to save files to. Path should include a basename for the files.
            For instance, path='~/abc' will create files like ~/abc.pickle,
            ~/abc_anchor1.bed, etc.
            
        name : str
            Descriptive name used for bed file trackline.
        """
        t = 'track type=bed name="{}_anchor1"'.format(name)
        self.bt1.saveas(path + '_anchor1.bed', trackline=t)
        self._bt1_path = path + '_anchor1.bed'
        t = 'track type=bed name="{}_anchor2"'.format(name)
        self.bt2.saveas(path + '_anchor2.bed', trackline=t)
        self._bt2_path = path + '_anchor2.bed'
        t = 'track type=bed name="{}_loop"'.format(name)
        self.bt_loop.saveas(path + '_loop.bed', trackline=t)
        self._bt_loop_path = path + '_loop.bed'
        t = 'track type=bed name="{}_loop_inner"'.format(name)
        self.bt_loop_inner.saveas(path + '_loop_inner.bed', trackline=t)
        self._bt_loop_inner_path = path + '_loop_inner.bed'
        import cPickle
        cPickle.dump(self, open(path + '.pickle', 'w'))   
    
    def annotate_bed(
        self,
        bt,
        name,
        col_name,
        complete=None,
        df_col=None,
    ):
        """
        Annotate the input bed file using one of the annotation beds.
        Parameters
        ----------
        bt : pybedtools.BedTool
            BedTool for either one of the anchors, the loops,
            or the loop inners.
        
        name : str
            The key for the annoation bed file in annot_beds. 
        
        col_name : str
            Used to name the columns that will be made.
            
        complete : bool
            If True, this method will check whether the features in the
            annotation bed are completely contained by the features in the input
            bed.
            
        df_col : str
            If the name for bt isn't the index of self.df, this specifies
            which column of self.df contains the names for bt. For instance,
            if bt is the anchor1 BedTool, the df_col='anchor11'.
        
        """
        import numpy as np
        import pandas as pd
        has_name_col = len(self.annot_beds[name][0].fields) > 3
        print('one')
        if complete:
            res = bt.intersect(self.annot_beds[name], sorted=True, wo=True, F=1)
        else:
            res = bt.intersect(self.annot_beds[name], sorted=True, wo=True)
        print('two')
        try:
            df = res.to_dataframe(names=range(len(res[0].fields)))
            ind = df[3].values
            if df_col is None:
                self.df[col_name] = False
                self.df.ix[set(ind), col_name] = True
            else:
                tdf = pd.DataFrame(True, index=ind, columns=[col_name])
                self.df = self.df.merge(tdf, left_on=df_col, right_index=True,
                                        how='outer')
                self.df[col_name] = self.df[col_name].fillna(False)
                #self.df.ix[self.df[col_name].isnull(), col_name] = False
            print('a')
            if has_name_col:
                vals = df[7].values
            else:
                vals = list(df[4].astype(str) + ':' +
                            df[5].astype(str) + '-' +
                            df[6].astype(str))
            print('b')
            df.index = vals
            gb = df.groupby(3)
            t = pd.Series(gb.groups)
            print('c')
            t = pd.DataFrame(t.apply(lambda x: set(x)))
            print('d')
            t.columns = ['{}_features'.format(col_name)]
            self.df = self.df.merge(t, left_on=df_col, right_index=True,
                                    how='outer')
            print('e')
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
            self.df.ix[ind, 'chrom1'] + ':' + 
            self.df.ix[ind, ['start1', 'start2']].min(axis=1).astype(str) + 
            '-' + self.df.ix[ind, ['end1', 'end2']].max(axis=1).astype(str))
        self.df['loop_length'] = (self.df[['end1', 'end2']].max(axis=1) - 
                                  self.df[['start1', 'start2']].min(axis=1))
        ind = ind[(self.df.ix[ind, ['start1', 'start2']].max(axis=1) >
                   self.df.ix[ind, ['end1', 'end2']].min(axis=1))]
        self.df['loop_inner'] = np.nan
        self.df.ix[ind, 'loop_inner'] = (
            self.df.ix[ind, 'chrom1'] + ':' + 
            self.df.ix[ind, ['end1', 'end2']].min(axis=1).astype(str) + '-' +
            self.df.ix[ind, ['start1', 'start2']].max(axis=1).astype(str))
        self.df['loop_inner_length'] = (
            self.df[['start1', 'start2']].max(axis=1) - 
            self.df[['end1', 'end2']].min(axis=1))
        
    def bts_from_df(self):         
        import pybedtools as pbt
        s = '\n'.join(list(set(
            self.df.chrom1.astype(str) + '\t' + self.df.start1.astype(str) +
            '\t' + self.df.end1.astype(str) + '\t' + self.df.chrom1.astype(str)
            + ':' + self.df.start1.astype(str) + '-' +
            self.df.end1.astype(str)))) + '\n'
        self.bt1 = pbt.BedTool(s, from_string=True).sort()
        s = '\n'.join(list(set(
            self.df.chrom2.astype(str) + '\t' + self.df.start2.astype(str) +
            '\t' + self.df.end2.astype(str) + '\t' + self.df.chrom2.astype(str)
            + ':' + self.df.start2.astype(str) + '-' +
            self.df.end2.astype(str)))) + '\n'
        self.bt2 = pbt.BedTool(s, from_string=True).sort()
        ind = self.df[self.df.intra].index
        s = '\n'.join(
            self.df.ix[ind, 'chrom1'].astype(str) + '\t' + 
            self.df.ix[ind, ['start1', 'start2']].min(axis=1).astype(str) + 
            '\t' + self.df.ix[ind, ['end1', 'end2']].max(axis=1).astype(str) +
            '\t' + self.df.ix[ind, 'name']) + '\n'
        self.bt_loop = pbt.BedTool(s, from_string=True).sort()
        ind = ind[(self.df.ix[ind, ['start1', 'start2']].max(axis=1) >
                   self.df.ix[ind, ['end1', 'end2']].min(axis=1))]
        s = '\n'.join(
            self.df.ix[ind, 'chrom1'].astype(str) + '\t' + 
            self.df.ix[ind, ['end1', 'end2']].min(axis=1).astype(str) + '\t' +
            self.df.ix[ind, ['start1', 'start2']].max(axis=1).astype(str)  +
            '\t' + self.df.ix[ind, 'name']) + '\n'
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
