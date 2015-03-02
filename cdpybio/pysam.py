import re

import pandas as pd
import pybedtools as pbt
import pysam

from general import parse_region

def get_region_nt_counts(region, bam, stranded=False):
    """
    Get counts of each nucleotide from a bam file for a given region. If R1 and
    R2 reads both overlap a position, only one count will be added. If the R1
    and R2 reads disagree at a position they both overlap, that read pair is not
    used for that position.  Can optionally output strand-specific counts.

    Parameters
    ----------
    region : str or list
        Region of type chrom:start-end, chrom:start-end:strand, or [chrom,
        start, end]. The strand is ignored for chrom:start-end:strand. For
        chrom:start-end, the coordinates are one-based inclusive. For example,
        the query chr1:10-11 will give you the counts for the 10th and 11th
        bases of chr1. For [chrom, start, end], the coordinates are zero-based
        and end exclusive (like a bed file). The query [chr1, 9, 11] will give
        you the coverage of the 10th and 11th bases of chr1. The region value is
        passed directly to pysam's pileup function. 

    bam : pysam.calignmentfile.AlignmentFile or str
        Bam file opened with pysam or path to bam file (must be sorted and
        indexed).

    stranded : boolean
        Boolean indicating whether read data is stranded and stranded nucleotide
        counts should be returned. Assumes R1 read on reverse strand implies +
        strand coverage etc.

    Returns
    -------
    counts : pandas.DataFrame
        Data frame with the counts for each base in the region. The index of
        this data frame is one-based for compatibility with VCF files.

    """
    # TODO: I should figure out what the different possible values are that
    # pysam could give me back (so far I only have ATCGN). Can I get deletions
    # and insertions? 
    # TODO: This could probably be parallelized.
    if type(bam) == str:
        bam = pysam.AlignmentFile(bam, 'rb')

    if type(region) is str:
        r = parse_region(region)
        if len(r) == 3:
            chrom, start, end = r
        elif len(r) == 4:
            chrom, start, end, strand = r
        start = int(start)
        end = int(end)
        ind = ['{}:{}'.format(chrom, x) for 
               x in range(start, end + 1)]
        pp = bam.pileup(region=region, truncate=True)
    
    elif type(region) is (list or tuple):
        chrom, start, end = region
        ind = ['{}:{}'.format(chrom, x) for 
               x in range(int(start) + 1, int(end) + 1)]
        pp = bam.pileup(chrom, start, end, truncate=True)
   
    cols = ['A', 'T', 'C', 'G', 'N']
    if stranded:
        cols = ['{}+'.format(x) for x in cols] + ['{}-'.format(x) for x in cols]
    counts = pd.DataFrame(0, index=ind, columns=cols)
    
    for pc in pp: 
        # Most of this code deals with R1 and R2 reads that overlap so that we
        # don't get two counts from one fragment.
        pos = pc.reference_pos + 1
        r1_qnames = []
        r1_nts = []
        r2_qnames = []
        r2_nts = []
        for pr in pc.pileups:
            qnames = [r1_qnames, r2_qnames][pr.alignment.is_read2]
            nts = [r1_nts, r2_nts][pr.alignment.is_read2]
            nt = _pos_nt(pr, pc.reference_pos, stranded)
            if nt:
                qnames.append(pr.alignment.qname)
                nts.append(nt)
        r1 = pd.Series(r1_nts, index=r1_qnames)
        r2 = pd.Series(r2_nts, index=r2_qnames)
        df = pd.DataFrame([r1, r2], index=['R1', 'R2']).T
        singles = df[df.isnull().sum(axis=1) == 1]
        doubles = df.dropna()
        vcs = []
        vcs.append(singles['R1'].value_counts())
        vcs.append(singles['R2'].value_counts())
        doubles = doubles[doubles.R1 == doubles.R2]
        vcs.append(doubles.R1.value_counts())
        for vc in vcs:
            counts.ix['{}:{}'.format(chrom, pos), vc.index] += vc
    return counts

def _pos_nt(pr, pos, stranded=False):
    """
    Given a pileup read and a position, return the base that is covered by the
    read at the given position if the position is covered.

    Parameters
    ----------
    pr : pysam.calignmentfile.PileupRead
        Region of type chrom:start-end, chrom:start-end:strand, or [chrom,
    
    pos : int
        Zero-based position of the nucleotide of interest in genomic
        coordinates.

    stranded : boolean
        Boolean indicating whether data is stranded and stranded nucleotide
        should be returned. Assumes R1 read on reverse strand implies + strand
        coverage etc.

    Returns
    -------
    nt : str or None
        If None, then the read did not cover the position. If not None, returns
        the nucleotide at that position (with + or - appended to indicate strand
        if desired).

    """
    nt = None
    bases = dict(zip(pr.alignment.get_reference_positions(), 
                     list(pr.alignment.seq.upper())))
    if pos in bases.keys():
        nt = bases[pos]
    if nt and stranded:
        strand = None
        if pr.alignment.is_read1 and pr.alignment.is_reverse:
            strand = '+'
        if pr.alignment.is_read2 and not pr.alignment.is_reverse:
            strand = '+'
        if pr.alignment.is_read1 and not pr.alignment.is_reverse:
            strand = '-'
        if pr.alignment.is_read2 and pr.alignment.is_reverse:
            strand = '-'
        nt = '{}{}'.format(nt, strand)
    return nt

def nt_counts(bam, positions, stranded=False, vcf=False, bed=False):
    """
    Find the number of nucleotides covered at all positions in a bed or vcf
    file.

    Parameters
    ----------
    bam : str or pysam.calignmentfile.AlignmentFile 
        Bam file opened with pysam or path to bam file (must
        be sorted and indexed).
    
    positions : str or pybedtools.BedTool
        Path to bed or vcf file or pybedtools.BedTool object. The extension is
        used to determine whether the file is a bed or vcf (.bed vs .vcf).

    stranded : boolean
        Boolean indicating whether read data is stranded and stranded nucleotide
        counts should be returned. Assumes R1 read on reverse strand implies +
        strand coverage etc.

    vcf : boolean
        Set to True if you are providing a vcf file that doesn't have a .vcf
        suffix.

    bed : boolean
        Set to True if you are providing a bed file that doesn't have a .bed
        suffix.

    Returns
    -------
    counts : pandas.DataFrame
        Data frame with the counts for each base in the region. The index of
        this data frame is one-based for compatibility with VCF files.

    """
    if not bed and not vcf:
        if type(positions) == pbt.bedtool.BedTool:
            df = positions.to_dataframe()
        assert type(positions) is str, ('positions must be BedTool, bed file, '
                                        'or vcf file')
        if positions[-4:] == '.bed':
            bed = True
        elif positions[-4:] == '.vcf':
            vcf = True

    if bed:
        df = pbt.BedTool(positions).to_dataframe()
    elif vcf:
        from variants import vcf_as_df
        tdf = vcf_as_df(positions)
        df = pd.DataFrame(index=tdf.index)
        df['chrom'] = tdf.CHROM
        df['start'] = tdf.POS - 1
        df['end'] = tdf.POS

    res = []
    for i in df.index:
        region = [df.ix[i, 'chrom'], df.ix[i, 'start'], df.ix[i, 'end']]
        res.append(get_region_nt_counts(region, bam, stranded))
    res = pd.concat(res)
    return res
