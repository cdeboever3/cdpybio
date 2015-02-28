import re

import pandas
import pysam

from general import parse_region

def get_region_counts(region, bam, stranded=False):
    """
    Get counts of each nucleotide from a bam file for a given region.

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
        Boolean indicating whether data is stranded and stranded counts should 
        be returned. Assumes R1 read on reverse strand implies + strand coverage
        etc.

    Returns
    -------
    out : pandas.DataFrame
        Data frame with the counts for each base in the region. The index of
        this data frame is one-based for compatibility with VCF files.

    """
    # TODO: I should figure out what the different values are that are possible
    # to get back (so far I only have ATCGN). Can I get deletions and
    # insertions? 
    # TODO: Strand-specific counting, don't double count for overlapping reads.
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
        pos = pc.reference_pos + 1
        for pr in pc.pileups:
            nt = _pos_nt(pr, pc.reference_pos, stranded)
            if nt:
                counts.ix['{}:{}'.format(chrom, pos), nt] += 1
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

def get_counts(bam, vcf, bed):
    # I'll probably need to just make a few big lists and make them into a
    # reasonable format at the end and save to a file. 
    
    
    #base_counts = pd.DataFrame(index=sample_hets.index, 
    #                           columns=['A', 'T', 'C', 'G'])

    if vcf:

    # iterate through regions
    # If POS is the position in the VCF file, I need to use samfile.pileup on
    # chr1:POS - 1:POS - 1


for s in allele_counts.index:
    allele_counts.ix[s] = pd.Series(nd)

