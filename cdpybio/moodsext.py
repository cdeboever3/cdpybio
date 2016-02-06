import pandas as pd

def _filter_variant_motif_res(
    motif_res, 
    variant_start, 
    variant_end, 
    motif_length, 
    seq,
):
    """
    Remove MOODS motif hits that don't overlap the variant of interest.
    
    Parameters
    ----------
    motif_res : list
        Result from MOODS search like [(21, 3.947748787969809), (-38,
        3.979759977155675)].
        
    variant_start : int
        Relative start position of the allele string in seq (not genomic
        coordinates).  In other words, seq[variant_start:variant_end] should
        give the allele.
        
    variant_end : int
        Relative end position of the allele string in seq (not genomic
        coordinates).  In other words, seq[variant_start:variant_end] should
        give the allele.
        
    motif_length : int
        Length of the motif.
        
    seq : str
        Sequence searched for motifs for this allele.

    Returns
    -------
    motif_res : list
        List in the same format as motif_res input but with entries that don't
        overlap the variant removed. 

    """
    import MOODS
    remove = []
    for r in motif_res:
        motif_start = r[0]
        motif_end = r[0] + motif_length
        if r[0] < 0:
            motif_start += len(seq)
            motif_end += len(seq)

        if motif_end <= variant_start or motif_start >= variant_end:
            remove.append(r)
    motif_res = list(set(motif_res) - set(remove))
    return motif_res

def find_motif_disruptions(
    position, 
    ref, 
    alt, 
    genome_fasta, 
    matrices,
):
    """
    Determine whether there is a difference between the ref and alt
    alleles for TF binding. Requires samtools in your path.
    
    Parameters
    ----------
    position : str
        Zero based genomic coordinates of the reference allele of the form
        chrom:start-end (chr5:100-101 for a SNV for instance). The value end -
        start should equal the length of the ref allele.

    ref : str
        Reference allele. This should match the reference sequence at "position"
        in genome_fasta.

    alt : str
        Alternate allele.

    genome_fasta : str
        Path to genome fasta file. This file should be indexed.
    
    matrices : dict
        Dict whose keys are motif names and whose values are pandas data frames 
        or numpy arrays containing PWMs with columns ACGT.

    Returns
    -------
    out : pandas.DataFrame
        Pandas data frame with motifs whose best matches that overlapped the
        variant differed between the reference and alternate sequences. A score
        of zero and a strand of '' indicates that there was not a match for the
        motif on the given allele.

    """
    import subprocess
    import MOODS
    # import pybedtools as pbt
    max_motif_length = max([x.shape[0] for x in matrices.values()])
    chrom, coords = position.split(':')
    start,end = [int(x) for x in coords.split('-')]
    s = '{}:{}-{}'.format(chrom, start - max_motif_length + 1, end +
                          max_motif_length - 1)
    c = 'samtools faidx {} {}'.format(genome_fasta, s)
    seq_lines = subprocess.check_output(c, shell=True).strip().split()
    ref_seq = seq_lines[1]
    alt_seq = ref_seq[0:max_motif_length - 1] + alt + ref_seq[max_motif_length + len(ref) - 1:]

    ref_variant_start = max_motif_length - 1
    ref_variant_end = max_motif_length - 1 + len(ref)
    alt_variant_start = max_motif_length - 1
    alt_variant_end = max_motif_length - 1 + len(alt)

    ms = [matrices[x].T.values.tolist() for x in matrices.keys()]
    ref_res = MOODS.search(ref_seq, ms, 0.001, both_strands=True, 
                           bg=[0.25, 0.25, 0.25, 0.25])
    ref_res = dict(zip(matrices.keys(), ref_res))
    alt_res = MOODS.search(alt_seq, ms, 0.001, both_strands=True, 
                           bg=[0.25, 0.25, 0.25, 0.25])
    alt_res = dict(zip(matrices.keys(), alt_res))

    # First we'll remove any motif matches that don't overlap the variant of interest (and thus
    # can't be affected by the variant and will be the same for ref and alt). Then we'll get the 
    # best match for each motif for ref and alt.
    rows = []
    for motif in ref_res.keys():
        ref_res[motif] = _filter_variant_motif_res(ref_res[motif], ref_variant_start, ref_variant_end, 
                                           matrices[motif].shape[0], ref_seq)
        alt_res[motif] = _filter_variant_motif_res(alt_res[motif], alt_variant_start, alt_variant_end, 
                                           matrices[motif].shape[0], alt_seq)

        if len(ref_res[motif]) > 0:
            ref_pos, ref_score = sorted(ref_res[motif], key=lambda x: x[1], reverse=True)[0]
            ref_strand = {True:'+', False:'-'}[ref_pos > 0]
        else:
            ref_score = 0
            ref_strand = ''
        if len(alt_res[motif]) > 0:
            alt_pos, alt_score = sorted(alt_res[motif], key=lambda x: x[1], reverse=True)[0]
            alt_strand = {True:'+', False:'-'}[alt_pos > 0]
        else:
            alt_score = 0
            alt_strand = ''
        if ref_score > 0 or alt_score > 0:
            diff = ref_score - alt_score
            rows.append([motif, ref_score, ref_strand, alt_score, alt_strand, diff])
    out = pd.DataFrame(rows, columns=['motif', 'ref_score', 'ref_strand', 'alt_score', 
                                      'alt_strand', 'score_diff'])
    out.index = out.motif
    out = out.drop('motif', axis=1)
    out = out[out.score_diff != 0]
    return out
