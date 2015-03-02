import argparse
import os

import cdpybio as cpb

def main():
    parser = argparse.ArgumentParser(
        description=('Count number of reads supporting each base for different '
                     'genomic positions. You can provide a VCF file or a bed '
                     'file to specifiy the positions to count for.'))
    parser.add_argument('bam', help=('Coordinate sorted and indexed bam file '
                                     'to use for counting. Index is assumed to '
                                     'be named the same as the bam but with '
                                     '.bai appended to the end of the filename '
                                     '(e.g. XXX.bam and XXX.bam.bai.'))
    parser.add_argument('counts', help=('Output tsv file for nucleotide '
                                        'counts.'))
    parser.add_argument('--vcf', 
                        help='VCF file defining positions to be counted.')
    parser.add_argument('--bed', 
                        help='Bed file defining positions to be counted.')
    parser.add_argument('--stranded', 
                        help='Count nucleotides in a strand-specific manner.',
                        action='store_true')
    args = parser.parse_args()

    bam = args.bam
    counts = args.counts
    vcf = args.vcf
    bed = args.bed
    stranded = args.stranded

    # assert os.path.exists(bam + '.bai'), 'can\'t find bam index file'
    assert vcf or bed, 'must specify VCF or bed file'
    assert not (vcf and bed), 'cannot provide both VCF and bed files'

    if vcf:
        positions = vcf
    elif bed:
        positions = bed

    df = cpb.pysam.nt_counts(bam, positions, stranded=stranded, vcf=vcf,
                             bed=bed)
    df.to_csv(counts, sep='\t')

if __name__ == '__main__':
    main()
