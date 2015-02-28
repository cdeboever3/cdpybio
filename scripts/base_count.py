import argparse
import os

import pysam

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
    parser.add_argument('--vcf', 
                        help='VCF file defining positions to be counted.')
    parser.add_argument('--bed', 
                        help='Bed file defining positions to be counted.')
    args = parser.parse_args()

    bam = args.bam
    vcf = args.vcf
    bed = args.bed

    # assert os.path.exists(bam + '.bai'), 'can\'t find bam index file'
    assert vcf or bed, 'must specify VCF or bed file'
    assert not (vcf and bed), 'cannot provide both VCF and bed files'

if __name__ == '__main__':
    main()
