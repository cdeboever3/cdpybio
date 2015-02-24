import os

import pandas as pd
import vcf

def record_variant_id(record):
    """Get variant ID from vcf.model._Record"""
    if record.ID:
        return record.ID
    else:
        return record.CHROM + ':' + str(record.POS)
    
def df_variant_id(row):
    """Get variant ID from vcf in DataFrame"""
    if row['ID'] != '.':
        return row['ID']
    else:
        return row['CHROM'] + ':' + str(row['POS'])
    
def bed_variant_id(interval):
    """Get variant ID from BedTool interval"""
    if interval.fields[2] != '.':
        return interval.fields[2]
    else:
        return bed_variant_position(interval)
    
def bed_variant_position(interval):
    """Get variant position from BedTool interval"""
    return interval.fields[0] + ':' + interval.fields[1]

def vcf_as_df(fn):
    """
    Read VCF file into pandas DataFrame.

    Parameters:
    -----------
    fn : str
        Path to VCF file.

    Returns
    -------
    df : pandas.DataFrame
        The VCF file as a data frame. Note that all header information is thrown
        away.

    """
    header_lines = 0
    with open(fn, 'r') as f:
        line = f.readline().strip()
        header_lines += 1
        while line[0] == '#':
            line = f.readline().strip()
            header_lines += 1
    
    header_lines -= 2
    df = pd.read_table(fn, skiprows=header_lines, header=0)
    df.columns = ['CHROM'] + list(df.columns[1:])
    return df

def make_het_matrix(fn):
    """
    Make boolean matrix of samples by variants. One indicates that the sample is
    heterozygous for that variant.

    Parameters:
    -----------
    vcf : str
        Path to VCF file.

    """
    # TODO: parallelize?
    vcf_df = vcf_as_df(fn)
    variant_ids = vcf_df.apply(lambda x: df_variant_id(x), axis=1)
    vcf_reader = vcf.Reader(open(fn, 'r'))
    record = vcf_reader.next()
    hets = pd.DataFrame(0, index=variant_ids,
                        columns=[x.sample for x in record.samples])
    
    vcf_reader = vcf.Reader(open(fn, 'r'))
    for record in vcf_reader:
        h = record.get_hets()
        i = record_variant_id(record)
        hets.ix[i, [x.sample for x in h]] = 1
    
    return hets
