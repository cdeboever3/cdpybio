import pdb
import sys

import pandas as pd

def _gencode_donor(chrom, start, end, strand):
    if strand == '+':
        return '{}:{}:{}'.format(chrom, start, strand)
    if strand == '-':
        return '{}:{}:{}'.format(chrom, end, strand)

def _gencode_acceptor(chrom, start, end, strand):
    if strand == '+':
        return '{}:{}:{}'.format(chrom, end, strand)
    if strand == '-':
        return '{}:{}:{}'.format(chrom, start, strand)

def make_transcript_gene_se(fn):
    """
    Make a Pandas Series with transcript ID's as the index and values as the
    gene ID containing that transcript.

    Parameters
    ----------
    fn : str of filename 
        Filename of the Gencode gtf file

    Returns
    -------
    se : pandas.Series
        Make a Pandas Series with transcript ID's as the index and values as the
        gene ID containing that transcript.

    """
    import itertools as it

    import HTSeq

    gtf = it.islice(HTSeq.GFF_Reader(fn), None)
    transcripts = []
    genes = []
    line = gtf.next()
    while line != '':
        if line.type == 'transcript':
            transcripts.append(line.attr['transcript_id'])
            genes.append(line.attr['gene_id'])
        try:
            line = gtf.next()
        except StopIteration:
            line = ''
   
    return pd.Series(genes, index=transcripts)

def make_gene_info_df(fn):
    """ 
    Make a Pandas dataframe with gene information

    Parameters
    ----------
    fn : str of filename 
        Filename of the Gencode gtf file

    Returns
    -------
    df : pandas.DataFrame
        Pandas dataframe indexed by gene id with the following columns:
        gene_type, gene_status, gene_name.

    """
    import itertools as it

    import HTSeq

    gff_iter = it.islice(HTSeq.GFF_Reader(fn), None)
    convD = dict()
    eof = False
    while not eof:
        try:
            entry = gff_iter.next()
            if entry.type == 'gene':
                convD[entry.attr['gene_id']] = [entry.attr['gene_type'], 
                                                entry.attr['gene_status'],
                                                entry.attr['gene_name']]
        except StopIteration:
            eof = True

    df = pd.DataFrame(convD, index=['gene_type', 'gene_status', 'gene_name']).T
    df.index.name = 'gene_id'
    return df

def make_splice_junction_df(fn, type='gene'):
    """Read the Gencode gtf file and make a pandas dataframe describing the 
    splice junctions

    Parameters
    ----------
    filename : str of filename 
        Filename of the Gencode gtf file

    Returns
    -------
    df : pandas.DataFrame
        Dataframe of splice junctions with the following columns 'gene',
        'chrom', 'start', 'end', 'strand', 'chrom:start', 'chrom:end', 'donor',
        'acceptor', 'intron'

    """
    import itertools as it

    import HTSeq
    import numpy as np

    # GFF_Reader has an option for end_included. However, I think it is
    # backwards.  So if your gtf is end-inclusive, you want the default
    # (end_included=False).  With this, one will NOT be subtracted from the end
    # coordinate.
    gffI = it.islice(HTSeq.GFF_Reader(fn), None)    
    juncL = []
    eof = False
    entry = gffI.next()
    count = 1
    last_count = 1
    while not eof:
        if entry.type == 'transcript':
            exonL = []
            entry = gffI.next()
            count += 1
            gene = entry.attr['gene_id']
            strand = entry.iv.strand
            while not eof and entry.type != 'transcript':
                if entry.type == 'exon':
                    exonL.append(entry)
                try:
                    entry = gffI.next()
                    count += 1
                except StopIteration:
                    eof = True
            # The gencode gtf file has one based, end inclusive coordinates for
            # exons.  HTSeq represents intervals as zero based, end exclusive.
            # We need one-based, end inclusive to compare with STAR output.
            if len(exonL) > 1:
                chrom = exonL[0].iv.chrom
                # On the minus strand, order of exons in gtf file is reversed.
                if strand == '-':
                    exonL.reverse() 
                # We take the exclusive end of the exon intervals and add one to
                # make the one-based start of the intron.
                startL = [ x.iv.end + 1 for x in exonL[:-1] ]   
                # The zero-based inclusive start of the exon is the one-based
                # inclusive end of the intron.
                endL = [ x.iv.start for x in exonL[1:] ]
                for i in range(len(startL)):
                    start = startL[i]
                    end = endL[i]
                    jxn = '{0}:{1}-{2}:{3}'.format(chrom, start, end, strand)
                    chrstart = '{}:{}'.format(chrom, start)
                    chrend = '{}:{}'.format(chrom, end)
                    donor = _gencode_donor(chrom, start, end, strand)
                    acceptor = _gencode_acceptor(chrom, start, end, strand)
                    intron = '{}:{}-{}'.format(chrom, start, end)
                    juncL.append([jxn, gene, chrom, str(start), str(end), 
                                  strand, chrstart, chrend, donor, acceptor,
                                  intron])
        else:
            try:
                entry = gffI.next()
                count += 1
            except StopIteration:
                eof = True
            last_count += 1
    header = ['gene', 'chrom', 'start', 'end', 'strand', 'chrom:start',
              'chrom:end', 'donor', 'acceptor', 'intron']
    juncA = np.array(juncL)
    df = pd.DataFrame(juncA[:,1:], index=juncA[:,0], 
                      columns=header).drop_duplicates() 
    df['start'] = df.start.astype(int)
    df['end'] = df.end.astype(int)
    return df
