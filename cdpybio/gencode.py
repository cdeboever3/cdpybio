import pdb
import sys

import pandas as pd
import pybedtools as pbt

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

def load_gffutils_db(f):
    """
    Load database for gffutils. 

    Parameters
    ----------
    f : str
        Path to database.

    Returns
    -------
    db : gffutils.FeatureDB
        gffutils feature database.

    """
    import gffutils
    db = gffutils.FeatureDB(f, keep_order=True)
    return db

def make_gffutils_db(gtf, db):
    """
    Make database for gffutils. 

    Parameters
    ----------
    gtf : str
        Path to Gencode gtf file.

    db : str
        Path to save database to. 

    Returns
    -------
    out_db : gffutils.FeatureDB
        gffutils feature database.

    """
    import gffutils
    out_db = gffutils.create_db(gtf,
                                db,
                                keep_order=True,
                                infer_gene_extent=False)
    return out_db 

def make_promoter_bed(gtf, up=2000, down=200, feature='transcript',
                      use_gene_id=False, merge_by_gene=False, out=None):
    """
    Make a bed file with promoters for transcripts or genes from the Gencode GTF
    file.

    Parameters
    ----------
    gtf : str
        Filename of the Gencode gtf file.

    up : int
        Number of bases to add upstream of the transcription start site.

    down : int
        Number of bases to add downstream of the transcription start site.

    feature : str
        Either transcript or gene. If transcript, promoters, for each transcript
        will be included. If gene, a promoter for each gene entry in the GTF
        will be included.

    use_gene_id : bool
        If True and feature='transcript', the gene ID will be used to name each
        promoter rather than the transcript ID.

    merge_by_gene : bool
        If True and feature='transcript', promoters from all of the transcripts
        of the same gene will be merged and named using the gene ID. This is
        useful for knowing what sequences are promoters for a given gene.

    out : str
        If provided, the bed file will be written to a file with this name.

    Returns
    -------
    bed : pybedtools.BedTool
        A sorted pybedtools BedTool object. 

    """
    import HTSeq
    import itertools as it

    plus_feats = []
    minus_feats = []
    if feature == 'gene' or use_gene_id:
        name_id = 'gene_id'
    elif feature == 'transcript':
        name_id = 'transcript_id'

    gtf = it.islice(HTSeq.GFF_Reader(gtf), None)
    line = gtf.next()
    while line != '':
        if line.type == feature:
            if line.iv.strand == '+':
                plus_feats.append(
                    ('\t'.join([line.iv.chrom, str(line.iv.start - 1),
                                str(line.iv.start - 1),
                                '{}_promoter'.format(line.attr[name_id]),
                                line.iv.strand])))
            elif line.iv.strand == '-':
                minus_feats.append(
                    ('\t'.join([line.iv.chrom, str(line.iv.end - 1),
                                str(line.iv.end - 1),
                                '{}_promoter'.format(line.attr[name_id]),
                                line.iv.strand])))
        try:
            line = gtf.next()
        except StopIteration:
            line = ''

    plus = pbt.BedTool('\n'.join(plus_feats) + '\n', from_string=True)
    minus = pbt.BedTool('\n'.join(minus_feats) + '\n', from_string=True)
    plus = plus.slop(l=up, r=down, g=pbt.chromsizes('hg19'))
    minus = minus.slop(l=down, r=up, g=pbt.chromsizes('hg19'))

    bt  = plus.cat(minus, postmerge=False)
    # We'll sort so bedtools operations can be done faster.
    bt = bt.sort()
    if out:
        bt.saveas(out)
    return bt

def make_feature_bed(gtf, feature, out=None):
    """
    Make a bed file with the start and stop coordinates for all of a particular
    geature in Gencode. Valid features are the features present in the third
    column of the Gencode GTF file.

    Parameters
    ----------
    gtf : str
        Filename of the Gencode gtf file.

    feature : str
        Feature from third column of Gencode GTF file. As of v19, these include
        CDS, exon, gene, Selenocysteine, start_codon, stop_codon, transcript,
        and UTR.

    out : str
        If provided, the bed file will be written to a file with this name.

    Returns
    -------
    bed : pybedtools.BedTool
        A sorted pybedtools BedTool object. 

    """
    bed_lines = []
    with open(gtf) as f:
        line = f.readline().strip()
        while line != '':
            if line[0] != '#':
                line = line.split('\t')
                if line[2] == feature:
                    chrom = line[0]
                    start = str(int(line[3]) - 1)
                    end = line[4]
                    if feature == 'gene':
                        name = line[8].split(';')[0].split(' ')[1].strip('"')
                    else:
                        # TODO: I may want to have some smarter naming for
                        # things that aren't genes or transcripts.
                        name = line[8].split(';')[1].split(' ')[2].strip('"')
                    strand = line[6]
                    bed_lines.append('\t'.join([chrom, start, end, name, '.',
                                                strand]) + '\n')
            line = f.readline().strip()
    bt = pbt.BedTool(''.join(bed_lines), from_string=True)
    # We'll sort so bedtools operations can be done faster.
    bt = bt.sort()
    if out:
        bt.saveas(out)
    return bt

def make_gene_bed(fn, out=None):
    """
    Make a bed file with the start and stop coordinates for each gene. Since
    each gene has a single interval, the introns are by definition included in 
    the interval.

    Parameters
    ----------
    fn : str
        Filename of the Gencode gtf file.

    out : str
        If provided, the bed file will be written to a file with this name.

    Returns
    -------
    bed : pybedtools.BedTool
        A pybedtools BedTool object.

    """
    import pybedtools as pbt
    bed_lines = []
    with open(fn) as f:
        line = f.readline().strip()
        while line != '':
            if line[0] != '#':
                line = line.split('\t')
                if line[2] == 'gene':
                    chrom = line[0]
                    start = str(int(line[3]) - 1)
                    end = line[4]
                    name = line[8].split(';')[0].split(' ')[1].strip('"')
                    strand = line[6]
                    bed_lines.append('\t'.join([chrom, start, end, name, '.',
                                                strand]) + '\n')
            line = f.readline().strip()
    bt = pbt.BedTool(''.join(bed_lines), from_string=True)
    # We'll sort so bedtools operations can be done faster.
    bt = bt.sort()
    if out:
        bt.saveas(out)
    return bt

def make_transcript_gene_se(fn):
    """
    Make a Pandas Series with transcript ID's as the index and values as the
    gene ID containing that transcript.

    Parameters
    ----------
    fn : str
        Filename of the Gencode gtf file.

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
