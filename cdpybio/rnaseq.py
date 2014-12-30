import os

def _cbarrett_paired_dup_removal(r1_fastqs, r2_fastqs, r1_nodup, r2_nodup,
                                 temp_dir):
    """
    Remove duplicates from paired fastq files using UNIX sort and uniq. Read 
    pairs with exactly the same sequences are removed such that every read pair
    has a different sequence. 

    Parameters
    ----------
    r1_fastqs : str
        R1 fastq file(s). If multiple files, each file should be separated by s
        space and should be ordered the same as the R2 files.

    r2_fastqs : str
        R2 fastq file(s). If multiple files, each file should be separated by s
        space and should be ordered the same as the R1 files.

    r1_nodup : str
        Path to write gzipped R1 fastq file with duplicates removed.

    r2_nodup : str
        Path to write gzipped R2 fastq file with duplicates removed.

    temp_dir : str
        Path to temporary directory where fastq files will be copied to.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    """
    lines = []
    lines.append('paste \\\n')
    lines.append('<(zcat {} | '.format(r1_fastqs) + 
                 'awk \'0==(NR+3)%4{ORS=" "; split($0,a," "); ' + 
                 'print substr(a[1],2)}0==(NR+2)%4{print} (NR!=1 && 0==NR%4)' + 
                 '{ORS="\\n";print}\') \\\n')
    lines.append('<(zcat *_R2_* | '.format(r2_fastqs) + 
                 'awk \'0==(NR+3)%4{ORS=" "; split($0,a," "); ' + 
                 'print substr(a[1],2)}0==(NR+2)%4{print} (NR!=1 && 0==NR%4)' + 
                 '{ORS="\\n";print}\') | \\\n')
    lines.append('awk \'{if ($2 < $5) printf "%s %s %s %s %s %s\\n",'
                 '$1,$3,$4,$6,$2,$5; else printf "%s %s %s %s %s %s\\n",'
                 '$1,$6,$4,$3,$5,$2}\' | \\\n')
    lines.append('sort -k 5,5 -k 6,6 -T {0} -S 30G --parallel=8 | '
                 'uniq -f 4 | \\\n'.format(temp_dir))
    lines.append('awk \'{printf "@%s\\n%s\\n+\\n%s\\n",$1,$5,$2 | '
                 '"gzip -c > ' + r1_nodup + 
                 '"; printf "@%s\\n%s\\n+\\n%s\\n",$3,$6,$4 | "gzip -c > ' + 
                  r2_nodup + '"}\'\n')
    return ''.join(lines)

def _picard_coord_sort(picard_path, picard_memory, temp_dir):
    # Coordinate sort using Picard Tools.
    line = (' \\\n'.join(['java -Xmx{}g -jar '.format(picard_memory),
                          '\t-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
                          '\t-Djava.io.tmpdir={} -jar'.format(temp_dir),
                          '\t{}/SortSam.jar'.format(picard_path),
                          '\tVALIDATION_STRINGENCY=SILENT',
                          '\tI=Aligned.out.bam',
                          '\tO=Aligned.out.coord.sorted.bam',
                          '\tSO=coordinate\n']))
    return line

def _star_align(sample, RGPL, RGPU, star_index, star_path, threads):
    # STAR alignment command.
    line = (' \\\n'.join([star_path,
                          '\t--runThreadN {}'.format(threads - 2),
                          '\t--genomeDir {}'.format(star_index),
                          '\t--genomeLoad NoSharedMemory',
                          '\t--readFilesCommand zcat', 
                          '\t--readFilesIn {} {}'.format('nodup_R1.fastq.gz',
                                                         'nodup_R2.fastq.gz'),
                          '\t--outSAMtype BAM Unsorted',
                          '\t--outSAMattributes All',
                          '\t--outSAMunmapped Within',
                          ('\t--outSAMattrRGline ID:1 ' + 
                           'PL:{} '.format(RGPL) + 
                           'PU:{} '.format(RGPU) + 
                           'LB:{0} SM:{0}'.format(sample)),
                          '\t--outFilterMultimapNmax 20',
                          '\t--outFilterMismatchNmax 999',
                          '\t--outFilterMismatchNoverLmax 0.04',
                          ('\t--outFilterIntronMotifs '
                           'RemoveNoncanonicalUnannotated'),
                          '\t--outSJfilterOverhangMin 6 6 6 6',
                          '\t--seedSearchStartLmax 20',
                          '\t--alignSJDBoverhangMin 1',
                          '\t--quantMode TranscriptomeSAM']) + '\n\n')
    return line

def _picard_index(picard_memory, picard_path, temp_dir):
    # Index bam file using Picard Tools.
    line = (' \\\n'.join(['java -Xmx{}g -jar'.format(picard_memory),
                          '\t-XX:-UseGCOverheadLimit -XX:-UseParallelGC',
                          '\t-Djava.io.tmpdir={} -jar'.format(temp_dir),
                          '\t{}/BuildBamIndex.jar'.format(picard_path),
                          '\tI=Aligned.out.coord.sorted.bam',
                          '\tO=Aligned.out.coord.sorted.bam.bai\n\n']))
    return line

def _pbs_header(out, err, name, threads, queue='high'):
    """
    Write header for PBS script

    Parameters
    ----------
    out : str
        Path to file for recording standard out.

    err : str
        Path to file for recording standard error.

    Returns
    -------
    lines : str
        Lines to be printed to shell/PBS script.

    """
    lines = ('\n'.join(['#PBS -q {}'.format(queue),
                        '#PBS -N {}'.format(name),
                        '#PBS -l nodes=1:ppn={}'.format(threads),
                        '#PBS -o {}'.format(out),
                        '#PBS -e {}\n\n'.format(err)]))
    return lines

def _bigwig_files():
    # Make bigwig coverage files.
    if strand_specific_cov:
        f.write(_coverage(TODO_bam, bedtools, sample_name, TODO_out_fn,
                          strand='+'))
        f.write(_coverage(TODO_bam, bedtools, sample_name, TODO_out_fn,
                          strand='-'))
    else:
        f.write(_coverage(TODO_bam, bedtools, sample_name, TODO_out_fn))
    f.write('wait\n\n')
    # TODO: move bedgraph to bigwig into functions
    f.write( ' '.join(['{} plus.bg'.format(ppy.bedGraphToBigWig),
                       '{}/genomes/human.hg19.genome'.format(
        os.path.split(os.path.split(ppy.bedtools)[0])[0]),
                       'plus.bw &\n']))
    f.write( ' '.join(['{} minus.bg'.format(ppy.bedGraphToBigWig),
                       '{}/genomes/human.hg19.genome'.format(
        os.path.split(os.path.split(ppy.bedtools)[0])[0]),
                       'minus.bw\n\n']))

def _process_fastqs(fastqs, temp_dir):
    """
    Create list of temporary fastq paths.

    Parameters
    ----------
    fastqs : list or str
        Either a list of paths to gzipped fastq files or path to a single
        gzipped fastq file.

    temp_dir : str
        Path to temporary directory where fastq files will be copied to.

    Returns
    -------
    fastqs : str
        Paths to original fastq files (concatenated with a space if multiple,
        e.g. 'path/fq1.fastq path/fq2.fastq').

    temp_fastqs : str
        Paths to temporary fastq files (concatenated with a space if multiple,
        e.g. 'tempdir/fq1.fastq tempdir/fq2.fastq').

    """
    if type(fastqs) == list:
        fns = [os.path.split(x)[1] for x in fastqs]
        temp_fastqs = [os.path.join(temp_dir, x) for x in fns]
        fastqs = ' '.join(fastqs)
    elif type(fastqs) == str:
        temp_fastqs = os.path.join(temp_dir, os.path.split(fastqs)[1])
    return fastqs, temp_fastqs

def align_and_sort(
    r1_fastqs, 
    r2_fastqs, 
    out_dir, 
    sample_name, 
    star_index,
    RGPL='',
    RGPU='',
    picard_path='',
    bedtools_path='',
    temp_dir='/scratch', 
    threads=32, 
    picard_memory=58, 
    remove_dup=True, 
    strand_specific_cov=False, 
    shell=False
):
    """
    Make a PBS or shell script for aligning RNA-seq reads with STAR. The
    defaults are set for use on the Frazer lab's PBS scheduler on FLC.

    Parameters
    ----------
    r1_fastqs : list or str
        Either a list of paths to gzipped fastq files with R1 reads or path to a
        single gzipped fastq file with R1 reads.

    r2_fastqs : list or str
        Either a list of paths to gzipped fastq files with R2 reads or path to a
        single gzipped fastq file with R2 reads.

    out_dir : str
        Directory to store PBS/shell file and aligment results.

    sample_name : str
        Sample name used for naming files etc.

    star_index : str
        Path to STAR index.

    RGPL : str
        Read Group platform (e.g. illumina, solid). RGPL and RGPU are
        required to add read group information.

    RGPU : str
        Read Group platform unit (eg. run barcode). RGPL and RGPU are required
        to add read group information.

    star_path : str
        Path to STAR aligner. If not provided, assumed to be in your path.

    picard_path : str
        Path to Picard tools. If not provided, assumed to be in your path.

    bedtools_path : str
        Path to bedtools. If not provided, assumed to be in your path.

    temp_dir : str
        Directory to store files as STAR runs.

    threads : int
        Number of threads to reserve using PBS scheduler. This number of threads
        minus 2 will be used by STAR, so this must be at least 3.

    picard_memory : int
        Amount of memory (in gb) to give Picard Tools.

    remove_dup : boolean
        Whether to remove duplicate reads prior to alignment.

    strand_specific_cov : boolean
        If true, make strand specific bigwig files. 

    shell : boolean
        If true, make a shell script rather than a PBS script.

    Returns
    -------
    fn : str
        Path to PBS/shell script.

    """
    assert threads >= 3

    if shell:
        pbs = False
    else: 
        pbs = True

    # I'm going to define some file names used later.
    r1_fastqs, temp_r1_fastqs = _process_fastqs(r1_fastqs, temp_dir)
    r2_fastqs, temp_r2_fastqs = _process_fastqs(r2_fastqs, temp_dir)
    r1_nodup = os.path.join(temp_dir, 'nodup_R1.fastq.gz')
    r2_nodup = os.path.join(temp_dir, 'nodup_R2.fastq.gz')
    
    # TODO: update these lists
    # Files to copy to output directory.
    files_to_copy = [out_bamN,'Log.out','Log.final.out','Log.progress.out','SJ.out.tab']
    # Temporary files that can be deleted at the end of the job. We may not want
    # to delete the temp directory if the temp and output directory are the
    # same.
    files_to_remove = [temp_r1_fastqs, temp_r2_fastqs, r1_nodup, r2_nodup]

    try:
        os.makedirs(out_dir)
    except OSError:
        pass

    fn = os.path.join(out_dir, '{}_alignment.pbs'.format(sample_name))
    f = open(fn, 'w')
    f.write('#!/bin/bash\n\n')
    if pbs:
        out = os.path.join(out_dir, '{}_alignment.out'.format(sample_name))
        err = os.path.join(out_dir, '{}_alignment.err'.format(sample_name))
        job_name = '{}_align'.format(sample)
        _pbs_header(out, err, job_name, threads)

    f.write('mkdir -p {}\n'.format(temp_dir))
    f.write('cd {}\n'.format(temp_dir))
    f.write('rsync -avz {} {} .\n\n'.format(r1_fastqs, r2_fastqs))

    if remove_dup:
        lines = _cbarrett_paired_dup_removal(temp_r1_fastqs, temp_r2_fastqs,
                                             r1_nodup, r2_nodup, temp_dir)
        f.write(lines)

    lines = _star_align(sample, RGPL, RGPU, star_index, star_path, threads)
    f.write(lines)
    # _picard_coord_sort(picard_path, picard_memory, temp_dir)
    # _picard_index(picard_memory, picard_path, temp_dir)

    # f.write('wait\n\n')
    # f.write('cp {} {}\n\n'.format(' '.join(['Aligned.out.coord.sorted.bam',
    #                                       'Aligned.out.coord.sorted.bam.bai',
    #                                       'Log.progress.out',
    #                                       'SJ.out.tab',
    #                                       'Aligned.toTranscriptome.out.bam',
    #                                       'Log.final.out',
    #                                       'Log.out',
    #                                       'plus.bw',
    #                                       'minus.bw']),
    #                             alignment))
    # f.write('cp {} {}\n\n'.format('SJ.out.tab',
    #                               os.path.join(ppy.root, 'data', 
    #                                            'SJ.out.tab', 
    #                                            '{}.SJ.out.tab'.format(sample))))
    # f.write('cp {} {}\n\n'.format('Log.final.out',
    #                               os.path.join(ppy.root, 'data', 
    #                                            'Log.final.out', 
    #                                            '{}.Log.final.out'.format(sample))))
    # f.write('rm -r {}\n\n'.format(temp_dir))
    f.close()
    return os.path.join(alignment, 'alignment.pbs')

def _coverage_bg(bam, bedtools, sample_name, out_fn, strand='.'):
    """
    Make lines that create a coverage bedgraph file.

    Parameters
    ----------
    bam : str
        Bam file to calculate coverage for.

    bedtools : str
        Path to bedtools.

    sample_name : str
        Sample name for naming files etc.

    out_fn : str
        Output coverage file.
    
    strand : str
        If '+' or '-', calculate strand-specific coverage. Otherwise, calculate
        coverage using all reads.

    Returns
    -------
    lines : str
        Lines to be written to PBS/shell script.

    """
    if strand == '+' or strand == '-':
        if strand == '+':
            name = '{}_plus'.format(sample_name)
        else:
            name = '{}_minus'.format(sample_name)
        lines = ' \\\n'.join(['{} genomecov -ibam'.format(bedtools),
                              '\t{}'.format(bam),
                              '\t-g hg19.genome -split -bg ',
                              '\t-strand {} -trackline'.format(strand),
                              '\t-trackopts \'name="{}"\''.format(name),
                              '\t> {}.bg &\n\n'.format(name)])
    else:
        name = sample_name
        lines = ' \\\n'.join(['{} genomecov -ibam'.format(bedtools),
                              '\t{}'.format(bam),
                              '\t-g hg19.genome -split -bg ',
                              '\t-trackline'.format(strand),
                              '\t-trackopts \'name="{}"\''.format(name),
                              '\t> {}.bg &\n\n'.format(name)])
    return lines

def ucsc_files():
    f.write('ln -s {} {}\n'.format(
        os.path.join(alignment, 'Aligned.out.coord.sorted.bam'),
        os.path.join(local_ucsc_bam_path, '{}.bam'.format(sample))))
    f.write('ln -s {} {}\n'.format(
        os.path.join(alignment, 'Aligned.out.coord.sorted.bam.bai'),
        os.path.join(local_ucsc_bam_path, '{}.bam.bai'.format(sample))))
    f.write('ln -s {} {}\n'.format(
        os.path.join(alignment, 'plus.bw'),
        os.path.join(local_ucsc_bigwig_path, '{}_plus.bw'.format(sample))))
    f.write('ln -s {} {}\n'.format(
        os.path.join(alignment, 'minus.bw'),
        os.path.join(local_ucsc_bigwig_path, '{}_minus.bw'.format(sample))))
    f.close()
    # Write file with UCSC track lines
    f = open(os.path.join(alignment, 'tracklines.txt'), 'w')
    f.write(' '.join(['track', 'type=bigWig', 'name="{}_plus"'.format(sample),
                      'description="Plus strand coverage for {}"'.format(sample),
                      'bigDataUrl={}/{}_plus.bw\n'.format(web_ucsc_bigwig_path,
                                                          sample)]))
    f.write(' '.join(['track', 'type=bigWig', 'name="{}_minus"'.format(sample),
                      'description="Minus strand coverage for {}"'.format(sample),
                      'bigDataUrl={}/{}_minus.bw\n'.format(web_ucsc_bigwig_path,
                                                           sample)]))
    f.write(' '.join(['track', 'type=bam', 'name="{}_bam"'.format(sample),
                      'description="Bam file for {}"'.format(sample),
                      'bigDataUrl={}/{}.bam\n'.format(web_ucsc_bam_path,
                                                      sample)]))
