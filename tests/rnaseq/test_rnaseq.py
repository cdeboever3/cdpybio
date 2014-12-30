import pytest

import cdpybio as cpb

class TestPbsHeader:
    def test_run(self):
        """Test to make sure the function at least runs"""
        out = 'run.out'
        err = 'run.err'
        name = 'run'
        threads = 30
        lines = cpb.rnaseq._pbs_header(out, err, name, threads)

class TestCbarrettPairedDupRemoval:
    def test_run(self):
        """Test to make sure the function at least runs"""
        r1_fastqs = 'r1.fastq.gz'
        r2_fastqs = 'r2.fastq.gz'
        r1_nodup = 'r1_nodup.fastq.gz'
        r2_nodup = 'r2_nodup.fastq.gz'
        temp_dir = 'temp_dir'
        lines = cpb.rnaseq._cbarrett_paired_dup_removal(r1_fastqs, r2_fastqs,
                                                        r1_nodup, r2_nodup,
                                                        temp_dir)
    
    def test_list_run(self):
        """Test to make sure the function at least runs with lists of files as
        input"""
        r1_fastqs = ['r1_1.fastq.gz', 'r1_2.fastq.gz']
        r2_fastqs = ['r2_1.fastq.gz', 'r2_2.fastq.gz']
        r1_nodup = 'r1_nodup.fastq.gz'
        r2_nodup = 'r2_nodup.fastq.gz'
        temp_dir = 'temp_dir'
        lines = cpb.rnaseq._cbarrett_paired_dup_removal(r1_fastqs, r2_fastqs,
                                                        r1_nodup, r2_nodup,
                                                        temp_dir)

class TestStarAlign:
    def test_run(self):
        """Test to make sure the function at least runs"""
        r1_fastqs = 'r1.fastq.gz'
        r2_fastqs = 'r2.fastq.gz'
        sample = 's1'
        rgpl = 'ILLUMINA'
        rgpu = 'flowcell_barcode'
        star_index = 'path/to/star/index'
        star_path = 'path/to/star/executable'
        threads = 30
        lines = cpb.rnaseq._star_align(r1_fastqs, r2_fastqs, sample, rgpl, rgpu,
                                       star_index, star_path, threads)

class TestPicardCoordSort:
    def test_run(self):
        """Test to make sure the function at least runs"""
        in_bam = 'test.bam'
        out_bam = 'test.sorted.bam'
        picard_path = 'picard'
        picard_memory = '58G'
        temp_dir = 'temp_dir'
        lines = cpb.rnaseq._picard_coord_sort(in_bam, out_bam, picard_path,
                                              picard_memory, temp_dir)

class TestPicardIndex:
    def test_run(self):
        """Test to make sure the function at least runs"""
        picard_memory = '58G'
        picard_path = 'path/to/picard'
        temp_dir = 'path/to/temp/dir'
        lines = cpb.rnaseq._picard_index(picard_memory, picard_path, temp_dir)
