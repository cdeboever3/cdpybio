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
        in_bam = 'test.bam'
        index = 'test.bam.bai'
        picard_memory = '58G'
        picard_path = 'path/to/picard'
        temp_dir = 'path/to/temp/dir'
        lines = cpb.rnaseq._picard_index(in_bam, index, picard_memory,
                                         picard_path, temp_dir)

class TestBedgraphToBigwig:
    def test_run(self):
        """Test to make sure the function at least runs"""
        bedgraph = 'test.bg'
        bigwig = 'out.bw'
        bedgraph_to_bigwig_path = 'path/to/BedGraphToBigWig'
        bedtools_path = 'path/to/bedtools'
        lines = cpb.rnaseq._bedgraph_to_bigwig(bedgraph, bigwig,
                                               bedgraph_to_bigwig_path,
                                               bedtools_path)

class TestCoverageBedgraph:
    def test_run(self):
        """Test to make sure the function at least runs"""
        bam = 'test.bam'
        bedgraph = 'test.bg'
        bedtools_path = 'path/to/bedtools'
        sample_name = 's1'
        lines = cpb.rnaseq._coverage_bedgraph(bam, bedgraph, bedtools_path,
                                              sample_name)

    def test_run_plus(self):
        """Test to make sure the function at least runs for plus strand"""
        bam = 'test.bam'
        bedgraph = 'test.bg'
        bedtools_path = 'path/to/bedtools'
        sample_name = 's1'
        strand = '+'
        lines = cpb.rnaseq._coverage_bedgraph(bam, bedgraph, bedtools_path,
                                              sample_name, strand)

    def test_run_minus(self):
        """Test to make sure the function at least runs for minus strand"""
        bam = 'test.bam'
        bedgraph = 'test.bg'
        bedtools_path = 'path/to/bedtools'
        sample_name = 's1'
        strand = '-'
        lines = cpb.rnaseq._coverage_bedgraph(bam, bedgraph, bedtools_path,
                                              sample_name, strand)

class TestBigwigFiles:
    def test_run(self):
        """Test to make sure the function at least runs"""
        in_bam = 'test.bam'
        out_bigwig = 'test.bw'
        sample_name = 's1'
        bedgraph_to_bigwig_path = '/path/to/BedGraphToBigWig'
        bedtools_path = 'path/to/bedtools'
        lines = cpb.rnaseq._bigwig_files(in_bam, out_bigwig, sample_name,
                                         bedgraph_to_bigwig_path, bedtools_path)

    def test_run_stranded(self):
        """Test to make sure the function at least runs for stranded output"""
        in_bam = 'test.bam'
        out_bigwig = 'plus.bw'
        sample_name = 's1'
        bedgraph_to_bigwig_path = '/path/to/BedGraphToBigWig'
        bedtools_path = 'path/to/bedtools'
        out_bigwig_minus = 'minus.bw'
        lines = cpb.rnaseq._bigwig_files(in_bam, out_bigwig, sample_name,
                                         bedgraph_to_bigwig_path, bedtools_path,
                                         out_bigwig_minus=out_bigwig_minus)

class TestAlignAndSort:
    def test_run(self):
        """Test to make sure the function at least runs"""
        r1_fastqs = 'r1.fastq.gz'
        r2_fastqs = 'r2.fastq.gz'
        out_dir = 'path/to/outdir'
        sample_name = 's1'
        star_index = 'path/to/index'
        remove_dup=True, 
        strand_specific_cov=False, 
        shell=False
        lines = cpb.rnaseq.align_and_sort(
            r1_fastqs, 
            r2_fastqs, 
            out_dir, 
            sample_name, 
            star_index,
            remove_dup=True, 
            strand_specific_cov=False, 
            shell=False
        )
    
    def test_run_no_remove_dup(self):
        """Test to make sure the function at least runs"""
        r1_fastqs = 'r1.fastq.gz'
        r2_fastqs = 'r2.fastq.gz'
        out_dir = 'path/to/outdir'
        sample_name = 's1'
        star_index = 'path/to/index'
        remove_dup = False,
        lines = cpb.rnaseq.align_and_sort(
            r1_fastqs, 
            r2_fastqs, 
            out_dir, 
            sample_name, 
            star_index,
            remove_dup=remove_dup
        )
    
    def test_run_no_strand_specific_cov(self):
        """Test to make sure the function at least runs"""
        r1_fastqs = 'r1.fastq.gz'
        r2_fastqs = 'r2.fastq.gz'
        out_dir = 'path/to/outdir'
        sample_name = 's1'
        star_index = 'path/to/index'
        strand_specific_cov = True,
        lines = cpb.rnaseq.align_and_sort(
            r1_fastqs, 
            r2_fastqs, 
            out_dir, 
            sample_name, 
            star_index,
            strand_specific_cov=strand_specific_cov
        )
    
    def test_run_shell(self):
        """Test to make sure the function at least runs"""
        r1_fastqs = 'r1.fastq.gz'
        r2_fastqs = 'r2.fastq.gz'
        out_dir = 'path/to/outdir'
        sample_name = 's1'
        star_index = 'path/to/index'
        shell = True 
        lines = cpb.rnaseq.align_and_sort(
            r1_fastqs, 
            r2_fastqs, 
            out_dir, 
            sample_name, 
            star_index,
            shell=True
        )