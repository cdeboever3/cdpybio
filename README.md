cdpybio
==========

Python module with various submodules tailored to specific tools, analyses,
etc. that I perform across projects. Documentation below is only an
introduction, most functions are documented with docstrings.

# Dependencies

* Biopython
* HTSeq
* numpy
* pandas
* pybedtools
* pysam (v0.6 or greater)
* PyVCF
* py.test for testing
* seaborn

All of these are available through conda (though you may need the [bioconda](https://bioconda.github.io/) channel):

```
conda install biopython htseq numpy pandas pybedtools pysam pyvcf seaborn
```

# Submodules

## `bedtools`

Methods for working with bed files. I make heavy use of pybedtools.

## `express`

This submodule has some methods that are useful for dealing with the output
from the RNA-seq expression estimation tool eXpress. 

 * `combine_express_output`: Combine multiple eXpress output files into a
   single pandas dataframe. You can choose which column to combine. You can
also aggregate values by gene ID (eXpress estimates transcript expression) if
you provide a mapping from transcript IDs to gene IDs. See
`gencode.make_transcript_gene_se`.

## `gencode`

Functions for parsing the Gencode gene annotation into various files that are
easier to work with.

 * `make_transcript_gene_se`: Make a file with a simple mapping from transcript
   IDs to gene IDs.
 * `make_gene_info_df`: Make a file indexed by gene ID that has some simple
   information about each gene.
 * `make_splice_junction_df`: Make a file indexed by splice site that has
   information about each splice site such as gene, strand, acceptor, donor,
etc.

## `general`

Some methods that are generally useful.

## `mutect`

Methods for working with MuTect output.

## `pysamext`

Provides extended functionality on top of pysam.

## `star`

 * `def read_sj_out_tab`: Read `sj.out.tab` file from STAR and parse it.
 * `def read_external_annotation`: Read file with junctions from some database.
   This does not have to be the same splice junction database used with STAR.
The file must have some specific columns (see docstring). Compatible with
output from the `gencode.make_splice_junction_df`.
 * `def combine_sj_out`: Combine SJ.out.tab files from STAR by filtering based
   on coverage and comparing to an external annotation to discover novel
junctions.
 * `def make_logs_df`: Make pandas DataFrame from multiple STAR Log.final.out
   files.

## `variants`

Useful for tools for working with DNA variants.
