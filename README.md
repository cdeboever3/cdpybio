cdpybio
==========

Python module with various submodules tailored to specific tools, analyses,
etc. that I perform across projects. Documentation below is only an
introduction, most functions are documented with docstrings.

# Dependencies
* HTSeq
* numpy
* pandas
* pybedtools
* py.test for testing

# Submodules

## `cghub`

This submodule contains some useful methods and classes for dealing with data
through CGHub. I'll go through the classes and methods.

### `GTFuseBam`

A `GTFuseBam` object is a single bam file from CGHub mounted with GTFuse. You
can `mount` and `unmount` the bam file as you'd like.

### `ReadsFromIntervalsBam`

This object takes a GTFuseBam object and a set of intervals and obtains the
reads from those intervals in the CGHub bam file and writes them to a local bam
file.

* `reads_from_intervals`: Obtain the reads from the given intervals.

### `ReadsFromIntervalsEngine`

This is an engine that runs in the background and obtains reads from intervals
for a given set of samples. The main process that runs the engine shares the
thread with your python session but I use the multiprocessing module to farm out
different bam files to different threads so you can obtain reads from multiple
bam files simultaneously. This class can be extended.

### `FLCVariantCallingEngine(ReadsFromIntervalsEngine)`

This engine extends the `ReadsFromIntervalsEngine` and performs variant calling
after obtaining the reads. Currently implemented to work in the Frazer lab
computing environment although it would be easy to change for a different
computing environment.

### `TumorNormalVariantCall`

Class that wraps the results of a variant calling job for a tumor/normal pair.

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
