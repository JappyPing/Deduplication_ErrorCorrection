Release Notes
================================================================================

Version 0.5.2:  March 8, 2016
-------------------------------------------------------------------------------

Fixed a bug with installation on Windows due to old file paths lingering in
presto.egg-info/SOURCES.txt.

Improvements to commandline usage help messages.

Updated license from CC BY-NC-SA 3.0 to CC BY-NC-SA 4.0.

AssemblePairs:

+ Added the flag ``--fill`` to the reference subcommand to allow insertion of 
  the reference sequence into the non-overlapping region of assembled 
  sequences. Use caution when using this flag, as this may lead to chimeric 
  sequences.
+ Changed default ``--minlen`` to 8 in align subcommand.


Version 0.5.1:  December 4, 2015
-------------------------------------------------------------------------------

ClusterSets:

+ Fixed bug wherein ``--failed`` flag did not work.


Version 0.5.0:  September 7, 2015
-------------------------------------------------------------------------------

Conversion to a proper Python package which uses pip and setuptools for 
installation.

The package now requires Python 3.4. Python 2.7 is not longer supported.

The required dependency versions have been bumped to numpy 1.8, scipy 0.14,
pandas 0.15, and biopython 1.65.

IgCore:

+ Divided IgCore functionality into the separate modules: Annotation, 
  Commandline, Defaults, IO, Multiprocessing and Sequence.


Version 0.4.8:  September 7, 2015
-------------------------------------------------------------------------------

Added support for additional input FASTA (.fna, .fa), FASTQ (.fq) and 
tab-delimited (.tsv) file extensions.

ParseHeaders:

+ Fixed a bug in the rename subcommand wherein renaming to an existing field
  deleted the old annotation, but did not merge the renamed annotation into
  the existing field.
+ Added the copy subcommand which will copy annotations into new field names
  or merge the annotations of existing fields.
+ Added the ``--act`` argument to the copy and rename subcommands allowing 
  collapse following the copy or rename operation.
+ Added a commandline check to ensure that the ``-f``, ``-k`` and ``--act`` 
  arguments contain the same number of fields for both the rename and copy 
  subcommands.


Version 0.4.7:  June 5, 2015
-------------------------------------------------------------------------------

IgCore:

+ Modified scoring functions to permit asymmetrical scores for N and gap 
  characters.
  
AssemblePairs:

+ Added support for SRA style coordinate information where the where the read 
  number has been appended to the spot number.
+ Altered scoring so gap characters are counted as mismatches in the error 
  rate and identity calculations.

BuildConsensus:

+ Altered scoring so gap characters are counted as mismatches in the diversity 
  and error rate calculations.

ConvertHeaders:

+ Added support for SRA style sequence headers where the read number has been 
  appended to the spot number; eg, output from 
  ``fastq-dump -I --split-files file.sra``.

ClusterSets:

+ Added missing OUTPUT console log field.
+ Changed ``--bf`` and ``--cf`` arguments to ``-f`` and ``-k``, respectively.

MaskPrimers:

+ Altering scoring behavior for N characters such that Ns in the input sequence 
  are always counted as a mismatch, while Ns in the primer sequence are counted 
  as a match, with priority given to the input sequence score.
+ Added ``--gap`` argument to the align subcommand which allows users to 
  specify the gap open and gap extension penalties for aligning primers. 
  Note:  gap penalties reduce the match count for purposes of calculating ERROR.

PairSeq:

+ Added support for SRA style coordinate information where the where the read 
  number has been appended to the spot number.


Version 0.4.6:  May 13, 2015
-------------------------------------------------------------------------------

BuildConsensus:

+ Changed ``--maxmiss`` argument to ``--maxgap`` and altered the behavior to 
  only perform deletion of positions based on gap characters (only "-" or "."
  and not "N" characters).
+ Added an error rate (``--maxerror``) calculation based on mismatches from 
  consensus. The ``--maxerror`` argument is mutually exclusive with the 
  ``--maxdiv`` argument and provides similar functionality. However, the 
  calculations are not equivalent, and ``--maxerror`` should be considerably 
  faster than ``--maxdiv``.
+ Added exclusion of positions from the error rate calculation that are deleted
  due to exceeding the ``--maxgap`` threshold .
+ Fixed misalignment of consensus sequence against input sequences when
  positions are deleted due to exceeding the ``--maxgap`` threshold.

ClusterSets:

+ New script to cluster read groups by barcode field (eg, UID barcodes) into
  clustering within the read group.

ConvertHeaders:

+ New script to handle conversion of different sequence description formats 
  to the pRESTO format.
  
FilterSeq:

+ Added count of masked characters to log output of maskqual subcommand.
+ Changed repeats subcommand log field REPEAT to REPEATS.

PairSeq:

+ Changed ``-f`` argument to ``--1f`` argument.
+ Added ``--2f`` argument to copy file 2 annotations to file 1.

ParseHeaders:

+ Moved convert subcommand to the generic subcommand of the new ConvertHeaders 
  script and modified the conversion behavior.


Version 0.4.5:  March 20, 2015
-------------------------------------------------------------------------------

Added details to the usage documentation for each tool which describes both
the output files and annotation fields.

Renamed ``--clean`` argument to ``--failed`` argument with opposite behavior, 
such that the default behavior of all scripts is now clean output.

IgCore:

+ Features added for Change-O compatibility.
+ Features added for PairSeq performance improvements.
+ Added custom help formatter.
+ Modifications to internals of multiprocessing code.
+ Fixed a few typos in error messages.

AssemblePairs:

+ Added reference subcommand which uses V-region germline alignments from
  ublast to assemble paired-ends.
+ Removed mate-pair matching operation to increase performance. Now requires
  both input files to contain matched and uniformly ordered reads. If files
  are not synchronized then PairSeq must be run first. AssemblePairs will
  check that coordinate info matches and error if the files are not
  synchronized. Unpaired reads are no longer output.
+ Added support for cases where one mate pair is the subsequence of the other.
+ Added ``--scanrev`` argument to allow for head sequence to overhang end of 
  tail.
+ Removed truncated (quick) error calculation in align subcommand.
+ Changed default values of the ``--maxerror`` and ``--alpha`` arguments of 
  the align subcommand to better tuned parameters.
+ Changed internal selection of top scoring alignment to use Z-score
  approximation rather than a combination of error rate and binomial
  mid-p value.
+ Internal changes to multiprocessing structure.
+ Changed inserted gap character from - to N in join subcommand for better
  compatibility with the behavior of IMGT/HighV-QUEST.
+ Changed PVAL log field to PVALUE.
+ Changed HEADSEQ and TAILSEQ log fields to SEQ1 and SEQ2.
+ Changed HEADFIELDS and TAILFIELDS log fields to FIELDS1 and FIELDS2.
+ Changed precision of ERROR and PVALUE log fields.
+ Added more verbose logging.

BuildConsensus:

+ Fixed bug where low quality positions where not being masked in single
  sequence barcode groups.
+ Added copy field (``--cf``) and copy action (``--act``) arguments to generate
  consensus annotations for barcode read groups.
+ Changed maximum consensus quality score from 93 to 90.

CollapseSeq:

+ Added ``--keep`` argument to allow retention of sequences with high missing 
  character counts in unique sequence output file.
+ Removed case insensitivity for performance reasons. Now requires all 
  sequences to have matching case.
+ Removed ``first`` and ``last`` from ``--act`` choices to avoid unexpected 
  behavior.

MaskPrimers:

+ Changed behavior of N characters in primer identification. Ns now count as a
  match against any character, rather than a mismatch.
+ Changed behavior of mask mode such that positions masked with Ns are now
  assigned quality scores of 0, rather than retaining their previous scores.
+ Fixed a bug with the align subcommand where deletions within the input
  sequence (gaps in the alignment) were causing an incorrect barcode start
  position.

PairSeq:

+ Performance improvements. The tool should now be considerably faster on very
  large files.
+ Specifying the ``--failed`` argument to request output of sequences which 
  do not have a mate pair will increase run time and memory usage.

ParseHeaders:

+ Add 'cat' action to collapse subcommand which concatenates strings into
  a single annotation.

SplitSeq:

+ Removed ``--clean`` (and ``--failed``) flag from all subcommands.
+ Added progress updates to sample and samplepair subcommands.
+ Performance improvements to samplepair subcommand.


Version 0.4.4:  June 10, 2014
-------------------------------------------------------------------------------

SplitSeq:

+ Removed a linux-specific dependency, allowing SplitSeq to work on Windows.

Version 0.4.3:  April 7, 2014
-------------------------------------------------------------------------------

CollapseSeq:

+ Fixed bug that occurs with Python 2.7.5 on OS X.

SplitSeq:

+ Fixed bug in samplepairs subcommand that occurs with Python 2.7.5 on OS X.


Version 0.4.2:  March 20, 2014
-------------------------------------------------------------------------------

Increased verbosity of exception reporting.

IgCore:

+ Updates to consensus functions to support changes to BuildConsensus.

AssemblePairs:

+ Set default alpha to 0.01.

BuildConsensus:

+ Added support for ``--freq value`` parameter to quality consensus method
  and set default value to 0.6.
+ Fixed a bug in the frequency consensus method where missing values were
  contributing to the total character count at each position.
+ Added the parameter ``--maxmiss value`` which provides a cut-off for 
  removal of positions with too many N or gap characters .

MaskPrimers:

+ Renamed the ``--reverse`` parameter to ``--revpr``.

SplitSeq:

+ Removed convert subcommand.


Version 0.4.1:  January 27, 2014
-------------------------------------------------------------------------------

Changes to the internals of multiple tools to provide support for 
multiprocessing in Windows environments.
  
Changes to the internals of multiple tools to provide clean exit of
child processes upon kill signal or exception in sibling process. 

Fixed unexpected behavior of ``--outname`` and ``--log`` arguments with 
multiple input files.

IgCore:

+ Added reporting of unknown exceptions when reading sequence files
+ Fixed scoring of lowercase sequences.

AlignSets:

+ Fixed a typo in the log output.

BuildConsensus:

+ Fixed a typo in the log output.

EstimateError:

+ Fixed bug where tool would improperly exit if no sets passed threshold
  criteria.
+ Fixed typo in console output.

MaskPrimers:

+ Added ``trim`` mode which will cut the region before primer alignment, but 
  leave primer region unmodified.
+ Fixed a bug with lowercase sequence data.
+ Fixed bug in the console and log output.
+ Added support for primer matching when setting ``--maxerr 1.0``.

ParseHeaders:

+ Added count of sequences without any valid fields (FAIL) to console output.

ParseLog:

+ Added count of records without any valid fields (FAIL) to console output.

SplitSeq:

+ Fixed typo in console output of samplepair subcommand.
+ Added increase of the open file limit to the group subcommand to allow for 
  a large number of groups.


Version 0.4.0:  September 30, 2013
-------------------------------------------------------------------------------

Minor name changes were made to multiple scripts, functions, parameters,
and output files.

AlignSets, AssemblePairs, BuildConsensus, EstimateError, FilterSeq, and 
MaskPrimers are now multithreaded.  The number of simultaneous processes
may be specified using ``--nproc value``. Note this means file ordering
is no longer preserved between the input and output sequence files.

Performance improvements were made to several tools.

The universal ``--verbose`` parameter was replaced with ``--log file_name``
which specifies a log file for verbose output, and disables verbose logging
if not specified.  

The report of input parameters and sequence counts is now separate from the 
log and is always printed to standard output.

Added a progress bar to the standard output of most tools.
  
Added a universal ``--outname file_prefix`` parameter which changes the leading
portion of the output file name.  If not specified, the current file name 
is used (excluding the file extension, as per the previous behavior).

Added a universal ``--clean`` parameter which if specified forces the tool 
not to create an output file of sequences which failed processing.
  
IgCore:

+ Changes to parameters and internals of multiple functions.
+ Added functions to support multithreading for single-end reads, paired-end 
  reads, and barcode sets.
+ Added safe annotation field renaming.
+ Added progress bar, logging and output file name conversion support.
+ Moved reusable AssemblePairs, BuildConsensus, PairSeq, and SplitSeq.
  operations into IgCore.

AssemblePairs:

+ Coordinate information is now specified by a coordinate type, rather than a 
  delimiter, using the ``--coord header_type`` parameter, where the header type
  may be one of ``illumina``, ``solexa``, ``sra``, ``454``, ``presto``.

CollapseSeq:

+ Sequences with a missing character count exceeding the user limit defined
  by ``-n maximum_missing_count`` are now exported to a separate 
  ``collapse-undetermined`` output file, rather than included in the 
  ``collapse-unique`` sequence output.

EstimateError:

+ Now outputs error estimations for positions, quality scores, nucleotide 
  pairs, and annotation sets.  
+ Machine reported quality scores and empirical quality scores have been added
  to all output tables.

FilterSeq:

+ Added ``length`` subcommand to filter sequences by minimum length.

PairSeq:

+ Coordinate information has been redefined as per AssemblePairs.

ParseHeaders:

+ Added new subcommand ``convert`` which attempts to reformat sequence headers 
  into the pRESTO format.
+ The ``rename`` subcommand will now append entries if the new field name already
  exists in the sequence header, rather than replace the entry.


Version 0.3 (prerelease 6):  August 13, 2013
-------------------------------------------------------------------------------

Toolkit is now dependent upon pandas 0.12 for the estimateError tool.

alignSets:

+ Changed MUSCLE execution to faster settings (``-diags``, ``-maxiters 2``).

filterQuality:

+ Added ``repeat`` subcommand to filter sequences with ``-n (value)`` repetitions 
  of a single character and. 
+ Changed ``-n`` parameter of ``ambig`` subcommand from fractional value to a 
  raw count.

estimateError:

+ New tool which estimates error of sequence sets by comparison to a consensus.

maskPrimers:

+ Bug fixes to alignment position calculation of ``align`` subcommand when primer
  alignment begins before start of sequence.
+ Removed ``--ann`` parameter.



Version 0.3 (prerelease 5):  August 7, 2013
-------------------------------------------------------------------------------

License changed to Creative Commons Attribution-NonCommercial-ShareAlike 3.0 
Unported License.

IgPipeline Core:

+ Bug fixes to diversity calculation.
+ Added support for files where all sequences do not share the same annotation 
  fields.
+ Added support for alternate scoring of gap and N-valued nucleotides.

alignSets:

+ Added ``--mode`` parameter with options of ``pad`` and ``cut`` to specify whether 
  to extend or trim read groups to the same start position.
+ Fixed intermittent 'muscle' subcommand stdout pipe deadlock when 
  executing MUSCLE.

assemblePairs:

+ Added ``join`` subcommand to support library preps where paired-end reads 
  do not overlap.
+ Speed improvements to p-value calculations.

buildConsensus:

+ ``--div`` parameter converted to ``--maxdiv value`` to allow filtering of read 
  groups by diversity.
+ Bug fixes to nucleotide frequency consensus method.
+ ``-q`` parameter renamed to ``--qual``.

collapseSequences:

+ Added support for files where all sequences do not share the same annotation 
  fields.

splitSeqFile:

+ ``samplepair`` subcommand added to allow random sampling from paired-end 
  file sets.
+ The behavior of the ``-c`` parameter of the ``sample`` and ``samplepair`` 
  subcommands changed to allow multiple samplings with the same command.


Version 0.3 (prerelease 4):  May 18, 2013
-------------------------------------------------------------------------------

Initial public prerelease
