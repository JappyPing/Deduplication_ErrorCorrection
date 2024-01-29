#!/usr/bin/env python3
"""
Assembles paired-end reads into a single sequence
"""
# Info
__author__ = 'Jason Anthony Vander Heiden, Gur Yaari, Christopher Bolen'
from presto import __version__, __date__

# Imports
import os
import sys
import numpy as np
import pandas as pd
import scipy.stats as stats
from argparse import ArgumentParser
from collections import OrderedDict
from io import StringIO
from textwrap import dedent
from time import time
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Presto imports
from presto.Defaults import default_delimiter, default_coord_choices, \
                            default_coord_type, default_missing_chars, \
                            default_out_args, default_usearch_exec
from presto.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from presto.Annotation import parseAnnotation, flattenAnnotation, mergeAnnotation, \
                              getCoordKey
from presto.Applications import runUBlastAlignment
from presto.IO import getFileType, readSeqFile, countSeqFile, getOutputHandle, \
                      printLog, printProgress
from presto.Sequence import getDNAScoreDict, reverseComplement, scoreSeqPair
from presto.Multiprocessing import SeqData, SeqResult, manageProcesses, processSeqQueue

# Defaults
default_alpha = 1e-5
default_max_error = 0.3
default_min_ident = 0.5
default_min_len = 8
default_max_len = 1000
default_gap = 0
default_evalue = 1e-5
default_max_hits = 100

class AssemblyRecord:
    """
    A class defining a paired-end assembly result
    """
    # Instantiation
    def __init__(self, seq=None):
        self.seq = seq
        self.ref_seq = None
        self.head_pos = None
        self.tail_pos = None
        self.ref_pos = None
        self.gap = None
        self.zscore = float('-inf')
        self.pvalue = None
        self.evalue = None
        self.error = None
        self.ident = None
        self.valid = False

    # Set boolean evaluation to valid value
    def __bool__(self):
        return self.valid

    # Set length evaluation to length of SeqRecord
    def __len__(self):
        if self.seq is None:
            return 0
        else:
            return len(self.seq)

    # Set overlap length to head_pos difference
    @property
    def overlap(self):
        if self.head_pos is None:
            return None
        else:
            return self.head_pos[1] - self.head_pos[0]


class AssemblyStats:
    """
    Class containing p-value and z-score matrices for scoring assemblies
    """
    # Instantiation
    def __init__(self, n):
        self.p = AssemblyStats._getPMatrix(n)
        self.z = AssemblyStats._getZMatrix(n)
        #print self.z

    @staticmethod
    def _getPMatrix(n):
        """
        Generates a matrix of mid-p correct p-values from a binomial distribution

        Arguments:
        n = maximum trials

        Returns:
        a numpy.array of successes by trials p-values
        """
        p_matrix = np.empty([n, n], dtype=float)
        p_matrix.fill(np.nan)
        k = np.arange(n, dtype=float)
        for i, x in enumerate(k):
            p_matrix[x, i:] = 1 - stats.binom.cdf(x - 1, k[i:], 0.25) - stats.binom.pmf(x, k[i:], 0.25) / 2.0
        return p_matrix

    @staticmethod
    def _getZMatrix(n):
        """
        Generates a matrix of z-score approximations for a binomial distribution

        Arguments:
        n = maximum trials

        Returns:
        a numpy.array of successes by trials z-scores
        """
        z_matrix = np.empty([n, n], dtype=float)
        z_matrix.fill(np.nan)
        k = np.arange(0, n, dtype=float)
        for i, x in enumerate(k):
            j = i + 1 if i == 0 else i
            z_matrix[x, j:] = (x - k[j:]/4.0)/np.sqrt(3.0/16.0*k[j:])
        return z_matrix


def referenceAssembly(head_seq, tail_seq, ref_dict, ref_file, min_ident=default_min_ident,
                      evalue=default_evalue, max_hits=default_max_hits, fill=False,
                      usearch_exec=default_usearch_exec,
                      score_dict=getDNAScoreDict(mask_score=(1, 1), gap_score=(0, 0))):
    """
    Stitches two sequences together by aligning against a reference database

    Arguments:
    head_seq = the head SeqRecord
    head_seq = the tail SeqRecord
    ref_dict = a dictionary of reference SeqRecord objects
    ref_file = the path to the reference database file
    min_ident = the minimum identity for a valid assembly
    evalue = the E-value cut-off for ublast
    max_hits = the maxhits output limit for ublast
    fill = if False non-overlapping regions will be assigned Ns;
           if True non-overlapping regions will be filled with the reference sequence.
    usearch_exec = the path to the usearch executable
    score_dict = optional dictionary of character scores in the
                 form {(char1, char2): score}

    Returns:
    an AssemblyRecord object
    """
    # Define general parameters
    head_len = len(head_seq)
    tail_len = len(tail_seq)

    # Determine if quality scores are present
    has_quality = hasattr(head_seq, 'letter_annotations') and \
                  hasattr(tail_seq, 'letter_annotations') and \
                  'phred_quality' in head_seq.letter_annotations and \
                  'phred_quality' in tail_seq.letter_annotations

    # Align against reference
    head_df = runUBlastAlignment(head_seq, ref_file, evalue=evalue, max_hits=max_hits,
                                 usearch_exec=usearch_exec)
    tail_df = runUBlastAlignment(tail_seq, ref_file, evalue=evalue, max_hits=max_hits,
                                 usearch_exec=usearch_exec)

    # Subset results to matching reference assignments
    align_df = pd.merge(head_df, tail_df, on='target', how='inner', suffixes=('_head', '_tail'))

    # If no matching targets return failed results
    if len(align_df) < 1:
        return AssemblyRecord()

    # Select top alignment
    align_top = align_df.ix[0, :]
    ref_id = align_top['target']
    ref_seq = ref_dict[ref_id]

    # Get offset of target and reference positions
    head_shift = align_top['target_start_head'] - align_top['query_start_head']
    tail_shift = align_top['target_start_tail'] - align_top['query_start_tail']

    # Get positions of outer reference match in head (a, b) and tail (x, y) sequences
    outer_start = align_top[['target_start_head', 'target_start_tail']].min()
    outer_end = align_top[['target_end_head', 'target_end_tail']].max()
    a_outer = outer_start - head_shift
    b_outer = outer_end - head_shift
    x_outer = outer_start - tail_shift
    y_outer = outer_end - tail_shift

    # Get positions of inner reference match in head (a,b) and tail (x,y) sequences
    inner_start = align_top[['target_start_head', 'target_start_tail']].max()
    inner_end = align_top[['target_end_head', 'target_end_tail']].min()
    a_inner = inner_start - head_shift
    b_inner = inner_end - head_shift
    x_inner = inner_start - tail_shift
    y_inner = inner_end - tail_shift

    # Determine head (a, b) and tail (x, y) overlap positions
    a = max(0, a_inner - x_inner)
    b = min(b_inner + (tail_len - y_inner), head_len)
    x = max(0, x_inner - a_inner)
    y = min(y_inner + (head_len - b_inner), tail_len)

    # Join sequences if head and tail do not overlap, otherwise assemble
    if a > b and x > y:
        stitch = joinSeqPair(head_seq, tail_seq, gap=(a - b), insert_seq=None)
    else:
        stitch = AssemblyRecord()
        stitch.gap = 0

        # Define overlap sequence
        if has_quality:
            # Build quality consensus
            overlap_seq = overlapConsensus(head_seq[a:b], tail_seq[x:y])
        else:
            # Assign head sequence to conflicts when no quality information is available
            overlap_seq = head_seq[a:b]

        # Assemble sequence
        if a > 0 and y < tail_len:
            # Tail overlaps end of head
            stitch.seq = head_seq[:a] + overlap_seq + tail_seq[y:]
        elif b < head_len and x > 0:
            # Head overlaps end of tail
            stitch.seq = tail_seq[:x] + overlap_seq + head_seq[b:]
        elif a == 0 and b == head_len:
            # Head is a subsequence of tail
            stitch.seq = tail_seq[:x] + overlap_seq + tail_seq[y:]
        elif x == 0 and y == tail_len:
            # Tail is a subsequence of head
            stitch.seq = head_seq[:a] + overlap_seq + head_seq[b:]
        else:
            sys.stderr.write('ERROR:  Invalid overlap condition for %s\n' % head_seq.id)

        # Define stitch ID
        stitch.seq.id = head_seq.id if head_seq.id == tail_seq.id \
                                    else '+'.join([head_seq.id, tail_seq.id])
        stitch.seq.name = stitch.seq.id
        stitch.seq.description = ''

    # Assign position info
    stitch.head_pos = (a, b)
    stitch.tail_pos = (x, y)

    # Assign reference info
    stitch.ref_seq = ref_seq[outer_start:outer_end]
    stitch.ref_pos = (max(a_outer, x_outer), max(b_outer, y_outer))
    stitch.evalue = tuple(align_top[['evalue_head', 'evalue_tail']])

    # Calculate assembly error
    score, weight, error = scoreSeqPair(stitch.seq.seq[stitch.ref_pos[0]:stitch.ref_pos[1]],
                                        ref_seq.seq[outer_start:outer_end],
                                        score_dict=score_dict)
    stitch.ident = 1 - error
    stitch.valid = bool(stitch.ident >= min_ident)

    # Fill gap with reference if required
    if a > b and x > y and fill:
        insert_seq = ref_seq.seq[(b + head_shift):(a + head_shift)]
        insert_rec = joinSeqPair(head_seq, tail_seq, gap=(a - b), insert_seq=insert_seq)
        stitch.seq = insert_rec.seq

    return stitch


def overlapConsensus(head_seq, tail_seq, ignore_chars=default_missing_chars):
    """
    Creates a consensus overlap sequences from two segments

    Arguments: 
    head_seq = the overlap head SeqRecord
    tail_seq = the overlap tail SeqRecord
    ignore_chars = list of characters which do not contribute to consensus
    
    Returns:
    A SeqRecord object with consensus characters and quality scores
    """
    # Initialize empty overlap character and quality score list
    seq_cons, score_cons = [], []
    # Define character and quality tuple iterators
    chars = list(zip(head_seq, tail_seq))
    quals = list(zip(head_seq.letter_annotations['phred_quality'], 
                 tail_seq.letter_annotations['phred_quality']))

    # Iterate over character and quality tuples and build consensus
    for c, q in zip(chars, quals):
        # Equivalent character case
        if c[0] == c[1]:
            c_cons = c[0]
            q_cons = max(q)
        # All ambiguous characters case
        elif all([x in ignore_chars for x in c]):
            c_cons = 'N'
            q_cons = max(q)
        # Some ambiguous characters case
        elif any([x in ignore_chars for x in c]):
            c_cons = [x for x in c if x not in ignore_chars][0]
            q_cons = q[c.index(c_cons)]
        # Conflicting character case        
        else:
            q_max = max(q)
            c_cons = c[q.index(q_max)]
            try:
                q_cons = int(q_max**2 / sum(q))
            except ZeroDivisionError:
                q_cons = 0
        # Append sequence and quality lists with consensus values
        seq_cons.append(c_cons)
        score_cons.append(q_cons)

    # Define overlap SeqRecord
    record = SeqRecord(Seq(''.join(seq_cons), IUPAC.ambiguous_dna), 
                       id='', 
                       name='', 
                       description='', 
                       letter_annotations={'phred_quality':score_cons})
    
    return record


def joinSeqPair(head_seq, tail_seq, gap=default_gap, insert_seq=None):
    """
    Concatenates two sequences 

    Arguments: 
    head_seq = the head SeqRecord
    tail_seq = the tail SeqRecord
    gap = number of gap characters to insert between head and tail
          ignored if insert_seq is not None.
    insert_seq = a string or Bio.Seq.Seq object, to insert between the head and tail;
                 if None insert with N characters

    Returns: 
    an AssemblyRecord object
    """
    # Define joined ID
    join_id = head_seq.id if head_seq.id == tail_seq.id \
              else '+'.join([head_seq.id, tail_seq.id])

    # Join sequences
    if insert_seq is None:
        join_seq = str(head_seq.seq) + 'N' * gap + str(tail_seq.seq)
    else:
        gap = len(insert_seq)
        join_seq = str(head_seq.seq) + str(insert_seq) + str(tail_seq.seq)
    
    # Define return record
    record = SeqRecord(Seq(join_seq, IUPAC.ambiguous_dna), 
                       id=join_id, 
                       name=join_id, 
                       description='')
    
    # Join quality score if present
    has_quality = hasattr(head_seq, 'letter_annotations') and \
                  hasattr(tail_seq, 'letter_annotations') and \
                  'phred_quality' in head_seq.letter_annotations and \
                  'phred_quality' in tail_seq.letter_annotations
    if has_quality:
        join_quality = head_seq.letter_annotations['phred_quality'] + \
                       [0] * gap + \
                       tail_seq.letter_annotations['phred_quality']
        record.letter_annotations = {'phred_quality':join_quality}

    stitch = AssemblyRecord(record)
    stitch.valid = True
    stitch.gap = gap

    return stitch


def alignAssembly(head_seq, tail_seq, alpha=default_alpha, max_error=default_max_error,
                  min_len=default_min_len, max_len=default_max_len, scan_reverse=False,
                  assembly_stats=None, score_dict=getDNAScoreDict(mask_score=(1, 1), gap_score=(0, 0))):
    """
    Stitches two sequences together by aligning the ends

    Arguments:
    head_seq = the head SeqRecord
    head_seq = the tail SeqRecord
    alpha = the minimum p-value for a valid assembly
    max_error = the maximum error rate for a valid assembly
    min_len = minimum length of overlap to test
    max_len = maximum length of overlap to test
    scan_reverse = if True allow the head sequence to overhang the end of the tail sequence
                   if False end alignment scan at end of tail sequence or start of head sequence
    assembly_stats = optional successes by trials numpy.array of p-values
    score_dict = optional dictionary of character scores in the 
                 form {(char1, char2): score}
                     
    Returns: 
    an AssemblyRecord object
    """
    # Set alignment parameters
    if assembly_stats is None:  assembly_stats = AssemblyStats(max_len + 1)

    # Define general parameters
    head_str = str(head_seq.seq)
    tail_str = str(tail_seq.seq)
    head_len = len(head_str)
    tail_len = len(tail_str)

    # Determine if quality scores are present
    has_quality = hasattr(head_seq, 'letter_annotations') and \
                  hasattr(tail_seq, 'letter_annotations') and \
                  'phred_quality' in head_seq.letter_annotations and \
                  'phred_quality' in tail_seq.letter_annotations

    # Determine if sub-sequences are allowed and define scan range
    if scan_reverse and max_len >= min(head_len, tail_len):
        scan_len = head_len + tail_len - min_len
    else:
        scan_len = min(max(head_len, tail_len), max_len)

    # Iterate and score overlap segments
    stitch = AssemblyRecord()
    for i in range(min_len, scan_len + 1):
        a = max(0, head_len - i)
        b = head_len - max(0, i - tail_len)
        x = max(0, i - head_len)
        y = min(tail_len, i)
        score, weight, error = scoreSeqPair(head_str[a:b], tail_str[x:y], score_dict=score_dict)
        z = assembly_stats.z[score, weight]
        # Save stitch as optimal if z-score improves
        if z > stitch.zscore:
           stitch.head_pos = (a, b)
           stitch.tail_pos = (x, y)
           stitch.zscore = z
           stitch.pvalue = assembly_stats.p[score, weight]
           stitch.error = error

    # Build stitched sequences and assign best_dict values
    if stitch.head_pos is not None:
        # Correct quality scores and resolve conflicts
        a, b = stitch.head_pos
        x, y = stitch.tail_pos
        if has_quality:
            # Build quality consensus
            overlap_seq = overlapConsensus(head_seq[a:b], tail_seq[x:y])
        else:
            # Assign head sequence to conflicts when no quality information is available
            overlap_seq = head_seq[a:b]

        if a > 0 and y < tail_len:
            # Tail overlaps end of head
            stitch.seq = head_seq[:a] + overlap_seq + tail_seq[y:]
        elif b < head_len and x > 0:
            # Head overlaps end of tail
            stitch.seq = tail_seq[:x] + overlap_seq + head_seq[b:]
        elif a == 0 and b == head_len:
            # Head is a subsequence of tail
            stitch.seq = tail_seq[:x] + overlap_seq + tail_seq[y:]
        elif x == 0 and y == tail_len:
            # Tail is a subsequence of head
            stitch.seq = head_seq[:a] + overlap_seq + head_seq[b:]
        else:
            sys.stderr.write('ERROR:  Invalid overlap condition for %s\n' % head_seq.id)


        # Define best stitch ID
        stitch.seq.id = head_seq.id if head_seq.id == tail_seq.id \
                              else '+'.join([head_seq.id, tail_seq.id])
        stitch.seq.name = stitch.seq.id
        stitch.seq.description = ''

    stitch.valid = bool(stitch.pvalue <= alpha and stitch.error <= max_error)

    return stitch


def feedPairQueue(alive, data_queue, seq_file_1, seq_file_2,
                  coord_type=default_coord_type, delimiter=default_delimiter):
    """
    Feeds the data queue with sequence pairs for processQueue processes

    Arguments:
    alive = a multiprocessing.Value boolean controlling whether processing 
            continues; when False function returns
    data_queue = an multiprocessing.Queue to hold data for processing
    seq_file_1 = the name of sequence file 1
    seq_file_2 = the name of sequence file 2
    coord_type = the sequence header format
    delimiter =  a tuple of delimiters for (fields, values, value lists)

    Returns: 
    None
    """
    # Function to get coordinate info
    def _key_func(x):
        return getCoordKey(x, coord_type=coord_type, delimiter=delimiter)

    # Generator function to read and check files
    def _read_pairs(seq_file_1, seq_file_2):
        iter_1 = readSeqFile(seq_file_1, index=False)
        iter_2 = readSeqFile(seq_file_2, index=False)
        for seq_1, seq_2 in zip(iter_1, iter_2):
            key_1 = getCoordKey(seq_1.description, coord_type=coord_type,
                                delimiter=delimiter)
            key_2 = getCoordKey(seq_2.description, coord_type=coord_type,
                                delimiter=delimiter)
            if key_1 == key_2:
                yield (key_1, [seq_1, seq_2])
            else:
                raise Exception('Coordinates for sequences %s and %s do not match' \
                                 % (key_1, key_2))

    try:
        # Open and parse input files
        data_iter = _read_pairs(seq_file_1, seq_file_2)

        # Iterate over data_iter and feed data queue 
        while alive.value:
            # Get data from queue
            if data_queue.full():  continue
            else:  data = next(data_iter, None)
            # Exit upon reaching end of iterator
            if data is None:  break

            # Feed queue
            data_queue.put(SeqData(*data))
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        raise

    return None


def processAssembly(data, assemble_func, assemble_args={}, rc=None,
                   fields_1=None, fields_2=None, delimiter=default_delimiter):
    """
    Performs assembly of a sequence pair

    Arguments:
    data = a SeqData object with a list of exactly two SeqRecords
    assemble_func = the function to use to assemble paired ends
    assemble_args = a dictionary of arguments to pass to the assembly function
    rc = Defines which sequences ('head','tail','both') to reverse complement
         before assembly; if None do not reverse complement sequences
    fields_1 = list of annotations in head SeqRecord to copy to assembled record;
               if None do not copy an annotation
    fields_2 = list of annotations in tail SeqRecord to copy to assembled record;
               if None do not copy an annotation
    delimiter = a tuple of delimiters for (fields, values, value lists)

    Returns:
    a SeqResult object
    """
    # Reverse complement sequences if required
    head_seq = data.data[0] if rc not in ('head', 'both') \
               else reverseComplement(data.data[0])
    tail_seq = data.data[1] if rc not in ('tail', 'both') \
               else reverseComplement(data.data[1])

    # Define result object
    result = SeqResult(data.id, [head_seq, tail_seq])

    # Define stitched sequence annotation
    stitch_ann = OrderedDict([('ID', data.id)])
    if fields_1 is not None:
        head_ann = parseAnnotation(head_seq.description, fields_1,
                                   delimiter=delimiter)
        stitch_ann = mergeAnnotation(stitch_ann, head_ann, delimiter=delimiter)
        result.log['FIELDS1'] = '|'.join(['%s=%s' % (k, v)
                                             for k, v in head_ann.items()])
    if fields_2 is not None:
        tail_ann = parseAnnotation(tail_seq.description, fields_2,
                                   delimiter=delimiter)
        stitch_ann = mergeAnnotation(stitch_ann, tail_ann, delimiter=delimiter)
        result.log['FIELDS2'] = '|'.join(['%s=%s' % (k, v)
                                             for k, v in tail_ann.items()])

    # Assemble sequences
    stitch = assemble_func(head_seq, tail_seq, **assemble_args)
    ab = stitch.head_pos
    xy = stitch.tail_pos
    result.valid = stitch.valid

    # Add reference to log
    if stitch.ref_seq is not None and stitch.ref_pos is not None:
        result.log['REFID'] = stitch.ref_seq.id
        result.log['REFSEQ'] = ' ' * stitch.ref_pos[0] + stitch.ref_seq.seq

    if ab is not None and xy is not None:
        result.log['SEQ1'] = ' ' * xy[0] + head_seq.seq
        result.log['SEQ2'] = ' ' * ab[0] + tail_seq.seq
    else:
        result.log['SEQ1'] = head_seq.seq
        result.log['SEQ2'] = ' ' * (len(head_seq) + (stitch.gap or 0)) + tail_seq.seq

    # Define stitching log
    if stitch.seq is not None:
        # Update stitch annotation
        stitch.seq.id = flattenAnnotation(stitch_ann, delimiter=delimiter)
        stitch.seq.name = stitch.seq.id
        stitch.seq.description = ''
        result.results = stitch.seq
        # Add assembly to log
        result.log['ASSEMBLY'] = stitch.seq.seq
        if 'phred_quality' in stitch.seq.letter_annotations:
            result.log['QUALITY'] = ''.join([chr(q+33) for q in
                                             stitch.seq.letter_annotations['phred_quality']])
        result.log['LENGTH'] = len(stitch)
        result.log['OVERLAP'] = stitch.overlap
    else:
        result.log['ASSEMBLY'] = None

    # Add mode specific log results
    if stitch.gap is not None:
        result.log['GAP'] = stitch.gap
    if stitch.error is not None:
        result.log['ERROR'] = '%.4f' % stitch.error
    if stitch.pvalue is not None:
        result.log['PVALUE'] = '%.4e' % stitch.pvalue
    if stitch.evalue is not None:
        result.log['EVALUE1'] = '%.4e' % stitch.evalue[0]
        result.log['EVALUE2'] = '%.4e' % stitch.evalue[1]
    if stitch.ident is not None:
        result.log['IDENTITY'] = '%.4f' % stitch.ident


    return result


def collectPairQueue(alive, result_queue, collect_queue, result_count,
                     seq_file_1, seq_file_2, out_args):
    """
    Pulls from results queue, assembles results and manages log and file IO

    Arguments: 
    alive = a multiprocessing.Value boolean controlling whether processing 
            continues; when False function returns
    result_queue = a multiprocessing.Queue holding worker results
    collect_queue = a multiprocessing.Queue holding collector return values
    result_count = the number of expected assembled sequences
    seq_file_1 = the first sequence file name
    seq_file_2 = the second sequence file name
    out_args = common output argument dictionary from parseCommonArgs
    
    Returns: 
    None
    (adds a dictionary of {log: log object, out_files: output file names} to collect_queue)
    """
    try:
        # Count records and define output format 
        out_type = getFileType(seq_file_1) if out_args['out_type'] is None \
                   else out_args['out_type']
        
        # Defined valid assembly output handle
        pass_handle = getOutputHandle(seq_file_1, 
                                      'assemble-pass', 
                                      out_dir=out_args['out_dir'], 
                                      out_name=out_args['out_name'], 
                                      out_type=out_type)
        # Defined failed assembly output handles
        if out_args['failed']:
            # Define output name
            if out_args['out_name'] is None:
                out_name_1 = out_name_2 = None
            else:
                out_name_1 = '%s-1' % out_args['out_name']
                out_name_2 = '%s-2' % out_args['out_name']
            fail_handle_1 = getOutputHandle(seq_file_1,
                                            'assemble-fail',
                                            out_dir=out_args['out_dir'],
                                            out_name=out_name_1,
                                            out_type=out_type)
            fail_handle_2 = getOutputHandle(seq_file_2,
                                            'assemble-fail',
                                            out_dir=out_args['out_dir'],
                                            out_name=out_name_2,
                                            out_type=out_type)
        else:
            fail_handle_1 = None
            fail_handle_2 = None

        # Define log handle
        if out_args['log_file'] is None:
            log_handle = None
        else:
            log_handle = open(out_args['log_file'], 'w')
    except:
        alive.value = False
        raise
    
    try:
        # Iterator over results queue until sentinel object reached
        start_time = time()
        iter_count = pass_count = fail_count = 0
        while alive.value:
            # Get result from queue
            if result_queue.empty():  continue
            else:  result = result_queue.get()
            # Exit upon reaching sentinel
            if result is None:  break

            # Print progress for previous iteration
            printProgress(iter_count, result_count, 0.05, start_time)
    
            # Update counts for iteration
            iter_count += 1
    
            # Write log
            printLog(result.log, handle=log_handle)

            # Write assembled sequences
            if result:
                pass_count += 1
                SeqIO.write(result.results, pass_handle, out_type)
            else:
                fail_count += 1
                if fail_handle_1 is not None and fail_handle_2 is not None:
                    SeqIO.write(result.data[0], fail_handle_1, out_type)
                    SeqIO.write(result.data[1], fail_handle_2, out_type)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
        
        # Print total counts
        printProgress(iter_count, result_count, 0.05, start_time)
    
        # Update return values
        log = OrderedDict()
        log['OUTPUT'] = os.path.basename(pass_handle.name)
        log['PAIRS'] = iter_count
        log['PASS'] = pass_count
        log['FAIL'] = fail_count
        collect_dict = {'log':log, 'out_files': [pass_handle.name]}
        collect_queue.put(collect_dict)
        
        # Close file handles
        pass_handle.close()
        if fail_handle_1 is not None:  fail_handle_1.close()
        if fail_handle_2 is not None:  fail_handle_2.close()
        if log_handle is not None:  log_handle.close()
    except:
        alive.value = False
        raise
    
    return None


def assemblePairs(head_file, tail_file, assemble_func, assemble_args={}, 
                  coord_type=default_coord_type, rc=None, 
                  head_fields=None, tail_fields=None,  
                  out_args=default_out_args, nproc=None, queue_size=None):
    """
    Generates consensus sequences

    Arguments: 
    head_file = the head sequence file name
    tail_file = the tail sequence file name
    assemble_func = the function to use to assemble paired ends
    assemble_args = a dictionary of arguments to pass to the assembly function
    coord_type = the sequence header format
    rc = Defines which sequences ('head','tail','both') to reverse complement before assembly;
         if None do not reverse complement sequences
    head_fields = list of annotations in head_file records to copy to assembled record;
                  if None do not copy an annotation
    tail_fields = list of annotations in tail_file records to copy to assembled record;
                  if None do not copy an annotation
    out_args = common output argument dictionary from parseCommonArgs
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc
                 
    Returns: 
    a list of successful output file names
    """
    # Define subcommand label dictionary
    cmd_dict = {alignAssembly:'align', joinSeqPair:'join', referenceAssembly:'reference'}

    # Print parameter info
    log = OrderedDict()
    log['START'] = 'AssemblePairs'
    log['COMMAND'] = cmd_dict.get(assemble_func, assemble_func.__name__)
    log['FILE1'] = os.path.basename(head_file) 
    log['FILE2'] = os.path.basename(tail_file)
    log['COORD_TYPE'] = coord_type
    if 'ref_file' in assemble_args:  log['REFFILE'] = assemble_args['ref_file']
    if 'alpha' in assemble_args:  log['ALPHA'] = assemble_args['alpha']
    if 'max_error' in assemble_args:  log['MAX_ERROR'] = assemble_args['max_error']
    if 'min_len' in assemble_args:  log['MIN_LEN'] = assemble_args['min_len']
    if 'max_len' in assemble_args:  log['MAX_LEN'] = assemble_args['max_len']
    if 'scan_reverse' in assemble_args:  log['SCAN_REVERSE'] = assemble_args['scan_reverse']
    if 'gap' in assemble_args:  log['GAP'] = assemble_args['gap']
    if 'min_ident' in assemble_args:  log['MIN_IDENT'] = assemble_args['min_ident']
    if 'evalue' in assemble_args:  log['EVALUE'] = assemble_args['evalue']
    if 'max_hits' in assemble_args:  log['MAX_HITS'] = assemble_args['max_hits']
    if 'fill' in assemble_args:  log['FILL'] = assemble_args['fill']
    log['NPROC'] = nproc
    printLog(log)

    # Count input files
    head_count = countSeqFile(head_file)
    tail_count = countSeqFile(tail_file)
    if head_count != tail_count:
        sys.exit('Error: FILE1 (n=%i) and FILE2 (n=%i) must have the same number of records' \
                 % (head_count, tail_count))

    # Define feeder function and arguments
    feed_func = feedPairQueue
    # feed_args = {'seq_file_1': head_file,
    #              'seq_file_2': tail_file,
    #              'index_dict': index_dict}
    feed_args = {'seq_file_1': head_file,
                 'seq_file_2': tail_file,
                 'coord_type': coord_type,
                 'delimiter': out_args['delimiter']}
    # Define worker function and arguments
    process_args = {'assemble_func': assemble_func,
                    'assemble_args': assemble_args,
                    'rc': rc,
                    'fields_1': head_fields,
                    'fields_2': tail_fields,
                    'delimiter': out_args['delimiter']}
    work_func = processSeqQueue
    work_args = {'process_func': processAssembly,
                 'process_args': process_args}
    # Define collector function and arguments
    collect_func = collectPairQueue
    # collect_args = {'result_count': pair_count,
    #                 'seq_file_1': head_file,
    #                 'seq_file_2': tail_file,
    #                 'out_args': out_args}
    collect_args = {'result_count': head_count,
                    'seq_file_1': head_file,
                    'seq_file_2': tail_file,
                    'out_args': out_args}

    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func, 
                             feed_args, work_args, collect_args, 
                             nproc, queue_size)
        
    # Print log
    log = OrderedDict()
    log['OUTPUT'] = result['log'].pop('OUTPUT')
    for k, v in result['log'].items():  log[k] = v
    log['END'] = 'AssemblePairs'
    printLog(log)
    
    return result['out_files']
        
        
def getArgParser():
    """
    Defines the ArgumentParser

    Arguments: 
    None
                      
    Returns: 
    an ArgumentParser object
    """
    # Define output file names and header fields
    fields = dedent(
             '''
             output files:
                 assemble-pass
                     successfully assembled reads.
                 assemble-fail
                     raw reads failing paired-end assembly.

             output annotation fields:
                 <user defined>
                     annotation fields specified by the --1f or --2f arguments.
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            formatter_class=CommonHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    subparsers = parser.add_subparsers(title='subcommands', dest='command', metavar='',
                                       help='Assembly method')
    # TODO:  This is a temporary fix for Python issue 9253
    subparsers.required = True

    # Parent parser    
    parser_parent = getCommonArgParser(paired=True, multiproc=True)
    parser_parent.add_argument('--coord', action='store', dest='coord_type', 
                               choices=default_coord_choices, default=default_coord_type,
                               help='The format of the sequence identifier which defines shared coordinate \
                                     information across paired ends')
    parser_parent.add_argument('--rc', action='store', dest='rc', choices=('head', 'tail', 'both'),
                               default=None, help='Specify to reverse complement sequences before stitching')
    parser_parent.add_argument('--1f', nargs='+', action='store', dest='head_fields', type=str, default=None, 
                               help='Specify annotation fields to copy from head records into assembled record')
    parser_parent.add_argument('--2f', nargs='+', action='store', dest='tail_fields', type=str, default=None, 
                               help='Specify annotation fields to copy from tail records into assembled record')
    
    # Paired end overlap alignment mode argument parser
    parser_align = subparsers.add_parser('align', parents=[parser_parent],
                                         formatter_class=CommonHelpFormatter,
                                         help='Assembled pairs by aligning ends')
    parser_align.add_argument('--alpha', action='store', dest='alpha', type=float,
                              default=default_alpha, help='Significance threshold for sequence assemble')
    parser_align.add_argument('--maxerror', action='store', dest='max_error', type=float,
                              default=default_max_error, help='Maximum allowable error rate')
    parser_align.add_argument('--minlen', action='store', dest='min_len', type=int,
                              default=default_min_len, help='Minimum sequence length to scan for overlap')
    parser_align.add_argument('--maxlen', action='store', dest='max_len', type=int,
                              default=default_max_len, help='Maximum sequence length to scan for overlap')
    parser_align.add_argument('--scanrev', action='store_true', dest='scan_reverse',
                              help='''If specified, scan past the end of the tail sequence to allow
                                      the head sequence to overhang the end of the tail sequence.''')
    parser_align.set_defaults(assemble_func=alignAssembly)
    
    # Paired end concatenation mode argument parser
    parser_join = subparsers.add_parser('join', parents=[parser_parent],
                                         formatter_class=CommonHelpFormatter,
                                         help='Assembled pairs by concatenating ends')
    parser_join.add_argument('--gap', action='store', dest='gap', type=int, default=default_gap, 
                             help='Number of N characters to place between ends')
    parser_join.set_defaults(assemble_func=joinSeqPair)    

    # Reference alignment mode argument parser
    parser_ref = subparsers.add_parser('reference', parents=[parser_parent],
                                        formatter_class=CommonHelpFormatter,
                                        help='''Assembled pairs by aligning reads against a
                                             reference database''')
    parser_ref.add_argument('-r', action='store', dest='ref_file', required=True,
                            help='''A FASTA file containing the reference sequence database.''')
    parser_ref.add_argument('--minident', action='store', dest='min_ident', type=float,
                            default=default_min_ident,
                            help='''Minimum identity of the assembled sequence required to call a
                                 valid assembly (between 0 and 1).''')
    parser_ref.add_argument('--evalue', action='store', dest='evalue', type=float,
                            default=default_evalue,
                            help='''Minimum E-value for the ublast reference alignment for both
                                 the head and tail sequence.''')
    parser_ref.add_argument('--maxhits', action='store', dest='max_hits', type=int,
                            default=default_max_hits,
                            help='''Maximum number of hits from ublast to check for matching
                                 head and tail sequence reference alignments.''')
    parser_ref.add_argument('--fill', action='store_true', dest='fill',
                            help='''Specify to insert change the behavior of inserted characters
                                  when the head and tail sequences do not overlap. If specified
                                  this will result in inserted of the V region reference sequence
                                  instead of a sequence of Ns in the non-overlapping region.
                                  Warning, you could end up making chimeric sequences by using
                                  this option.''')
    parser_ref.add_argument('--exec', action='store', dest='usearch_exec',
                            default=default_usearch_exec,
                            help='''The path to the usearch executable file.''')
    parser_ref.set_defaults(assemble_func=referenceAssembly)


    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    # Parse arguments
    parser = getArgParser()
    args = parser.parse_args()
    args_dict = parseCommonArgs(args, in_arg='ref_file')
    
    # Convert case of fields
    if args_dict['head_fields']:  args_dict['head_fields'] = list(map(str.upper, args_dict['head_fields'])) 
    if args_dict['tail_fields']:  args_dict['tail_fields'] = list(map(str.upper, args_dict['tail_fields'])) 
    
    # Define assemble_args dictionary to pass to maskPrimers
    if args_dict['assemble_func'] is alignAssembly:
        args_dict['assemble_args'] = {'alpha': args_dict['alpha'],
                                      'max_error': args_dict['max_error'],
                                      'min_len': args_dict['min_len'],
                                      'max_len': args_dict['max_len'],
                                      'scan_reverse': args_dict['scan_reverse'],
                                      'assembly_stats': AssemblyStats(args_dict['max_len'] + 1)}
        del args_dict['alpha']
        del args_dict['max_error']
        del args_dict['min_len']
        del args_dict['max_len']
        del args_dict['scan_reverse']
    elif args_dict['assemble_func'] is joinSeqPair:
        args_dict['assemble_args'] = {'gap':args_dict['gap']}
        del args_dict['gap']
    elif args_dict['assemble_func'] is referenceAssembly:
        ref_dict = {s.id:s.upper() for s in readSeqFile(args_dict['ref_file'])}
        #ref_file = makeUsearchDb(args_dict['ref_file'], args_dict['usearch_exec'])
        args_dict['assemble_args'] = {'ref_file': args_dict['ref_file'],
                                      'ref_dict': ref_dict,
                                      'min_ident': args_dict['min_ident'],
                                      'evalue': args_dict['evalue'],
                                      'max_hits': args_dict['max_hits'],
                                      'fill': args_dict['fill'],
                                      'usearch_exec': args_dict['usearch_exec']}
        del args_dict['ref_file']
        del args_dict['min_ident']
        del args_dict['evalue']
        del args_dict['max_hits']
        del args_dict['fill']
        del args_dict['usearch_exec']

        # Check if a valid USEARCH executable was specified
        if not os.path.isfile(args.usearch_exec):
            parser.error('%s does not exist' % args.usearch_exec)

    # Call assemblePairs for each sample file
    del args_dict['command']
    del args_dict['seq_files_1']
    del args_dict['seq_files_2']
    for head, tail in zip(args.__dict__['seq_files_1'], 
                          args.__dict__['seq_files_2']):
        args_dict['head_file'] = head
        args_dict['tail_file'] = tail
        assemblePairs(**args_dict)
