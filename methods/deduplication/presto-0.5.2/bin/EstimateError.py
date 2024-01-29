#!/usr/bin/env python3
"""
Calculates annotation set error rates
"""
# Info
__author__ = 'Jason Anthony Vander Heiden, Namita Gupta'
from presto import __version__, __date__

# Imports
import os
import sys
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from collections import OrderedDict
from itertools import permutations
from textwrap import dedent
from time import time

# Presto imports
from presto.Defaults import default_barcode_field, default_missing_chars, \
                            default_min_freq, default_min_qual, default_out_args
from presto.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from presto.IO import getFileType, countSeqSets, getOutputHandle, printLog, printProgress
from presto.Sequence import getDNAScoreDict, calculateDiversity, qualityConsensus, \
                            frequencyConsensus, indexSeqSets
from presto.Multiprocessing import SeqResult, manageProcesses, feedSeqQueue

# Defaults
default_min_count = 10
        

def countMismatches(seq_list, ref_seq, ignore_chars=default_missing_chars, 
                    score_dict=getDNAScoreDict(mask_score=(1, 1), gap_score=(1, 1))):
    """
    Counts the occurrence of nucleotide mismatches in a set of sequences

    Arguments: 
    seq_list = a list of SeqRecord objects with aligned sequences
    ref_seq = a SeqRecord object containing the reference sequence to match against
    ignore_chars = list of characters to exclude from mismatch counts
    score_dict = optional dictionary of alignment scores as {(char1, char2): score}

    Returns: 
    a dictionary of pandas.DataFrame objects containing [mismatch, qsum, total] counts  
    for {pos:sequence position, nuc:nucleotide pairs, qual:quality score, set:sequence set} 
    """
    # Define position mismatch DataFrame
    pos_max = max([len(s) for s in seq_list])
    pos_df = pd.DataFrame(0, index=list(range(pos_max)), 
                          columns=['mismatch', 'q_sum', 'total'], dtype=float)
    # Define nucleotide mismatch DataFrame
    nuc_pairs = list(permutations(['A', 'C', 'G', 'T'], 2))
    nuc_df = pd.DataFrame(0, index=pd.MultiIndex.from_tuples(nuc_pairs, names=['obs', 'ref']), 
                          columns=['mismatch', 'q_sum', 'total'], dtype=float)
    # Define quality mismatch DataFrame
    qual_df = pd.DataFrame(0, index=list(range(94)), 
                           columns=['mismatch', 'q_sum', 'total'], dtype=float)    

    # Iterate over seq_list and count mismatches
    for seq in seq_list:
        qual = seq.letter_annotations['phred_quality']
        for i, b in enumerate(seq):
            a = ref_seq[i]
            q = qual[i]
            # Update total counts and qualities
            if {a, b}.isdisjoint(ignore_chars):  
                pos_df.ix[i, 'total'] += 1
                pos_df.ix[i, 'q_sum'] += q
                nuc_df.ix[b]['total'] += 1
                nuc_df.ix[b]['q_sum'] += q
                qual_df.ix[q, 'total'] += 1
                qual_df.ix[q, 'q_sum'] += q
            # Update mismatch counts
            if score_dict[(a, b)] == 0:
                pos_df.ix[i, 'mismatch'] += 1
                nuc_df.ix[(a, b), 'mismatch'] += 1
                qual_df.ix[q, 'mismatch'] += 1

    # Define set total mismatch DataFrame
    set_df = pd.DataFrame([pos_df.sum(axis=0)], index=[len(seq_list)], 
                          columns=['mismatch', 'q_sum', 'total'], dtype=float)

    return {'pos':pos_df, 'nuc':nuc_df, 'qual':qual_df, 'set':set_df}

    
def processEEQueue(alive, data_queue, result_queue, cons_func, cons_args={}, 
                   min_count=default_min_count, max_diversity=None):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments: 
    alive = a multiprocessing.Value boolean controlling whether processing 
            continues; when False function returns
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    cons_func = the function to use for consensus generation 
    cons_args = a dictionary of optional arguments for the consensus function
    min_count = threshold number of sequences to retain a set
    max_diversity = the minimum diversity score to retain a set;
                    if None do not calculate diversity
                        
    Returns: 
    None
    """
    try:
        # Iterator over data queue until sentinel object reached
        while alive.value:
            # Get data from queue
            if data_queue.empty():  continue
            else:  data = data_queue.get()
            # Exit upon reaching sentinel
            if data is None:  break
            
            # Define result dictionary for iteration
            result = SeqResult(data.id, data.data)
            result.results = {'pos':None, 
                              'nuc':None, 
                              'qual':None, 
                              'set':None}
            
            # Define sequences set
            seq_list = data.data
            seq_count = len(data)
            
            # Update log
            result.log['SET'] = data.id
            result.log['SEQCOUNT'] = seq_count
            for i, s in enumerate(seq_list):
                result.log['SEQ%i' % (i + 1)] = str(s.seq)
            
            # Check count threshold and continue if failed
            if len(data) < min_count:
                result_queue.put(result)
                continue
                
            # Calculate average pairwise error rate
            if max_diversity is not None:
                diversity = calculateDiversity(seq_list)
                result.log['DIVERSITY'] = diversity
                # Check diversity threshold and continue if failed
                if diversity > max_diversity:
                    result_queue.put(result)
                    continue
                
            # Define reference sequence by consensus
            ref_seq = cons_func(seq_list, **cons_args)
            
            # Count mismatches against consensus
            mismatch = countMismatches(seq_list, ref_seq)
            
            # Calculate average reported and observed error
            reported_q = mismatch['set']['q_sum'].sum() / mismatch['set']['total'].sum()
            error_rate = mismatch['set']['mismatch'].sum() / mismatch['set']['total'].sum()
    
            # Update log
            result.log['REFERENCE'] = str(ref_seq.seq)
            result.log['MISMATCH'] = ''.join(['*' if x > 0 else ' ' \
                                              for x in mismatch['pos']['mismatch']])
            result.log['ERROR'] = '%.6f' % error_rate
            result.log['REPORTED_Q'] = '%.2f' % reported_q
            result.log['EMPIRICAL_Q'] = '%.2f' % (-10 * np.log10(error_rate))
                
            # Update results and feed result queue
            result.valid = True
            result.results.update(mismatch)
            result_queue.put(result)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        sys.stderr.write('Error processing sequence set with ID: %s.\n' % data.id)
        raise

    return None


def collectEEQueue(alive, result_queue, collect_queue, seq_file, out_args, set_field):
    """
    Pulls from results queue, assembles results and manages log and file IO

    Arguments: 
    alive = a multiprocessing.Value boolean controlling whether processing 
            continues; when False function returns
    result_queue = a multiprocessing.Queue holding worker results
    collect_queue = a multiprocessing.Queue to store collector return values
    seq_file = the sample sequence file name
    out_args = common output argument dictionary from parseCommonArgs
    set_field = the field defining set membership

    Returns:
    None
    (adds a dictionary of {log: log object, out_files: output file names} to collect_queue)
    """
    try:
        # Count sets
        result_count = countSeqSets(seq_file, set_field, out_args['delimiter'])
        
        # Define empty DataFrames to store assembled results
        pos_df = pd.DataFrame(None, columns=['mismatch', 'q_sum', 'total'], dtype=float)
        qual_df = pd.DataFrame(None, columns=['mismatch', 'q_sum', 'total'], dtype=float)
        nuc_df = pd.DataFrame(None, columns=['mismatch', 'q_sum', 'total'], dtype=float)
        set_df = pd.DataFrame(None, columns=['mismatch', 'q_sum', 'total'], dtype=float)
        
        # Open log file
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
        set_count = seq_count = pass_count = fail_count = 0
        while alive.value:
            # Get result from queue
            if result_queue.empty():  continue
            else:  result = result_queue.get()
            # Exit upon reaching sentinel
            if result is None:  break

            # Print progress for previous iteration
            printProgress(set_count, result_count, 0.05, start_time)
            
            # Update counts for iteration
            set_count += 1
            seq_count += result.data_count
            
            # Sum results
            if result:
                pass_count += 1
                pos_df = pos_df.add(result.results['pos'], fill_value=0)
                qual_df = qual_df.add(result.results['qual'], fill_value=0)
                nuc_df = nuc_df.add(result.results['nuc'], fill_value=0)
                set_df = set_df.add(result.results['set'], fill_value=0)
            else:
                fail_count += 1
                
            # Write log
            printLog(result.log, handle=log_handle)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
        
        # Print final progress
        printProgress(set_count, result_count, 0.05, start_time)
        
        # Generate log
        log = OrderedDict()
        for i in range(4): 
            log['OUTPUT%i' % (i + 1)] = None
        log['SETS'] = set_count
        log['SEQUENCES'] = seq_count
        log['PASS'] = pass_count
        log['FAIL'] = fail_count
        log['POSITION_ERROR'] = None 
        log['NUCLEOTIDE_ERROR'] = None
        log['QUALITY_ERROR'] = None
        log['SET_ERROR'] = None
     
        # Return if no mismatch data
        if pass_count == 0:
            collect_dict = {'log':log, 'out_files': None}
            collect_queue.put(collect_dict)
            return None
    
        # Calculate error rates
        pos_df['error'] = pos_df['mismatch'] / pos_df['total'] 
        nuc_df['error'] = nuc_df['mismatch'] / nuc_df['total']
        qual_df['error'] = qual_df['mismatch'] / qual_df['total']
        set_df['error'] = set_df['mismatch'] / set_df['total']
        
        # Convert error to empirical quality score
        pos_df['emp_q'] = -10 * np.log10(pos_df['error'])
        nuc_df['emp_q'] = -10 * np.log10(nuc_df['error'])
        qual_df['emp_q'] = -10 * np.log10(qual_df['error'])
        set_df['emp_q'] = -10 * np.log10(set_df['error'])
    
        # Calculate reported quality means
        pos_df['rep_q'] = pos_df['q_sum'] / pos_df['total'] 
        nuc_df['rep_q'] = nuc_df['q_sum'] / nuc_df['total']
        qual_df['rep_q'] = qual_df['q_sum'] / qual_df['total']
        set_df['rep_q'] = set_df['q_sum'] / set_df['total']
            
        # Calculate overall error rate
        pos_error = pos_df['mismatch'].sum() / pos_df['total'].sum() 
        qual_error = qual_df['mismatch'].sum() / qual_df['total'].sum() 
        nuc_error = nuc_df['mismatch'].sum() / nuc_df['total'].groupby(level='obs').mean().sum()
        set_error = set_df['mismatch'].sum() / set_df['total'].sum() 
    
        # Build results dictionary
        assembled = {'pos':pos_df, 'qual':qual_df, 'nuc':nuc_df, 'set':set_df}
        
        # Write assembled error counts to output files
        out_files = writeResults(assembled, seq_file, out_args)
        
        # Update log
        for i, f in enumerate(out_files): 
            log['OUTPUT%i' % (i + 1)] = os.path.basename(f)
        log['POSITION_ERROR'] = pos_error 
        log['NUCLEOTIDE_ERROR'] = nuc_error 
        log['QUALITY_ERROR'] = qual_error
        log['SET_ERROR'] = set_error
        
        # Update collector results
        collect_dict = {'log':log, 'out_files': out_files}
        collect_queue.put(collect_dict)
    except:
        alive.value = False
        raise
    
    return None


def writeResults(results, seq_file, out_args):
    """
    Formats results and writes to output files

    Arguments: 
    results = assembled results dictionary
    seq_file = the sample sequence file name
    out_args = common output argument dictionary from parseCommonArgs

    Returns:
    a tuple of (position error, nucleotide pairwise error, quality error, sequence set) file names
    """
    pos_df = results['pos']
    nuc_df = results['nuc']
    qual_df = results['qual']
    set_df = results['set']

    # Type conversion to int of mismatch and total columns
    pos_df[['mismatch', 'total']] = pos_df[['mismatch', 'total']].astype(int) 
    nuc_df[['mismatch', 'total']] = nuc_df[['mismatch', 'total']].astype(int) 
    qual_df[['mismatch', 'total']] = qual_df[['mismatch', 'total']].astype(int) 
    set_df[['mismatch', 'total']] = set_df[['mismatch', 'total']].astype(int) 
    
    # Write to tab delimited files
    file_args = {'out_dir':out_args['out_dir'], 'out_name':out_args['out_name'], 'out_type':'tab'}
    with getOutputHandle(seq_file, 'error-position', **file_args) as pos_handle, \
            getOutputHandle(seq_file, 'error-quality', **file_args) as qual_handle, \
            getOutputHandle(seq_file, 'error-nucleotide', **file_args) as nuc_handle, \
            getOutputHandle(seq_file, 'error-set', **file_args) as set_handle:
        pos_df.to_csv(pos_handle, sep='\t', na_rep='NA', index_label='POSITION', 
                      cols=['rep_q', 'mismatch', 'total', 'error', 'emp_q'],
                      header=['REPORTED_Q', 'MISMATCHES', 'OBSERVATIONS', 'ERROR', 'EMPIRICAL_Q'], float_format='%.6f')
        nuc_df.to_csv(nuc_handle, sep='\t', na_rep='NA', index_label=['OBSERVED', 'REFERENCE'], 
                      cols=['rep_q', 'mismatch', 'total', 'error', 'emp_q'],
                      header=['REPORTED_Q', 'MISMATCHES', 'OBSERVATIONS', 'ERROR', 'EMPIRICAL_Q'], float_format='%.6f')
        qual_df.to_csv(qual_handle, sep='\t', na_rep='NA', index_label='Q',
                       cols=['rep_q', 'mismatch', 'total', 'error', 'emp_q'], 
                       header=['REPORTED_Q', 'MISMATCHES', 'OBSERVATIONS', 'ERROR', 'EMPIRICAL_Q'], float_format='%.6f')
        set_df.to_csv(set_handle, sep='\t', na_rep='NA', index_label='SET_COUNT', 
                      cols=['rep_q', 'mismatch', 'total', 'error', 'emp_q'],
                      header=['REPORTED_Q', 'MISMATCHES', 'OBSERVATIONS', 'ERROR', 'EMPIRICAL_Q'], float_format='%.6f')

    return (pos_handle.name, qual_handle.name, nuc_handle.name, set_handle.name)


def estimateError(seq_file, cons_func=frequencyConsensus, cons_args={}, 
                  set_field=default_barcode_field, min_count=default_min_count, max_diversity=None, 
                  out_args=default_out_args, nproc=None, queue_size=None):
    """
    Calculates error rates of sequence sets

    Arguments: 
    seq_file = the sample sequence file name
    cons_func = the function to use for consensus generation 
    cons_args = a dictionary of arguments for the consensus function
    set_field = the annotation field containing set IDs
    min_count = threshold number of sequences to consider a set
    max_diversity = a threshold defining the average pairwise error rate required to retain a read group;
                    if None do not calculate diversity
    out_args = common output argument dictionary from parseCommonArgs
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc
                    
    Returns: 
    a list of tuples of (position error, quality error, nucleotide pairwise error) output file names
    """
    # Define subcommand label dictionary
    cmd_dict = {frequencyConsensus:'freq', qualityConsensus:'qual'}
    
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'EstimateError'
    log['FILE'] = os.path.basename(seq_file)
    log['MODE'] = cmd_dict.get(cons_func, cons_func.__name__)
    log['SET_FIELD'] = set_field
    log['MIN_COUNT'] = min_count
    log['MAX_DIVERSITY'] = max_diversity
    log['NPROC'] = nproc
    printLog(log)
    
    # Check input file type
    in_type = getFileType(seq_file)
    if in_type != 'fastq':  sys.exit('ERROR:  Input file must be FASTQ')
    
        # Define feeder function and arguments
    index_args = {'field': set_field, 'delimiter': out_args['delimiter']}
    feed_func = feedSeqQueue
    feed_args = {'seq_file': seq_file,
                 'index_func': indexSeqSets, 
                 'index_args': index_args}
    # Define worker function and arguments
    work_func = processEEQueue
    work_args = {'cons_func': cons_func, 
                 'cons_args': cons_args,
                 'min_count': min_count,
                 'max_diversity': max_diversity}
    # Define collector function and arguments
    collect_func = collectEEQueue
    collect_args = {'seq_file': seq_file,
                    'out_args': out_args,
                    'set_field': set_field}

    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func, 
                             feed_args, work_args, collect_args, 
                             nproc, queue_size)
        
    # Print log
    result['log']['END'] = 'EstimateError'
    printLog(result['log'])
        
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
                 error-position
                     estimated error by read position.
                 error-quality
                     estimated error by the quality score assigned within the input file.
                 error-nucleotide
                     estimated error by nucleotide.
                 error-set
                     estimated error by barcode read group size.

             output fields:
                 POSITION
                     read position with base zero indexing.
                 Q
                     Phred quality score.
                 OBSERVED
                     observed nucleotide value.
                 REFERENCE
                     consensus nucleotide for the barcode read group.
                 SET_COUNT
                     barcode read group size.
                 REPORTED_Q
                     mean Phred quality score reported within the input file for the given
                     position, quality score, nucleotide or read group.
                 MISMATCHES
                     count of observed mismatches from consensus for the given position,
                     quality score, nucleotide or read group.
                 OBSERVATIONS
                     total count of observed values for each position, quality score,
                     nucleotide or read group size.
                 ERROR
                     estimated error rate.
                 EMPIRICAL_Q
                     estimated error rate converted to a Phred quality score.
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[getCommonArgParser(seq_out=False,
                                                        failed=False,
                                                        multiproc=True)],
                            formatter_class=CommonHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    
    parser.add_argument('-f', action='store', dest='set_field', type=str, default=default_barcode_field, 
                        help='The name of the annotation field to group sequences by')
    parser.add_argument('-n', action='store', dest='min_count', type=int, default=default_min_count,
                        help='The minimum number of sequences needed to consider a set')
    parser.add_argument('--mode', action='store', dest='mode', choices=('freq', 'qual'), default='freq', 
                        help='''Specifies which method to use to determine the consensus
                             sequence. The "freq" method will determine the consensus by
                             nucleotide frequency at each position and assign the most
                             common value.  The "qual" method will weight values by their
                             quality scores to determine the consensus nucleotide at
                             each position.''')
    parser.add_argument('-q', action='store', dest='min_qual', type=float, default=default_min_qual,
                        help='''Consensus quality score cut-off under which an ambiguous '
                             character is assigned.''')
    parser.add_argument('--freq', action='store', dest='min_freq', type=float, default=default_min_freq,
                        help='''Fraction of character occurrences under which an ambiguous
                             character is assigned.''')
    parser.add_argument('--maxdiv', action='store', dest='max_diversity', type=float, default=None,
                        help='''Specify to calculate the nucleotide diversity of each read
                              group (average pairwise error rate) and exclude groups which
                              exceed the given diversity threshold.''')
    
    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    # Parse arguments
    parser = getArgParser()
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    
    # Convert case of fields
    if args_dict['set_field']:  args_dict['set_field'] = args_dict['set_field'].upper()
    
    # Define cons_func and cons_args
    if args_dict['mode'] == 'freq':
        args_dict['cons_func'] = frequencyConsensus
        args_dict['cons_args'] = {'min_freq':args_dict['min_freq']}
    elif args_dict['mode'] == 'qual':
        args_dict['cons_func'] = qualityConsensus
        args_dict['cons_args'] = {'min_qual':args_dict['min_qual'], 
                                  'min_freq':args_dict['min_freq'], 
                                  'dependent':False}
    del args_dict['mode']
    if 'min_freq' in args_dict:  del args_dict['min_freq']
    if 'min_qual' in args_dict:  del args_dict['min_qual']
    
    # Call estimateError for each sample file    
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        estimateError(**args_dict)
