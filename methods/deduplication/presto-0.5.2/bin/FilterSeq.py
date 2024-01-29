#!/usr/bin/env python3
"""
Filters sequences in FASTA/FASTQ files
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from presto import __version__, __date__

# Imports
import os
import sys
from argparse import ArgumentParser
from collections import OrderedDict
from itertools import groupby
from textwrap import dedent
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Presto imports
from presto.Defaults import default_missing_chars, default_min_qual, default_out_args
from presto.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from presto.IO import getFileType, printLog
from presto.Multiprocessing import SeqResult, manageProcesses, feedSeqQueue, \
                                   processSeqQueue, collectSeqQueue

# Defaults
default_max_missing = 10
default_max_repeat = 15
default_min_length = 250
default_window = 10


def filterLength(data, min_length=default_min_length, inner=True,
                 missing_chars=''.join(default_missing_chars)):
    """
    Filters sequences by length
    
    Arguments: 
    data = a SeqData object with a single SeqRecord to process
    min_length = The minimum length allowed
    inner = if True exclude outer missing characters from calculation
    missing_chars = a string of missing character values
    
    Returns:
    a SeqResult object
    """
    # Get SeqRecord
    seq = data.data

    # Remove outer missing characters if required
    if inner:  
        seq_str = str(seq.seq).strip(missing_chars)
        n = len(seq_str)
    else:
        n = len(seq)
    
    # Build result object
    valid = (n >= min_length)
    result = SeqResult(data.id, seq)
    if valid:
        result.results = seq
        result.valid = True
        
    # Update result log
    result.log['SEQ'] = seq.seq
    result.log['LENGTH'] = n
    
    return result


def filterMissing(data, max_missing=default_max_missing, inner=True,
                  missing_chars=''.join(default_missing_chars)):
    """
    Filters sequences by number of missing nucleotides
    
    Arguments: 
    data = a SeqData object with a single SeqRecord to process
    max_missing = The maximum number of allowed ambiguous characters
    inner = if True exclude outer missing characters from calculation
    missing_chars = a string of missing character values
    
    Returns:
    a SeqResult object
    """
    # Get SeqRecord
    seq = data.data
    seq_str = str(seq.seq)

    # Remove outer missing character if required
    if inner:  seq_str = seq_str.strip(missing_chars)
    # Count missing characters
    n = len([c for c in seq_str if c in missing_chars])
    
    # Build result object
    valid = (n <= max_missing)
    result = SeqResult(data.id, seq)
    if valid:
        result.results = seq
        result.valid = True
    
    # Update result log
    result.log['SEQ'] = seq.seq
    result.log['MISSING'] = n
    
    return result


def filterRepeats(data, max_repeat=default_max_repeat, include_missing=False, inner=True,
                  missing_chars=''.join(default_missing_chars)):
    """
    Filters sequences by fraction of ambiguous nucleotides
    
    Arguments: 
    data = a SeqData object with a single SeqRecord to process
    max_repeat = The maximum number of allowed repeating characters
    include_missing = if True count ambiguous character repeats;
                      if False do not consider ambiguous character repeats
    inner = if True exclude outer missing characters from calculation
    missing_chars = a string of missing character values
    
    Returns:
    a SeqResult object
    """
    # Get SeqRecord
    seq = data.data
    seq_str = str(seq.seq)

    # Remove outer missing character if required
    if inner:  seq_str = seq_str.strip(missing_chars)
    # Remove missing characters if required
    if not include_missing:
        seq_str = ''.join([c for c in seq_str if c not in missing_chars])
    
    groups = ((c, len(list(g))) for c, g in groupby(seq_str))
    __, n = max(groups, key=lambda x: x[1])
    
    # Build result object
    valid = (n <= max_repeat)
    result = SeqResult(data.id, seq)
    if valid:
        result.results = seq
        result.valid = True
    
    # Update result log
    result.log['SEQ'] = seq.seq
    result.log['REPEATS'] = n
    
    return result


def filterQuality(data, min_qual=default_min_qual, inner=True,
                  missing_chars=''.join(default_missing_chars)):
    """
    Filters sequences by quality score
    
    Arguments: 
    data = a SeqData object with a single SeqRecord to process
    min_qual = minimum mean quality score for retained sequences
    inner = if True exclude outer missing characters from calculation
    missing_chars = a string of missing character values
    
    Returns:
    a SeqResult object
    """
    # Get SeqRecord
    seq = data.data

    if inner:  
        seq_str = str(seq.seq)
        seq_cut = seq_str.strip(missing_chars)
        s = seq_str.find(seq_cut)
        quals =  seq.letter_annotations['phred_quality'][s:s + len(seq_cut)]
    else:
        quals = seq.letter_annotations['phred_quality']
    
    q = sum(quals) / len(quals)

    # Build result object
    valid = (q >= min_qual)
    result = SeqResult(data.id, seq)
    if valid:
        result.results = seq
        result.valid = True
    
    # Update result log
    result.log['SEQ'] = seq.seq
    result.log['QUALITY'] = q
    
    return result


def trimQuality(data, min_qual=default_min_qual, window=default_window, reverse=False):
    """
    Cuts sequences using a moving mean quality score
    
    Arguments: 
    data = a SeqData object with a single SeqRecord to process
    min_qual = minimum mean quality to define a cut point
    window = nucleotide window size
    reverse = if True cut the head of the sequence;
              if False cut the tail of the sequence
        
    Returns:
    a SeqResult object
    """
    # Get SeqRecord
    seq = data.data

    quals = seq.letter_annotations['phred_quality']
    # Reverse quality scores if required
    if reverse:  quals = quals[::-1]
    
    # Scan across quality scores for first quality drop-off
    end = len(quals)
    for s in range(0, end, window):
        q_win = quals[s:s + window]
        q = sum(q_win) / len(q_win)
        if q < min_qual:
            end = s
            break

    # Define trimmed sequence
    if not reverse:
        trim_seq = seq[:end]
        out_str = str(trim_seq.seq)
    else:
        trim_seq = seq[len(seq) - end:]
        out_str =  ' ' * (len(seq) - end) + str(trim_seq.seq)
        
    # Build result object
    valid = (len(trim_seq) > 0)
    result = SeqResult(data.id, seq)
    if valid:
        result.results = trim_seq
        result.valid = True
    
    # Update result log
    result.log['INSEQ'] = seq.seq
    result.log['OUTSEQ'] = out_str
    result.log['LENGTH'] = len(trim_seq)
    
    return result


def maskQuality(data, min_qual=default_min_qual):
    """
    Masks characters by in sequence by quality score
    
    Arguments: 
    data = a SeqData object with a single SeqRecord to process
    min_qual = minimum quality for retained characters
        
    Returns:
    a SeqResult object
    """
    # Get SeqRecord
    seq = data.data
    seq_str = str(seq.seq)
    quals = seq.letter_annotations['phred_quality']

    # Mask low quality nucleotides
    mask_chars = [seq_str[i] if q >= min_qual else 'N' for i, q in enumerate(quals)]
    mask_count = sum(1 for q in quals if q < min_qual)

    # Define masked SeqRecord
    mask_seq = SeqRecord(Seq(''.join(mask_chars), IUPAC.ambiguous_dna), 
                         id=seq.id, 
                         name=seq.name, 
                         description=seq.description,
                         letter_annotations=seq.letter_annotations)
    
    # Build result object
    result = SeqResult(data.id, seq)
    result.results = mask_seq
    result.valid = True
    
    # Update result log
    result.log['INSEQ'] = seq.seq
    result.log['OUTSEQ'] = mask_seq.seq
    result.log['MASKED'] = mask_count
    
    return result


def filterSeq(seq_file, filter_func, filter_args={}, out_args=default_out_args, 
              nproc=None, queue_size=None):
    """
    Filters sequences by fraction of ambiguous nucleotides
    
    Arguments: 
    seq_file = the sequence file to filter
    filter_func = the function to use for filtering sequences
    filter_args = a dictionary of arguments to pass to filter_func
    out_args = common output argument dictionary from parseCommonArgs
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc
                 
    Returns:
    a list of successful output file names
    """
    # Define output file label dictionary
    cmd_dict = {filterLength:'length', filterMissing:'missing', 
                filterRepeats:'repeats', filterQuality:'quality', 
                maskQuality:'maskqual', trimQuality:'trimqual'}
    
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'FilterSeq'
    log['COMMAND'] = cmd_dict.get(filter_func, filter_func.__name__)
    log['FILE'] = os.path.basename(seq_file)
    for k in sorted(filter_args):  log[k.upper()] = filter_args[k]
    log['NPROC'] = nproc
    printLog(log)
    
    # Check input type
    in_type = getFileType(seq_file)
    if in_type != 'fastq' and filter_func in (filterQuality, maskQuality, trimQuality):
        sys.exit('ERROR:  Input file must be FASTQ for %s mode' % cmd_dict[filter_func])
    
    # Define feeder function and arguments
    feed_func = feedSeqQueue
    feed_args = {'seq_file': seq_file}
    # Define worker function and arguments
    work_func = processSeqQueue
    work_args = {'process_func': filter_func, 
                 'process_args': filter_args}
    # Define collector function and arguments
    collect_func = collectSeqQueue
    collect_args = {'seq_file': seq_file,
                    'task_label': cmd_dict[filter_func],
                    'out_args': out_args}
    
    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func, 
                             feed_args, work_args, collect_args, 
                             nproc, queue_size)
        
    # Print log
    result['log']['END'] = 'FilterSeq'
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
                 <command>-pass
                     reads passing filtering operation and modified accordingly, where
                     <command> is the name of the filtering operation that was run.
                 <command>-fail
                     raw reads failing filtering criteria, where <command> is the name of
                     the filtering operation.

             output annotation fields:
                 None
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            formatter_class=CommonHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    subparsers = parser.add_subparsers(title='subcommands', metavar='',
                                       help='Filtering operation')
    # TODO:  This is a temporary fix for Python issue 9253
    subparsers.required = True

    # Parent parser
    parser_parent = getCommonArgParser(annotation=False, log=True, multiproc=True)
    
    # Length filter mode argument parser
    parser_length = subparsers.add_parser('length', parents=[parser_parent],
                                          formatter_class=CommonHelpFormatter, 
                                          help='Sequence length filtering mode')
    parser_length.add_argument('-n', action='store', dest='min_length', type=int, 
                               default=default_min_length, 
                               help='Minimum sequence length to retain.')
    parser_length.add_argument('--inner', action='store_true', dest='inner',
                               help='''If specified exclude consecutive missing characters
                                    at either end of the sequence.''')
    parser_length.set_defaults(filter_func=filterLength)
    
    # Missing character filter mode argument parser
    parser_missing = subparsers.add_parser('missing', parents=[parser_parent],
                                           formatter_class=CommonHelpFormatter, 
                                           help='Missing nucleotide filtering mode')
    parser_missing.add_argument('-n', action='store', dest='max_missing', type=int, 
                                default=default_max_missing, 
                                help='Threshold for fraction of gap or N nucleotides.')
    parser_missing.add_argument('--inner', action='store_true', dest='inner',
                                help='''If specified exclude consecutive missing characters
                                     at either end of the sequence.''')
    parser_missing.set_defaults(filter_func=filterMissing)
    
    # Continuous repeating character filter mode argument parser
    parser_repeats = subparsers.add_parser('repeats', parents=[parser_parent],
                                           formatter_class=CommonHelpFormatter, 
                                           help='Consecutive nucleotide repeating filtering mode')
    parser_repeats.add_argument('-n', action='store', dest='max_repeat', type=int, 
                                default=default_max_repeat, 
                                help='Threshold for fraction of repeating nucleotides.')
    parser_repeats.add_argument('--missing', action='store_true', dest='include_missing',
                                help='''If specified count consecutive gap and N characters '
                                     in addition to {A,C,G,T}.''')
    parser_repeats.add_argument('--inner', action='store_true', dest='inner',
                                help='''If specified exclude consecutive missing characters
                                     at either end of the sequence.''')
    parser_repeats.set_defaults(filter_func=filterRepeats)
    
    # Quality filter mode argument parser
    parser_quality = subparsers.add_parser('quality', parents=[parser_parent],
                                          formatter_class=CommonHelpFormatter, 
                                          help='Quality filtering mode')
    parser_quality.add_argument('-q', action='store', dest='min_qual', type=float, 
                                default=default_min_qual, help='Quality score threshold.')
    parser_quality.add_argument('--inner', action='store_true', dest='inner',
                                help='''If specified exclude consecutive missing characters
                                     at either end of the sequence.''')
    parser_quality.set_defaults(filter_func=filterQuality)

    # Mask mode argument parser
    parser_maskqual = subparsers.add_parser('maskqual', parents=[parser_parent], 
                                        formatter_class=CommonHelpFormatter,
                                        help='Character masking mode')
    parser_maskqual.add_argument('-q', action='store', dest='min_qual', type=float, 
                             default=default_min_qual, help='Quality score threshold.')
    parser_maskqual.set_defaults(filter_func=maskQuality)

    # Trim mode argument parser
    parser_trimqual = subparsers.add_parser('trimqual', parents=[parser_parent], 
                                            formatter_class=CommonHelpFormatter,
                                            help='Sequence trimming mode')
    parser_trimqual.add_argument('-q', action='store', dest='min_qual', type=float, 
                                 default=default_min_qual, help='Quality score threshold.')
    parser_trimqual.add_argument('--win', action='store', dest='window', type=int, 
                                 default=default_window, 
                                 help='Nucleotide window size for moving average calculation.')
    parser_trimqual.add_argument('--reverse', action='store_true', dest='reverse', 
                                 help='''Specify to trim the head of the sequence rather
                                      than the tail.''')
    parser_trimqual.set_defaults(filter_func=trimQuality)
    
    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    # Parse arguments
    parser = getArgParser()
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    
    # Create filter_args
    filter_keys = ['min_qual', 'max_repeat', 'max_missing', 'min_length', 'inner', 
                   'include_missing', 'window', 'reverse']
    args_dict['filter_args'] = dict((k, args_dict[k]) for k in args_dict if k in filter_keys)
    for k in args_dict['filter_args']:  del args_dict[k]
    
    # Calls quality processing function
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        filterSeq(**args_dict)
