#!/usr/bin/env python3
"""
Multiple aligns input sequences by group
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from presto import __version__, __date__

# Imports
import csv
import os
import sys
from argparse import ArgumentParser
from collections import deque, OrderedDict

from textwrap import dedent
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Presto imports
from presto.Defaults import default_delimiter, default_barcode_field, \
                            default_primer_field, default_out_args, default_muscle_exec
from presto.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from presto.Annotation import parseAnnotation
from presto.Applications import runMuscle
from presto.Sequence import calculateDiversity, indexSeqSets
from presto.IO import readPrimerFile, getOutputHandle, printLog
from presto.Multiprocessing import SeqResult, manageProcesses, feedSeqQueue, \
                                   collectSeqQueue


def offsetSeqSet(seq_list, offset_dict, field=default_primer_field, 
                 mode='pad', delimiter=default_delimiter):
    """
    Pads the head of a set of sequences with gaps according to an offset list

    Arguments: 
    seq_list = a list of SeqRecord objects to offset
    offset_dict = a dictionary of {set ID: offset values}
    field = the field in sequence description containing set IDs
    mode = defines the action taken; one of 'pad','cut'
    delimiter = a tuple of delimiters for (annotations, field/values, value lists)
        
    Returns: 
    a MultipleSeqAlignment object containing the alignment
    """
    ann_list = [parseAnnotation(s.description, delimiter=delimiter) for s in seq_list]
    tag_list = [a[field] for a in ann_list]

    # Pad sequences with offsets
    align_list = []
    if mode == 'pad':
        max_len = max([len(s) + offset_dict[t] 
                  for s, t in zip(seq_list, tag_list)])
        for rec, tag in zip(seq_list, tag_list):
            new_rec = rec[:]
            new_rec.letter_annotations = {}
            new_rec.seq = '-' * offset_dict[tag] + new_rec.seq
            new_rec.seq += '-' * (max_len - len(new_rec.seq))
            align_list.append(new_rec)
    # Cut sequences to common start position
    elif mode == 'cut':
        max_offset = max(offset_dict.values())
        cut_dict = {k:(max_offset - v) for k, v in offset_dict.items()}
        max_len = max([len(s) - cut_dict[t] 
                  for s, t in zip(seq_list, tag_list)])
        for rec, tag in zip(seq_list, tag_list):
            new_rec = rec[:]
            new_rec.letter_annotations = {}
            new_rec.seq = new_rec.seq[cut_dict[tag]:]
            new_rec.seq += '-' * (max_len - len(new_rec.seq))
            align_list.append(new_rec)
    else:
        exit('offsettSeqList error:  invalid offset mode')

    # Convert list to MultipleSeqAlignment object
    align = MultipleSeqAlignment(align_list)
    
    return align


def getOffsets(seq_list, align_func=runMuscle, align_args={}, reverse=False):
    """
    Create an offset dictionary for a list of sequences

    Arguments: 
    seq_list = a list of SeqRecord objects
    align_func = the function to use to align sequence sets
    align_args = a dictionary of arguments to pass to align_func
    reverse = if True count tail gaps; if False count head gaps
    
    Returns: 
    a dictionary of {sequence ID: offset value}
    """
    # Perform alignment
    align_list = align_func(seq_list, **align_args)
    
    # Create offset dictionary
    offsets = OrderedDict()
    for aln in align_list:
        # Count non-tail gaps
        seq = str(aln.seq)
        offsets[aln.id] = seq.rstrip('-').count('-') if not reverse \
                          else seq.lstrip('-').count('-')

    # Correct reverse mode offsets
    if reverse:  
        max_offset = max(offsets.values())
        for k in offsets:  offsets[k] = max_offset - offsets[k]
        
    return offsets


def readOffsetFile(offset_file):
    """
    Parses offset file

    Arguments: 
    offset_file = a tab delimited file of set IDs and offset values
    
    Returns: 
    a dictionary of {annotation values: offset values}
    """
    with open(offset_file, 'r') as offset_handle:
        offset_iter = csv.reader(offset_handle, delimiter='\t')
        offset_dict = {r[0]:int(r[1]) for r in offset_iter}
    
    return offset_dict


def writeOffsetFile(primer_file, align_func=runMuscle, align_args={},
                    reverse=False, out_args=default_out_args):
    """
    Generates an offset table from a sequence file

    Arguments: 
    primer_file = name of file containing primer sequences
    align_func = the function to use to align sequence sets
    align_args = a dictionary of arguments to pass to align_func
    reverse = if True count tail gaps; if False count head gaps
    out_args = common output argument dictionary from parseCommonArgs
        
    Returns: 
    the name of the offset output file
    """
    log = OrderedDict()
    log['START'] = 'AlignSets'
    log['COMMAND'] = 'table'
    log['FILE'] = os.path.basename(primer_file)
    log['REVERSE'] = reverse
    printLog(log)
    
    # Read primer file
    primers = readPrimerFile(primer_file)

    # Get offset dictionary
    seq_list = [SeqRecord(Seq(v, IUPAC.ambiguous_dna), id=k) for k, v in primers.items()]
    offset_dict = getOffsets(seq_list, align_func, align_args, reverse)

    # Print log and write offsets to file
    log = OrderedDict()
    for s in seq_list:
        log[s.id] = '%s %i' % (s.seq, offset_dict[s.id])
    printLog(log)

    # Write offset table
    out_tag = 'reverse' if reverse else 'forward'
    with getOutputHandle(primer_file, 'offsets-%s' % out_tag, out_dir=out_args['out_dir'], 
                         out_name=out_args['out_name'], out_type='tab') as out_handle:
        for k, v in offset_dict.items():
            out_handle.write('%s\t%i\n' % (k, v))
    
    # Print final log
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(out_handle.name)
    log['END'] = 'AlignSets'
    printLog(log)
    
    return out_handle.name


def processASQueue(alive, data_queue, result_queue, align_func, align_args={}, 
                      calc_div=False, delimiter=default_delimiter):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments: 
    alive = a multiprocessing.Value boolean controlling whether processing 
            continues; when False function returns
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    align_func = the function to use for alignment
    align_args = a dictionary of optional arguments for the alignment function
    calc_div = if True perform diversity calculation
    delimiter = a tuple of delimiters for (annotations, field/values, value lists)

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
            
            # Define result object
            result = SeqResult(data.id, data.data)
            result.log['BARCODE'] = data.id
            result.log['SEQCOUNT'] = len(data)
    
            # Perform alignment
            seq_list = data.data
            align_list = align_func(seq_list, **align_args)
    
            # Process alignment
            if align_list is not None:
                # Calculate diversity
                if calc_div:
                    diversity = calculateDiversity(align_list)
                    result.log['DIVERSITY'] = diversity
                
                # Restore quality scores
                has_quality = hasattr(seq_list[0], 'letter_annotations') and \
                              'phred_quality' in seq_list[0].letter_annotations
                if has_quality:
                    qual_dict = {seq.id:seq.letter_annotations['phred_quality'] \
                                 for seq in seq_list}
                    for seq in align_list:
                        qual = deque(qual_dict[seq.id])
                        qual_new = [0 if c == '-' else qual.popleft() for c in seq.seq]
                        seq.letter_annotations['phred_quality'] = qual_new
    
                # Add alignment to log
                if 'field' in align_args:
                    for i, seq in enumerate(align_list):
                        ann = parseAnnotation(seq.description, delimiter=delimiter)
                        primer = ann[align_args['field']]
                        result.log['ALIGN%i:%s' % (i + 1, primer)] = seq.seq
                else:
                    for i, seq in enumerate(align_list):  
                        result.log['ALIGN%i' % (i + 1)] = seq.seq
                
                # Add alignment to results
                result.results = align_list
                result.valid = True
                        
            # Feed results to result queue
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


def alignSets(seq_file, align_func, align_args, barcode_field=default_barcode_field,
              calc_div=False, out_args=default_out_args, nproc=None, queue_size=None):
    """
    Performs a multiple alignment on sets of sequences

    Arguments: 
    seq_file = the sample sequence file name
    align_func = the function to use to align sequence sets
    align_args = a dictionary of arguments to pass to align_func
    barcode_field = the annotation containing set IDs
    calc_div = if True calculate average pairwise error for each sequence set
    out_args = common output argument dictionary from parseCommonArgs
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc
                      
    Returns: 
    a tuple of (valid_file, invalid_file) names
    """
    # Define subcommand label dictionary
    cmd_dict = {runMuscle:'align', offsetSeqSet:'offset'}
    
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'AlignSets'
    log['COMMAND'] = cmd_dict[align_func]
    log['FILE'] = os.path.basename(seq_file)
    if 'mode' in align_args: log['MODE'] = align_args['mode']
    log['BARCODE_FIELD'] = barcode_field
    if 'field' in align_args: log['OFFSET_FIELD'] = align_args['field']
    log['CALC_DIV'] = calc_div
    log['NPROC'] = nproc
    printLog(log)
 
    # Define feeder function and arguments
    index_args = {'field': barcode_field, 'delimiter': out_args['delimiter']}
    feed_func = feedSeqQueue
    feed_args = {'seq_file': seq_file,
                 'index_func': indexSeqSets, 
                 'index_args': index_args}
    # Define worker function and arguments
    work_func = processASQueue
    work_args = {'align_func': align_func, 
                 'align_args': align_args,
                 'calc_div': calc_div,
                 'delimiter': out_args['delimiter']}
    # Define collector function and arguments
    collect_func = collectSeqQueue
    collect_args = {'seq_file': seq_file,
                    'task_label': 'align',
                    'out_args': out_args,
                    'index_field': barcode_field}
    
    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func, 
                             feed_args, work_args, collect_args, 
                             nproc, queue_size)
        
    # Print log
    result['log']['END'] = 'AlignSets'
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
             align-pass
                 multiple aligned reads.
             align-fail
                 raw reads failing multiple alignment.
             offsets-forward
                 5\' offset table for input into offset subcommand.
             offsets-reverse
                 3\' offset table for input into offset subcommand.

         output annotation fields:
             None
         ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            formatter_class=CommonHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    subparsers = parser.add_subparsers(title='subcommands', dest='command', metavar='',
                                       help='Alignment method')
    # TODO:  This is a temporary fix for Python issue 9253
    subparsers.required = True

    # Parent parser    
    parser_parent = getCommonArgParser(multiproc=True)
    parser_parent.add_argument('--bf', action='store', dest='barcode_field', type=str,
                               default=default_barcode_field, 
                               help='The annotation field containing barcode labels for sequence grouping')
    parser_parent.add_argument('--div', action='store_true', dest='calc_div',
                               help='Specify to calculate nucleotide diversity of each set (average pairwise error rate)')    
        
    # MUSCLE mode argument parser
    parser_muscle = subparsers.add_parser('muscle', parents=[parser_parent],
                                          formatter_class=CommonHelpFormatter,
                                          help='Align sequence sets using MUSCLE')
    parser_muscle.add_argument('--exec', action='store', dest='muscle_exec', default=default_muscle_exec,
                               help='The location of the MUSCLE executable')
    parser_muscle.set_defaults(align_func=runMuscle)

    # Primer offset mode argument parser
    parser_offset = subparsers.add_parser('offset', parents=[parser_parent],
                                          formatter_class=CommonHelpFormatter,
                                          help='Align sequence sets using predefined 5\' offset')
    parser_offset.add_argument('-d', action='store', dest='offset_table', default=None,
                               help='The tab delimited file of offset tags and values')
    parser_offset.add_argument('--pf', action='store', dest='primer_field', type=str, 
                               default=default_primer_field, 
                               help='The primer field to use for offset assignment')
    parser_offset.add_argument('--mode', action='store', dest='offset_mode', 
                               choices=('pad', 'cut'), default='pad', 
                               help='Specifies whether or align sequence by padding with gaps or by cutting the 5\' \
                                     sequence to a common start position')
    parser_offset.set_defaults(align_func=offsetSeqSet)

    # Offset table generation argument parser
    parser_table = subparsers.add_parser('table', parents=[getCommonArgParser(seq_in=False, seq_out=False, log=False, multiproc=False)],
                                         formatter_class=CommonHelpFormatter,
                                         help='Create a 5\' offset table by primer multiple alignment')
    parser_table.add_argument('-p', action='store', dest='primer_file', required=True, 
                               help='A FASTA or REGEX file containing primer sequences')
    parser_table.add_argument('--reverse', action='store_true', dest='reverse',  
                               help='If specified create a 3\' offset table instead')
    parser_table.add_argument('--exec', action='store', dest='muscle_exec', default=default_muscle_exec,
                               help='The location of the MUSCLE executable')
    parser_table.set_defaults(align_func=runMuscle)
    
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
    if 'barcode_field' in args_dict and args_dict['barcode_field']:
        args_dict['barcode_field'] = args_dict['barcode_field'].upper()
    if 'primer_field' in args_dict and args_dict['primer_field']:
        args_dict['primer_field'] = args_dict['primer_field'].upper()
    
    # Check if a valid MUSCLE executable was specified for muscle mode
    if args.command in ['muscle', 'table'] and not os.path.isfile(args.muscle_exec):
        parser.error('%s does not exist' % args.muscle_exec)
    
    # Define align_args
    if args_dict['align_func'] is runMuscle:
        args_dict['align_args'] = {'muscle_exec':args_dict['muscle_exec']}
        del args_dict['muscle_exec']
    elif args_dict['align_func'] is offsetSeqSet:
        args_dict['align_args'] = {'offset_dict':readOffsetFile(args_dict['offset_table']),
                                   'field':args_dict['primer_field'],
                                   'mode':args_dict['offset_mode']}
        del args_dict['offset_table']
        del args_dict['primer_field']
        del args_dict['offset_mode']
        
    # Call alignSets or createOffsetTable for each input file
    del args_dict['command']
    if args.command in ['muscle', 'offset']: 
        del args_dict['seq_files']
        for f in args.__dict__['seq_files']:
            args_dict['seq_file'] = f
            alignSets(**args_dict)

    elif args.command == 'table':
        writeOffsetFile(**args_dict)
