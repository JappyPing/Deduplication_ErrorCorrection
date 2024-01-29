#!/usr/bin/env python3
"""
Builds a consensus sequence for each set of input sequences
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from presto import __version__, __date__

# Imports
import os
import sys
from argparse import ArgumentParser
from collections import OrderedDict

from textwrap import dedent

# Presto imports
from presto.Defaults import default_delimiter, default_barcode_field, \
    default_min_freq, default_out_args
from presto.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from presto.Annotation import flattenAnnotation, mergeAnnotation, getAnnotationValues, \
                              annotationConsensus
from presto.IO import getFileType, printLog
from presto.Sequence import subsetSeqSet, calculateDiversity, \
                            qualityConsensus, frequencyConsensus, indexSeqSets, \
                            calculateSetError, deleteSeqPositions, findGapPositions
from presto.Multiprocessing import SeqResult, manageProcesses, feedSeqQueue, \
                                   collectSeqQueue

# Defaults
default_min_count = 1
default_min_qual = 0


def processBCQueue(alive, data_queue, result_queue, cons_func, cons_args={},
                   min_count=default_min_count, primer_field=None, primer_freq=None,
                   max_gap=None, max_error=None, max_diversity=None,
                   copy_fields=None, copy_actions=None, delimiter=default_delimiter):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments:
    alive = a multiprocessing.Value boolean controlling whether processing 
            continues; when False function returns
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    cons_func = the function to use for consensus generation 
    cons_args = a dictionary of optional arguments for the consensus function
    min_count = threshold number of sequences to define a consensus 
    primer_field = the annotation field containing primer names;
                   if None do not annotate with primer names
    primer_freq = the maximum primer frequency that must be meet to build a consensus;
                  if None do not filter by primer frequency
    max_gap = the maximum frequency of (., -) characters allowed before
              deleting a position; if None do not delete positions
    max_error = the minimum error rate to retain a set;
                if None do not calculate error rate
    max_diversity = a threshold defining the average pairwise error rate required to retain a read group;
                    if None do not calculate diversity
    copy_fields = a list of annotations to copy into consensus sequence annotations;
                  if None no additional annotations will be copied
    copy_actions = the list of actions to take for each copy_fields;
                   one of ['set', 'majority', 'min', 'max', 'sum']
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
            
            # Define result dictionary for iteration
            result = SeqResult(data.id, data.data)
            result.log['BARCODE'] = data.id
            result.log['SEQCOUNT'] = len(data)
    
            # Define primer annotations and consensus primer if applicable
            if primer_field is None:
                primer_ann = None
                seq_list = data.data
            else:
                # Calculate consensus primer            
                primer_ann = OrderedDict()
                prcons = annotationConsensus(data.data, primer_field, delimiter=delimiter)
                result.log['PRIMER'] = ','.join(prcons['set'])
                result.log['PRCOUNT'] = ','.join([str(c) for c in prcons['count']])
                result.log['PRCONS'] = prcons['cons']
                result.log['PRFREQ'] = prcons['freq']
                if primer_freq is None:
                    # Retain full sequence set if not in primer consensus mode
                    seq_list = data.data
                    primer_ann = mergeAnnotation(primer_ann, {'PRIMER':prcons['set']}, 
                                                 delimiter=delimiter)
                    primer_ann = mergeAnnotation(primer_ann, {'PRCOUNT':prcons['count']}, 
                                                 delimiter=delimiter)
                elif prcons['freq'] >= primer_freq:
                    # Define consensus subset
                    seq_list = subsetSeqSet(data.data, primer_field, prcons['cons'], 
                                            delimiter=delimiter)
                    primer_ann = mergeAnnotation(primer_ann, {'PRCONS':prcons['cons']}, 
                                                 delimiter=delimiter)
                    primer_ann = mergeAnnotation(primer_ann, {'PRFREQ':prcons['freq']}, 
                                                 delimiter=delimiter)
                else:
                    # If set fails primer consensus, feed result queue and continue
                    result_queue.put(result)
                    continue
    
            # Check count threshold
            cons_count = len(seq_list)
            result.log['CONSCOUNT'] = cons_count
            if cons_count < min_count:
                # If set fails count threshold, feed result queue and continue
                result_queue.put(result)
                continue

            # Update log with input sequences
            for i, s in enumerate(seq_list):
                result.log['INSEQ%i' % (i + 1)] = str(s.seq)

            # If primer and count filters pass, generate consensus sequence
            consensus = cons_func(seq_list, **cons_args)

            # Delete positions with gap frequency over max_gap and update log with consensus
            if max_gap is not None:
                gap_positions = set(findGapPositions(seq_list, max_gap))
                result.log['CONSENSUS'] = ''.join([' ' if i in gap_positions else x \
                                                   for i, x in enumerate(consensus.seq)])
                if 'phred_quality' in consensus.letter_annotations:
                    result.log['QUALITY'] = ''.join([' ' if i in gap_positions else chr(q + 33) \
                                                     for i, q in enumerate(consensus.letter_annotations['phred_quality'])])
                consensus = deleteSeqPositions(consensus, gap_positions)
            else:
                gap_positions = None
                result.log['CONSENSUS'] = str(consensus.seq)
                if 'phred_quality' in consensus.letter_annotations:
                    result.log['QUALITY'] = ''.join([chr(q + 33) for q in consensus.letter_annotations['phred_quality']])

            # Calculate nucleotide diversity
            if max_diversity is not None:
                diversity = calculateDiversity(seq_list)
                result.log['DIVERSITY'] = diversity

                # If diversity exceeds threshold, feed result queue and continue
                if diversity > max_diversity:
                    result_queue.put(result)
                    continue

            # Calculate set error against consensus
            if max_error is not None:
                # Delete positions if required and calculate error
                if gap_positions is not None:
                    seq_check = [deleteSeqPositions(s, gap_positions) for s in seq_list]
                else:
                    seq_check = seq_list
                error = calculateSetError(seq_check, consensus)
                result.log['ERROR'] = error

                # If error exceeds threshold, feed result queue and continue
                if error > max_error:
                    result_queue.put(result)
                    continue

            # TODO:  should move this into an improved annotationConsensus function with an action argument
            # Parse copy_field annotations and define consensus annotations
            if copy_fields is not None and copy_actions is not None:
                copy_ann = OrderedDict()
                for f, act in zip(copy_fields, copy_actions):
                    # Numeric operations
                    if act == 'min':
                        vals = getAnnotationValues(seq_list, f, delimiter=delimiter)
                        copy_ann[f] = '%.12g' % min([float(x or 0) for x in vals])
                    elif act == 'max':
                        vals = getAnnotationValues(seq_list, f, delimiter=delimiter)
                        copy_ann[f] = '%.12g' % max([float(x or 0) for x in vals])
                    elif act == 'sum':
                        vals = getAnnotationValues(seq_list, f, delimiter=delimiter)
                        copy_ann[f] = '%.12g' % sum([float(x or 0) for x in vals])
                    elif act == 'set':
                        vals = annotationConsensus(seq_list, f, delimiter=delimiter)
                        copy_ann[f] = vals['set']
                        copy_ann['%s_COUNT' % f] = vals['count']
                    elif act == 'majority':
                        vals = annotationConsensus(seq_list, f, delimiter=delimiter)
                        copy_ann[f] = vals['cons']
                        copy_ann['%s_FREQ' % f] = vals['freq']
            else:
                copy_ann = None

            # Define annotation for output sequence
            cons_ann = OrderedDict([('ID', data.id),
                                    ('CONSCOUNT', cons_count)])

            # Merge addition consensus annotations into output sequence annotations
            if primer_ann is not None:
                cons_ann = mergeAnnotation(cons_ann, primer_ann, delimiter=delimiter)
            if copy_ann is not None:
                cons_ann = mergeAnnotation(cons_ann, copy_ann, delimiter=delimiter)

            # Add output sequence annotations to consensus sequence
            consensus.id = consensus.name = flattenAnnotation(cons_ann, delimiter=delimiter)
            consensus.description = ''
            result.results = consensus
            result.valid = True

            # Feed results to result queue
            result_queue.put(result)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        sys.stderr.write('Error processing sequence set with ID: %s\n' % data.id)
        raise
    
    return None


def buildConsensus(seq_file, barcode_field=default_barcode_field, 
                   min_count=default_min_count, min_freq=default_min_freq,
                   min_qual=default_min_qual, primer_field=None, primer_freq=None,
                   max_gap=None, max_error=None, max_diversity=None,
                   copy_fields=None, copy_actions=None, dependent=False,
                   out_args=default_out_args, nproc=None, queue_size=None):
    """
    Generates consensus sequences

    Arguments: 
    seq_file = the sample sequence file name
    barcode_field = the annotation field containing set IDs
    min_count = threshold number of sequences to define a consensus 
    min_freq = the frequency cutoff to assign a base 
    min_qual = the quality cutoff to assign a base
    primer_field = the annotation field containing primer tags;
                   if None do not annotate with primer tags
    primer_freq = the maximum primer tag frequency that must be meet to build a consensus;
                  if None do not filter by primer frequency
    max_gap = the maximum frequency of (., -) characters allowed before
               deleting a position; if None do not delete positions
    max_error = a threshold defining the maximum allowed error rate to retain a read group;
                if None do not calculate error rate
    max_diversity = a threshold defining the average pairwise error rate required to retain a read group;
                    if None do not calculate diversity
    dependent = if False treat barcode group sequences as independent data
    copy_fields = a list of annotations to copy into consensus sequence annotations;
                  if None no additional annotations will be copied
    copy_actions = the list of actions to take for each copy_fields;
                   one of ['set', 'majority', 'min', 'max', 'sum']
    out_args = common output argument dictionary from parseCommonArgs
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc
                    
    Returns: 
    a list of successful output file names
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'BuildConsensus'
    log['FILE'] = os.path.basename(seq_file)
    log['BARCODE_FIELD'] = barcode_field
    log['MIN_COUNT'] = min_count
    log['MIN_FREQUENCY'] = min_freq
    log['MIN_QUALITY'] = min_qual
    log['MAX_GAP'] = max_gap
    log['PRIMER_FIELD'] = primer_field
    log['PRIMER_FREQUENCY'] = primer_freq
    log['MAX_ERROR'] = max_error
    log['MAX_DIVERSITY'] = max_diversity
    log['DEPENDENT'] = dependent
    log['COPY_FIELDS'] = ','.join(copy_fields) if copy_fields is not None else None
    log['COPY_ACTIONS'] = ','.join(copy_actions) if copy_actions is not None else None
    log['NPROC'] = nproc
    printLog(log)
    
    # Set consensus building function
    in_type = getFileType(seq_file)
    if in_type == 'fastq':
        cons_func = qualityConsensus
        cons_args = {'min_qual': min_qual, 
                     'min_freq': min_freq,
                     'dependent': dependent}
    elif in_type == 'fasta':  
        cons_func = frequencyConsensus
        cons_args = {'min_freq': min_freq}
    else:
        sys.exit('ERROR:  Input file must be FASTA or FASTQ')
    
    # Define feeder function and arguments
    index_args = {'field': barcode_field, 'delimiter': out_args['delimiter']}
    feed_func = feedSeqQueue
    feed_args = {'seq_file': seq_file,
                 'index_func': indexSeqSets, 
                 'index_args': index_args}
    # Define worker function and arguments
    work_func = processBCQueue
    work_args = {'cons_func': cons_func, 
                 'cons_args': cons_args,
                 'min_count': min_count,
                 'primer_field': primer_field,
                 'primer_freq': primer_freq,
                 'max_gap': max_gap,
                 'max_error': max_error,
                 'max_diversity': max_diversity,
                 'copy_fields': copy_fields,
                 'copy_actions': copy_actions,
                 'delimiter': out_args['delimiter']}
    # Define collector function and arguments
    collect_func = collectSeqQueue
    collect_args = {'seq_file': seq_file,
                    'task_label': 'consensus',
                    'out_args': out_args,
                    'index_field': barcode_field}
    
    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func, 
                             feed_args, work_args, collect_args, 
                             nproc, queue_size)
        
    # Print log
    result['log']['END'] = 'BuildConsensus'
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
                 consensus-pass
                     consensus reads.
                 consensus-fail
                     raw reads failing consensus filtering criteria.

             output annotation fields:
                 PRIMER
                     a comma delimited list of unique primer annotations found within the
                     barcode read group.
                 PRCOUNT
                     a comma delimited list of the corresponding counts of unique primer
                     annotations.
                 PRCONS
                     the majority primer within the barcode read group.
                 PRFREQ
                     the frequency of the majority primer.
                 CONSCOUNT
                     the count of reads within the barcode read group which contributed to
                     the consensus sequence. This is the total size of the read group,
                     minus sequence excluded due to user defined filtering criteria.
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[getCommonArgParser(multiproc=True)], 
                            formatter_class=CommonHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    
    parser.add_argument('-n', action='store', dest='min_count', type=int, default=default_min_count,
                        help='The minimum number of sequences needed to define a valid consensus')
    parser.add_argument('--bf', action='store', dest='barcode_field', type=str,
                        default=default_barcode_field, 
                        help='Position of description barcode field to group sequences by')
    parser.add_argument('-q', action='store', dest='min_qual', type=int, default=default_min_qual,
                        help='Consensus quality score cut-off under which an ambiguous character is assigned; \
                              does not apply when quality scores are unavailable')
    parser.add_argument('--freq', action='store', dest='min_freq', type=float, default=default_min_freq,
                        help='Fraction of character occurrences under which an ambiguous character is assigned.')
    parser.add_argument('--maxgap', action='store', dest='max_gap', type=float, default=None,
                        help='If specified, this defines a cut-off for the frequency of allowed \
                              gap values for each position. Positions exceeding the threshold \
                              are deleted from the consensus. If not defined, positions are always \
                              retained.')
    parser.add_argument('--pf', action='store', dest='primer_field', type=str, default=None, 
                        help='Specifies the field name of the primer annotations')
    parser.add_argument('--prcons', action='store', dest='primer_freq', type=float, default=None, 
                        help='Specify to define a minimum primer frequency required to assign a consensus primer, \
                              and filter out sequences with minority primers from the consensus building step')
    parser.add_argument('--cf', nargs='+', action='store', dest='copy_fields', type=str, default=None,
                        help='''Specifies a set of additional annotation fields to copy into
                             the consensus sequence annotations.''')
    parser.add_argument('--act', nargs='+', action='store', dest='copy_actions', default=None,
                        choices=['min', 'max', 'sum', 'set', 'majority'],
                        help='''List of actions to take for each copy field which defines how
                             each annotation will be combined into a single value. The actions
                             "min", "max", "sum" perform the corresponding mathematical
                             operation on numeric annotations. The action "set" combines
                             annotations into a comma delimited list of unique values and
                             adds an annotation named <FIELD>_COUNT specifying the count
                             of each item in the set. The action "majority" assigns the
                             most frequent annotation to the consensus annotation and adds
                             an annotation named <FIELD>_FREQ specifying the frequency
                             of the majority value.''')
    parser.add_argument('--dep', action='store_true', dest='dependent',
                        help='Specify to calculate consensus quality with a non-independence assumption')

    # Mutually exclusive argument group
    arg_group = parser.add_mutually_exclusive_group()
    arg_group.add_argument('--maxdiv', action='store', dest='max_diversity', type=float,
                           default=None,
                           help='''Specify to calculate the nucleotide diversity of each read
                                group (average pairwise error rate) and remove groups exceeding
                                the given diversity threshold. Diversity is calculate for all
                                positions within the read group, ignoring any character filtering
                                imposed by the -q, --freq and --maxgap arguments.
                                Mutually exclusive with --maxerror.''')
    arg_group.add_argument('--maxerror', action='store', dest='max_error', type=float,
                           default=None,
                           help='''Specify to calculate the error rate of each read group
                                (rate of mismatches from consensus) and remove groups exceeding
                                the given error threshold. The error rate is calculated against
                                the final consensus sequence, which may include masked positions
                                due to the -q and --freq arguments and may have deleted
                                positions due to the --maxgap argument.
                                Mutually exclusive with --maxdiv.''')

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
    if args_dict['barcode_field']:  args_dict['barcode_field'] = args_dict['barcode_field'].upper()
    if args_dict['primer_field']:  args_dict['primer_field'] = args_dict['primer_field'].upper()
    
    # Check prcons argument dependencies
    if args.primer_freq and not args.primer_field:
        parser.error('You must define a primer field with --prf to use the --prcons option')

    # Check copy field and action arguments
    if bool(args_dict['copy_fields']) ^ bool(args_dict['copy_actions']) or \
       len((args_dict['copy_fields'] or '')) != len((args_dict['copy_actions'] or '')):
            parser.error('You must specify exactly one copy action (--act) per copy field (--cf)')

    # Call buildConsensus for each sample file    
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        buildConsensus(**args_dict)
