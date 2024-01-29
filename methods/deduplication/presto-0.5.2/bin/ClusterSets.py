#!/usr/bin/env python3
"""
Cluster sequences by group
"""
# Info
__author__ = 'Christopher Bolen, Jason Anthony Vander Heiden'
from presto import __version__, __date__

# Imports
import os
import sys
import tempfile
from argparse import ArgumentParser
from collections import OrderedDict
from textwrap import dedent

# Presto imports
from presto.Defaults import default_delimiter, default_barcode_field, \
                            default_cluster_field, default_out_args, \
                            default_usearch_exec
from presto.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from presto.Annotation import parseAnnotation, flattenAnnotation, mergeAnnotation
from presto.Applications import runUClust
from presto.IO import printLog
from presto.Sequence import indexSeqSets
from presto.Multiprocessing import SeqResult, manageProcesses, feedSeqQueue, \
                                   collectSeqQueue

# Defaults
default_ident = 0.9


def processCSQueue(alive, data_queue, result_queue, cluster_field,
                  cluster_args={}, delimiter=default_delimiter):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments:
    alive = a multiprocessing.Value boolean controlling whether processing
            continues; when False function returns
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    cluster_args = a dictionary of optional arguments for the clustering function
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

            # Perform clustering
            cluster_dict = runUClust(data.data, **cluster_args)

            # Process failed result
            if cluster_dict is None:
                # Update log
                result.log['CLUSTERS'] = 0
                for i, seq in enumerate(data.data, start=1):
                    result.log['CLUST0-%i' % i] = str(seq.seq)

                # Feed results queue and continue
                result_queue.put(result)
                continue

            # Get number of clusters
            result.log['CLUSTERS'] = len(cluster_dict)

            # Update sequence annotations with cluster assignments
            results = list()
            seq_dict = {s.id: s for s in data.data}
            for clust, id_list in cluster_dict.items():
                for i, seq_id in enumerate(id_list, start=1):
                    # Add cluster annotation
                    seq = seq_dict[seq_id]
                    header = parseAnnotation(seq.description, delimiter=delimiter)
                    header = mergeAnnotation(header, {cluster_field:clust}, delimiter=delimiter)
                    seq.id = seq.name = flattenAnnotation(header, delimiter=delimiter)
                    seq.description = ''

                    # Update log and results
                    result.log['CLUST%i-%i' % (clust, i)] = str(seq.seq)
                    results.append(seq)

            # Check results
            result.results = results
            result.valid = (len(results) == len(seq_dict))

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


def clusterSets(seq_file, barcode_field=default_barcode_field,
                cluster_field=default_cluster_field,
                ident=default_ident, seq_start=None, seq_end=None,
                usearch_exec=default_usearch_exec,
                out_args=default_out_args, nproc=None,
                queue_size=None):
    """
    Performs clustering on sets of sequences

    Arguments:
    seq_file = the sample sequence file name
    barcode_field = the annotation containing set IDs
    ident = the identity threshold for clustering sequences
    seq_start = the start position to trim sequences at before clustering
    seq_end = the end position to trim sequences at before clustering
    usearch_exec = the path to the executable for usearch
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc

    Returns:
    the clustered output file name
    """
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'ClusterSets'
    log['FILE'] = os.path.basename(seq_file)
    log['BARCODE_FIELD'] = barcode_field
    log['CLUSTER_FIELD'] = cluster_field
    log['IDENTITY'] = ident
    log['SEQUENCE_START'] = seq_start
    log['SEQUENCE_END'] = seq_end
    log['NPROC'] = nproc
    printLog(log)

    # Define cluster function parameters
    cluster_args = {'usearch_exec':usearch_exec,
                    'ident':ident,
                    'seq_start':seq_start,
                    'seq_end':seq_end}

    # Define feeder function and arguments
    index_args = {'field': barcode_field, 'delimiter': out_args['delimiter']}
    feed_func = feedSeqQueue
    feed_args = {'seq_file': seq_file,
                 'index_func': indexSeqSets,
                 'index_args': index_args}
    # Define worker function and arguments
    work_func = processCSQueue
    work_args = {'cluster_field': cluster_field,
                 'cluster_args': cluster_args,
                 'delimiter': out_args['delimiter']}
    # Define collector function and arguments
    collect_func = collectSeqQueue
    collect_args = {'seq_file': seq_file,
                    'task_label': 'cluster',
                    'out_args': out_args,
                    'index_field': barcode_field}

    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func,
                             feed_args, work_args, collect_args,
                             nproc, queue_size)

    # Print log
    log = OrderedDict()
    log['OUTPUT'] = result['log'].pop('OUTPUT')
    for k, v in result['log'].items():  log[k] = v
    log['END'] = 'ClusterSets'
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
                 cluster-pass
                    clustered reads.
                 cluster-fail
                    raw reads failing clustering.

             output annotation fields:
                 CLUSTER
                    a numeric cluster identifier defining the within-group cluster.
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[getCommonArgParser(multiproc=True)],
                            formatter_class=CommonHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s-%s' %(__version__, __date__))

    parser.add_argument('-f', action='store', dest='barcode_field', type=str,
                        default=default_barcode_field,
                        help='''The annotation field containing annotations, such as UID
                             barcode, for sequence grouping.''')
    parser.add_argument('-k', action='store', dest='cluster_field', type=str,
                        default=default_cluster_field,
                        help='''The name of the output annotation field to add with the
                             cluster information for each sequence.''')
    parser.add_argument('--id', action='store', dest='ident', type=float,
                        default=default_ident,
                        help='The sequence identity threshold for the usearch algorithm.')
    parser.add_argument('--start', action='store', dest='seq_start', type=int,
                        help='''The start of the region to be used for clustering.
                             Together with --end, this parameter can be used to specify a
                             subsequence of each read to use in the clustering algorithm.''')
    parser.add_argument('--end', action='store', dest='seq_end', type=int,
                        help='The end of the region to be used for clustering.')
    parser.add_argument('--exec', action='store', dest='usearch_exec',
                        default=default_usearch_exec,
                        help='The location of the USEARCH executable.')

    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    # Parse arguments
    parser = getArgParser()
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)

    # Convert fields to uppercase
    if 'barcode_field' in args_dict and args_dict['barcode_field'] is not None:
        args_dict['barcode_field'] = args_dict['barcode_field'].upper()
    if 'cluster_field' in args_dict and args_dict['cluster_field'] is not None:
        args_dict['cluster_field'] = args_dict['cluster_field'].upper()
    
    # Check if a valid usearch executable was specified
    if not os.path.isfile(args.usearch_exec):
        parser.error('%s does not exist' % args.usearch_exec)
    
        
    # Call cluster for each input file
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        clusterSets(**args_dict)
