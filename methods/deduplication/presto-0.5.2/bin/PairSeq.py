#!/usr/bin/env python3
"""
Sorts and matches sequence records with matching coordinates across files
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from presto import __version__, __date__

# Imports
import os
from argparse import ArgumentParser
from collections import OrderedDict
from textwrap import dedent
from time import time
from Bio import SeqIO

# Presto imports
from presto.Defaults import default_coord_choices, default_coord_type, default_out_args
from presto.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from presto.Annotation import parseAnnotation, flattenAnnotation, mergeAnnotation, \
                              getCoordKey
from presto.IO import getFileType, readSeqFile, countSeqFile, getOutputHandle, \
                      printLog, printMessage, printProgress


def pairSeq(seq_file_1, seq_file_2, fields_1=None, fields_2=None,
            coord_type=default_coord_type,
            out_args=default_out_args):
    """
    Generates consensus sequences

    Arguments: 
    seq_file_1 = the file containing the grouped sequences and annotations
    seq_file_2 = the file to assign annotations to from seq_file_1
    fields_1 = list of annotations in seq_file_1 records to copy to seq_file_2 records;
               if None do not copy any annotations
    fields_2 = list of annotations in seq_file_2 records to copy to seq_file_1 records;
               if None do not copy any annotations
    coord_type = the sequence header format
    out_args = common output argument dictionary from parseCommonArgs
                    
    Returns: 
    a list of tuples holding successfully paired filenames for (seq_file_1, seq_file_2)
    """
    # Define private functions
    def _key_func(x):
        return getCoordKey(x, coord_type=coord_type, delimiter=out_args['delimiter'])

    log = OrderedDict()
    log['START'] = 'PairSeq'
    log['FILE1'] = os.path.basename(seq_file_1)
    log['FILE2'] = os.path.basename(seq_file_2)
    log['FIELDS_1'] = ','.join(fields_1) if fields_1 is not None else None
    log['FIELDS_2'] = ','.join(fields_2) if fields_2 is not None else None
    log['COORD_TYPE'] = coord_type
    printLog(log)

    # Define output type
    if out_args['out_type'] is None:
        out_type_1 = getFileType(seq_file_1)
        out_type_2 = getFileType(seq_file_2)
    else: 
        out_type_1 = out_type_2 = out_args['out_type']

    # Define output name
    if out_args['out_name'] is None:
        out_name_1 = out_name_2 = None
    else: 
        out_name_1 = '%s-1' % out_args['out_name']
        out_name_2 = '%s-2' % out_args['out_name']

    # Open and count files
    start_time = time()
    printMessage("Indexing files", start_time=start_time)
    # Index file 1
    seq_count_1 = countSeqFile(seq_file_1)
    seq_dict_1 = readSeqFile(seq_file_1, index=True, key_func=_key_func)
    # Define file 2 iterator
    seq_count_2 = countSeqFile(seq_file_2)
    seq_iter_2 = readSeqFile(seq_file_2, index=False)
    printMessage("Done", start_time=start_time, end=True)

    # Open output file handles
    pass_handle_1 = getOutputHandle(seq_file_1, 'pair-pass', out_args['out_dir'], 
                                    out_name=out_name_1, out_type=out_type_1)
    pass_handle_2 = getOutputHandle(seq_file_2, 'pair-pass', out_args['out_dir'], 
                                    out_name=out_name_2, out_type=out_type_2)

    if out_args['failed']:
        fail_handle_1 = getOutputHandle(seq_file_1, 'pair-fail', out_dir=out_args['out_dir'],
                                        out_name=out_name_1, out_type=out_type_1)
        fail_handle_2 = getOutputHandle(seq_file_2, 'pair-fail', out_dir=out_args['out_dir'],
                                        out_name=out_name_2, out_type=out_type_2)
        pass_keys = list()

    # Iterate over pairs and write to output files
    start_time = time()
    rec_count = pair_count = 0
    for seq_2 in seq_iter_2:
        # Print progress for previous iteration
        printProgress(rec_count, seq_count_2, 0.05, start_time)
        rec_count += 1

        # Check for file 2 mate pair in file 1
        coord_2 = getCoordKey(seq_2.id, coord_type=coord_type,
                              delimiter=out_args['delimiter'])
        seq_1 = seq_dict_1.get(coord_2, None)

        if seq_1 is not None:
            # Record paired keys
            pair_count += 1

            if fields_1 is not None or fields_2 is not None:
                ann_1 = parseAnnotation(seq_1.description,
                                        delimiter=out_args['delimiter'])
                ann_2 = parseAnnotation(seq_2.description,
                                        delimiter=out_args['delimiter'])

                # Prepend annotations from seq_1 to seq_2
                if fields_1 is not None:
                    copy_ann = OrderedDict([(k, v) for k, v in ann_1.items() \
                                            if k in fields_1])
                    merge_ann = mergeAnnotation(ann_2, copy_ann, prepend=True,
                                                delimiter=out_args['delimiter'])
                    seq_2.id = flattenAnnotation(merge_ann,
                                                 delimiter=out_args['delimiter'])
                    seq_2.description = ''

                # Append annotations from seq_2 to seq_1
                if fields_2 is not None:
                    copy_ann = OrderedDict([(k, v) for k, v in ann_2.items() \
                                            if k in fields_2])
                    merge_ann = mergeAnnotation(ann_1, copy_ann, prepend=False,
                                                delimiter=out_args['delimiter'])
                    seq_1.id = flattenAnnotation(merge_ann,
                                                 delimiter=out_args['delimiter'])
                    seq_1.description = ''

            # Write paired records
            SeqIO.write(seq_1, pass_handle_1, out_type_1)
            SeqIO.write(seq_2, pass_handle_2, out_type_2)

        # Write unpaired file 2 records and updated paired key list for finding unpaired file 1 records
        if out_args['failed']:
            if seq_1 is not None:  pass_keys.append(coord_2)
            else:  SeqIO.write(seq_2, fail_handle_2, out_type_2)

    # Print final progress
    printProgress(rec_count, seq_count_2, 0.05, start_time)

    # Find and write unpaired file 1 records
    if out_args['failed']:
        start_time = time()
        printMessage("Finding unpaired", start_time=start_time)

        # Find file 1 unpaired keys
        pass_keys = set(pass_keys)
        unpaired = set(seq_dict_1).difference(pass_keys)
        # Write unpaired file 1 records
        for k in unpaired:  SeqIO.write(seq_dict_1[k], fail_handle_1, out_type_1)

        printMessage("Done", start_time=start_time, end=True)

    # Print log
    log = OrderedDict()
    log['OUTPUT1'] = os.path.basename(pass_handle_1.name)
    log['OUTPUT2'] = os.path.basename(pass_handle_2.name)
    log['SEQUENCES1'] = seq_count_1
    log['SEQUENCES2'] = seq_count_2
    log['PASS'] = pair_count
    log['END'] = 'PairSeq'
    printLog(log)
   
    # Close file handles
    pass_handle_1.close()
    pass_handle_2.close()

    return [(pass_handle_1.name, pass_handle_2.name)]


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
                 pair-pass
                     successfully paired reads with modified annotations.
                 pair-fail
                     raw reads that could not be assigned to a mate-pair.

             output annotation fields:
                 <user defined>
                     annotation fields specified by the --1f or --2f arguments.
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[getCommonArgParser(paired=True, failed=True, log=False)],
                            formatter_class=CommonHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s-%s' %(__version__, __date__))

    parser.add_argument('--1f', nargs='+', action='store', dest='fields_1', type=str, default=None,
                        help='''The annotation fields to copy from file 1 records into file
                             2 records. If a copied annotation already exists in a file 2
                             record, then the annotations copied from file 1 will be added
                             to the front of the existing annotation.''')
    parser.add_argument('--2f', nargs='+', action='store', dest='fields_2', type=str, default=None,
                        help='''The annotation fields to copy from file 2 records into file
                             1 records. If a copied annotation already exists in a file 1
                             record, then the annotations copied from file 2 will be added
                             to the end of the existing annotation.''')
    parser.add_argument('--coord', action='store', dest='coord_type', 
                        choices=default_coord_choices, default=default_coord_type,
                        help='''The format of the sequence identifier which defines shared
                             coordinate information across mate pairs.''')
    
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
    if args_dict['fields_1'] is not None:
        args_dict['fields_1'] = list(map(str.upper, args_dict['fields_1']))
    if args_dict['fields_2'] is not None:
        args_dict['fields_2'] = list(map(str.upper, args_dict['fields_2']))

    # Call pairSeq
    del args_dict['seq_files_1']
    del args_dict['seq_files_2']
    for file_1, file_2 in zip(args.__dict__['seq_files_1'], 
                              args.__dict__['seq_files_2']):
        args_dict['seq_file_1'] = file_1
        args_dict['seq_file_2'] = file_2
        pairSeq(**args_dict)
