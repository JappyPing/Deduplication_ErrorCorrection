#!/usr/bin/env python3
"""
Parses pRESTO annotations in FASTA/FASTQ sequence headers
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from presto import __version__, __date__

# Imports
import csv
import os
from argparse import ArgumentParser
from collections import OrderedDict

from textwrap import dedent
from time import time
from Bio import SeqIO

# Presto imports
from presto.Defaults import default_delimiter, default_separator, default_out_args
from presto.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from presto.Annotation import parseAnnotation, flattenAnnotation, mergeAnnotation, \
                              renameAnnotation, collapseAnnotation
from presto.IO import getFileType, readSeqFile, countSeqFile, getOutputHandle, \
                      printLog, printProgress


def addHeader(header, fields, values, delimiter=default_delimiter):
    """
    Adds fields and values to a sequence header

    Arguments: 
    header = an annotation dictionary returned by parseAnnotation
    fields = the list of fields to add or append to
    values = the list of annotation values to add for each field
    delimiter = a tuple of delimiters for (fields, values, value lists)
                    
    Returns: 
    the modified header dictionary
    """
    for f, v in zip(fields, values):
        header = mergeAnnotation(header, {f:v}, delimiter=delimiter)
        
    return header


def collapseHeader(header, fields, actions, delimiter=default_delimiter):
    """
    Collapses a sequence header

    Arguments: 
    header = an annotation dictionary returned by parseAnnotation
    fields = the list of fields to collapse
    actions = the list of collapse action take;
              one of (max, min, sum, first, last, set, cat) for each field
    delimiter = a tuple of delimiters for (fields, values, value lists)
                    
    Returns: 
    the output file name
    """
    for f, a in zip(fields, actions):
        header = collapseAnnotation(header, a, f, delimiter=delimiter)
        
    return header


def copyHeader(header, fields, names, actions=None, delimiter=default_delimiter):
    """
    Copies fields in a sequence header

    Arguments:
    header = an annotation dictionary returned by parseAnnotation
    fields = a list of the field names to copy
    names = a list of the new field names
    actions = the list of collapse action take after the copy;
              one of (max, min, sum, first, last, set, cat) for each field
    delimiter = a tuple of delimiters for (fields, values, value lists)

    Returns:
    the modified header dictionary
    """
    old_header = header.copy()
    for f, n in zip(fields, names):
        header = mergeAnnotation(header, {n: old_header[f]}, delimiter=delimiter)

    if actions is not None:
        header = collapseHeader(header, names, actions, delimiter=delimiter)

    return header


def deleteHeader(header, fields, delimiter=default_delimiter):
    """
    Deletes fields from a sequence header

    Arguments: 
    header = an annotation dictionary returned by parseAnnotation
    fields = the list of fields to delete
    delimiter = a tuple of delimiters for (fields, values, value lists)
                        
    Returns: 
    the modified header dictionary
    """
    for f in fields:  del header[f]

    return header


def expandHeader(header, fields, separator=default_separator, 
                 delimiter=default_delimiter):
    """
    Splits and annotation value into separate fields in a sequence header

    Arguments: 
    header = an annotation dictionary returned by parseAnnotation
    fields = the field to split
    separator = the delimiter to split the values by
    delimiter = a tuple of delimiters for (fields, values, value lists)
                        
    Returns: 
    the modified header dictionary
    """
    for f in fields:
        values = header[f].split(separator)
        names = [f + str(i + 1) for i in range(len(values))]
        ann = OrderedDict([(n, v) for n, v in zip(names, values)])
        header = mergeAnnotation(header, ann, delimiter=delimiter)
        del header[f]
    
    return header


def renameHeader(header, fields, names, actions=None, delimiter=default_delimiter):
    """
    Renames fields in a sequence header

    Arguments: 
    header = an annotation dictionary returned by parseAnnotation
    fields = a list of the current field names
    names = a list of the new field names
    actions = the list of collapse action take after the rename;
              one of (max, min, sum, first, last, set, cat) for each field
    delimiter = a tuple of delimiters for (fields, values, value lists)
                            
    Returns: 
    the modified header dictionary
    """
    for f, n in zip(fields, names):
        header = renameAnnotation(header, f, n, delimiter=delimiter)

    if actions is not None:
        header = collapseHeader(header, names, actions, delimiter=delimiter)

    return header


def modifyHeaders(seq_file, modify_func, modify_args, out_args=default_out_args):
    """
    Modifies sequence headers

    Arguments: 
    seq_file = the sequence file name
    modify_func = the function defining the modification operation
    modify_args = a dictionary of arguments to pass to modify_func
    out_args = common output argument dictionary from parseCommonArgs
                    
    Returns: 
    the output file name
    """
    # Define subcommand label dictionary
    cmd_dict = {addHeader: 'add',
                copyHeader: 'copy',
                collapseHeader: 'collapse',
                deleteHeader: 'delete',
                expandHeader: 'expand',
                renameHeader: 'rename'}
    
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'ParseHeaders'
    log['COMMAND'] = cmd_dict.get(modify_func, modify_func.__name__)
    log['FILE'] = os.path.basename(seq_file)
    for k in sorted(modify_args):  
        v = modify_args[k]
        log[k.upper()] = ','.join(v) if isinstance(v, list) else v
    printLog(log)
    
    # Open file handles
    in_type = getFileType(seq_file)
    seq_iter = readSeqFile(seq_file)
    if out_args['out_type'] is None:  out_args['out_type'] = in_type
    out_handle = getOutputHandle(seq_file, 'reheader', out_dir=out_args['out_dir'],
                                 out_name=out_args['out_name'], out_type=out_args['out_type'])

    # Count records
    result_count = countSeqFile(seq_file)
    
    # Iterate over sequences
    start_time = time()
    seq_count = 0
    for seq in seq_iter:
        # Print progress for previous iteration
        printProgress(seq_count, result_count, 0.05, start_time)
        
        #Update counts
        seq_count += 1
        
        # Modify header
        header = parseAnnotation(seq.description, delimiter=out_args['delimiter'])
        header = modify_func(header, delimiter=out_args['delimiter'], **modify_args)
        
        # Write new sequence
        seq.id = seq.name = flattenAnnotation(header, delimiter=out_args['delimiter'])
        seq.description = ''
        SeqIO.write(seq, out_handle, out_args['out_type'])
        
    # print counts
    printProgress(seq_count, result_count, 0.05, start_time)    
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(out_handle.name)
    log['SEQUENCES'] = seq_count
    log['END'] = 'ParseHeaders'               
    printLog(log)

    # Close file handles
    out_handle.close()
 
    return out_handle.name


def tableHeaders(seq_file, fields, out_args=default_out_args):
    """
    Builds a table of sequence header annotations

    Arguments: 
    seq_file = the sequence file name
    fields = the list of fields to output
    out_args = common output argument dictionary from parseCommonArgs
                    
    Returns: 
    the output table file name
    """
    log = OrderedDict()
    log['START'] = 'ParseHeaders'
    log['COMMAND'] = 'table'
    log['FILE'] = os.path.basename(seq_file)
    printLog(log)
    
    # Open file handles
    seq_iter = readSeqFile(seq_file)
    out_handle = getOutputHandle(seq_file, out_label='headers', out_dir=out_args['out_dir'], 
                                 out_name=out_args['out_name'], out_type='tab')
    # Count records
    result_count = countSeqFile(seq_file)
    
    # Open csv writer and write header
    out_writer = csv.DictWriter(out_handle, extrasaction='ignore', restval='', 
                                delimiter='\t', fieldnames=fields)
    out_writer.writeheader()
    
    # Iterate over sequences
    start_time = time()
    seq_count = pass_count = fail_count = 0
    for seq in seq_iter:
        # Print progress for previous iteration
        printProgress(seq_count, result_count, 0.05, start_time)
        
        # Get annotations
        seq_count += 1
        ann = parseAnnotation(seq.description, fields, delimiter=out_args['delimiter'])

        # Write records
        if ann:
            pass_count += 1
            out_writer.writerow(ann)
        else:
            fail_count += 1
        
    # Print counts
    printProgress(seq_count, result_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(out_handle.name)
    log['SEQUENCES'] = seq_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['END'] = 'ParseHeaders'
    printLog(log)

    # Close file handles
    out_handle.close()
 
    return out_handle.name


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
                 reheader-pass
                     reads passing annotation operation and modified accordingly.
                 reheader-fail
                     raw reads failing annotation operation.
                 headers
                     tab delimited table of the selected annotations.

             output annotation fields:
                 <user defined>
                     annotation fields specified by the -f argument.
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            formatter_class=CommonHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    subparsers = parser.add_subparsers(title='subcommands', dest='command', metavar='',
                                       help='Annotation operation')
    # TODO:  This is a temporary fix for Python issue 9253
    subparsers.required = True

    # Subparser to add header fields
    parser_add = subparsers.add_parser('add', parents=[getCommonArgParser(log=False)],
                                       formatter_class=CommonHelpFormatter,
                                       help='Adds field/value pairs to header annotations')    
    parser_add.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                            help='List of fields to add.')
    parser_add.add_argument('-u', nargs='+', action='store', dest='values', required=True,
                            help='List of values to add for each field.')
    parser_add.set_defaults(func=modifyHeaders)
    parser_add.set_defaults(modify_func=addHeader)

    # Subparser to collapse header fields
    parser_collapse = subparsers.add_parser('collapse', parents=[getCommonArgParser(log=False)],
                                            formatter_class=CommonHelpFormatter,
                                            help='Collapses header annotations with multiple entries')    
    parser_collapse.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                                 help='List of fields to collapse.')
    parser_collapse.add_argument('--act', nargs='+', action='store', dest='actions', required=True,
                                 choices=['min', 'max', 'sum', 'first', 'last', 'set', 'cat'],
                                 help='''List of actions to take for each field defining how
                                      each annotation will be combined into a single value.
                                      The actions "min", "max", "sum" perform the corresponding
                                      mathematical operation on numeric annotations. The
                                      actions "first" and "last" choose the value from the
                                      corresponding position in the annotation. The action
                                      "set" collapses annotations into a comma delimited
                                      list of unique values. The action "cat" concatenates
                                      the values together into a single string.''')
    parser_collapse.set_defaults(func=modifyHeaders)
    parser_collapse.set_defaults(modify_func=collapseHeader)

    # Subparser to copy header fields
    parser_copy = subparsers.add_parser('copy', parents=[getCommonArgParser(log=False)],
                                        formatter_class=CommonHelpFormatter,
                                        help='Copies header annotation fields')
    parser_copy.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                               help='List of fields to copy.')
    parser_copy.add_argument('-k', nargs='+', action='store', dest='names', required=True,
                               help='''List of names for each copied field. If the new field
                                    is already present, the copied field will be merged into
                                    the existing field.''')
    parser_copy.add_argument('--act', nargs='+', action='store', dest='actions', required=False,
                                 choices=['min', 'max', 'sum', 'first', 'last', 'set', 'cat'],
                                 help='''List of collapse actions to take on each new field
                                      following the copy operation defining how each annotation
                                      will be combined into a single value. The actions
                                      "min", "max", "sum" perform the corresponding
                                      mathematical operation on numeric annotations. The
                                      actions "first" and "last" choose the value from the
                                      corresponding position in the annotation. The action
                                      "set" collapses annotations into a comma delimited
                                      list of unique values. The action "cat" concatenates
                                      the values together into a single string.''')
    parser_copy.set_defaults(func=modifyHeaders)
    parser_copy.set_defaults(modify_func=copyHeader)

    # Subparser to delete header fields
    parser_delete = subparsers.add_parser('delete', parents=[getCommonArgParser(log=False)],
                                          formatter_class=CommonHelpFormatter,
                                          help='Deletes fields from header annotations')
    parser_delete.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                               help='List of fields to delete.')
    parser_delete.set_defaults(func=modifyHeaders)
    parser_delete.set_defaults(modify_func=deleteHeader)

    # Subparser to expand header fields
    parser_expand = subparsers.add_parser('expand', parents=[getCommonArgParser(log=False)],
                                          formatter_class=CommonHelpFormatter,
                                          help='Expands annotation fields with multiple values')    
    parser_expand.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                               help='List of fields to expand.')
    parser_expand.add_argument('--sep', action='store', dest='separator', 
                               default=default_separator,
                               help='The character separating each value in the fields.')
    parser_expand.set_defaults(func=modifyHeaders)
    parser_expand.set_defaults(modify_func=expandHeader)

    # Subparser to rename header fields
    parser_rename = subparsers.add_parser('rename', parents=[getCommonArgParser(log=False)],
                                          formatter_class=CommonHelpFormatter,
                                          help='Renames header annotation fields')
    parser_rename.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                               help='List of fields to rename.')
    parser_rename.add_argument('-k', nargs='+', action='store', dest='names', required=True,
                               help='''List of new names for each field. If the new field is
                                    already present, the renamed field will be merged into
                                    the existing field and the old field will be deleted.''')
    parser_rename.add_argument('--act', nargs='+', action='store', dest='actions', required=False,
                               choices=['min', 'max', 'sum', 'first', 'last', 'set', 'cat'],
                               help='''List of collapse actions to take on each new field
                                    following the rename operation defining how each annotation
                                    will be combined into a single value. The actions
                                    "min", "max", "sum" perform the corresponding
                                    mathematical operation on numeric annotations. The
                                    actions "first" and "last" choose the value from the
                                    corresponding position in the annotation. The action
                                    "set" collapses annotations into a comma delimited
                                    list of unique values. The action "cat" concatenates
                                    the values together into a single string.''')
    parser_rename.set_defaults(func=modifyHeaders)
    parser_rename.set_defaults(modify_func=renameHeader)
            
    # Subparser to create a header table
    parser_table = subparsers.add_parser('table', parents=[getCommonArgParser(seq_out=False, log=False)],
                                         formatter_class=CommonHelpFormatter,
                                         help='Writes sequence headers to a table')
    parser_table.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                              help='''List of fields to collect. The sequence identifier may
                                   be specified using the hidden field name "ID".''')
    parser_table.set_defaults(func=tableHeaders)

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
    if 'fields' in args_dict and args_dict['fields']:  
        args_dict['fields'] = list(map(str.upper, args_dict['fields']))
    # Built modify_args dictionary
    if args.func == modifyHeaders:
        modify_args = {}
        if 'fields' in args_dict:  
            modify_args['fields'] = args_dict.pop('fields')
        if 'actions' in args_dict:
            modify_args['actions'] = args_dict.pop('actions') if args_dict['actions'] is None \
                                     else list(map(str.lower, args_dict.pop('actions')))
        if 'names' in args_dict:  
            modify_args['names'] = list(map(str.upper, args_dict.pop('names')))
        if 'values' in args_dict: 
            modify_args['values'] = args_dict.pop('values')
        if 'separator' in args_dict: 
            modify_args['separator'] = args_dict.pop('separator')
        args_dict['modify_args'] = modify_args
        
    # Check modify_args arguments
    if args.command == 'add' and len(modify_args['fields']) != len(modify_args['values']):
        parser.error('You must specify exactly one value (-u) per field (-f)')
    if args.command == 'collapse' and len(modify_args['fields']) != len(modify_args['actions']):
        parser.error('You must specify exactly one action (-a) per field (-f)')
    if args.command in ('copy', 'rename'):
        if len(modify_args['fields']) != len(modify_args['names']):
            parser.error('You must specify exactly one new name (-k) per field (-f)')
        if modify_args['actions'] is not None and len(modify_args['actions']) != len(modify_args['names']):
            parser.error('You must specify exactly one action (--act) per new name (-k)')

    # Calls header processing function
    del args_dict['command']
    del args_dict['func']
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        args.func(**args_dict)
