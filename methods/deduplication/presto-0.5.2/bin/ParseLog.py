#!/usr/bin/env python3
"""
Parses records in the console log of pRESTO modules
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

# Presto imports
from presto.Defaults import default_out_args
from presto.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from presto.IO import getOutputHandle, printLog, printProgress


def parseLogRecord(record):
    """
    Parses an pRESTO log record

    Arguments: 
    record = a string of lines representing a log record
                    
    Returns: 
    an OrderedDict of {field names: values}
    """
    record_dict = OrderedDict()
    for line in record.split('\n'):
        line = line.strip().split('> ')
        if len(line) == 2:
            record_dict[line[0]] = line[1]
    
    return record_dict


def tableLog(record_file, fields, out_args=default_out_args):
    """
    Converts a pRESTO log to a table of annotations

    Arguments: 
    record_file = the log file name
    fields = the list of fields to output
    out_args = common output argument dictionary from parseCommonArgs
                    
    Returns: 
    the output table file name
    """
    log = OrderedDict()
    log['START'] = 'ParseLog'
    log['FILE'] = os.path.basename(record_file)
    printLog(log)
    
    # Open file handles
    log_handle = open(record_file, 'rU')
    out_handle = getOutputHandle(record_file, 
                                 'table', 
                                 out_dir=out_args['out_dir'], 
                                 out_name=out_args['out_name'], 
                                 out_type='tab')
        
    # Open csv writer and write header
    out_writer = csv.DictWriter(out_handle, extrasaction='ignore', restval='', 
                                delimiter='\t', fieldnames=fields)
    out_writer.writeheader()
    
    # Iterate over log records
    start_time = time()
    record = ''
    rec_count = pass_count = fail_count = 0
    for line in log_handle:
        if line.strip() == '' and record:
            # Print progress for previous iteration
            printProgress(rec_count, None, 1e5, start_time)
            
            # Parse record block
            rec_count += 1
            record_dict = parseLogRecord(record)

            # Write records
            if any([f in fields for f in record_dict]):
                pass_count += 1
                out_writer.writerow(record_dict)
            elif record_dict:
                fail_count += 1
                
            # Empty record string
            record = ''
        else:
            # Append to record
            record += line
    else:
        # Write final record
        if record: 
            record_dict = parseLogRecord(record)
            if any([f in fields for f in record_dict]):
                pass_count += 1
                out_writer.writerow(record_dict)
            elif record_dict:
                fail_count += 1
    
    # Print counts
    printProgress(rec_count, None, 1e5, start_time, end=True)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(out_handle.name)
    log['RECORDS'] = rec_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['END'] = 'ParseLog'
    printLog(log)

    # Close file handles
    log_handle.close()
    out_handle.close()
 
    return log_handle.name


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
                 table
                     tab delimited table of the selected annotations.

             output annotation fields:
                 <user defined>
                     annotation fields specified by the -f argument.
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            parents=[getCommonArgParser(seq_in=False, seq_out=False,
                                                        failed=False, log=False)],
                            formatter_class=CommonHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s-%s' %(__version__, __date__))

    parser.add_argument('-l', nargs='+', action='store', dest='record_files', required=True,
                        help='List of log files to parse.')
    parser.add_argument('-f', nargs='+', action='store', dest='fields', required=True,
                        help='''List of fields to collect. The sequence identifier may
                             be specified using the hidden field name "ID".''')
    
    return parser

    
if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    # Parse arguments
    parser = getArgParser()
    args = parser.parse_args()
    args_dict = parseCommonArgs(args, in_arg='record_files')
    # Convert case of fields
    if args_dict['fields']:  args_dict['fields'] = list(map(str.upper, args_dict['fields'])) 
    
    # Call parseLog for each log file
    del args_dict['record_files']
    for f in args.__dict__['record_files']:
        args_dict['record_file'] = f
        tableLog(**args_dict)
