#!/usr/bin/env python3
"""
Converts sequence headers to the pRESTO format
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from presto import __version__, __date__

# Imports
import os
import re
from argparse import ArgumentParser
from collections import OrderedDict
from textwrap import dedent
from time import time
from Bio import SeqIO

# Presto imports
from presto.Defaults import default_delimiter, default_out_args
from presto.Annotation import parseAnnotation, flattenAnnotation
from presto.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from presto.IO import getFileType, readSeqFile, countSeqFile, getOutputHandle, \
                      printLog, printProgress


def convertGenericHeader(desc, delimiter=default_delimiter):
    """
    Converts any header to the pRESTO format

    Arguments:
    desc = a sequence description string
    delimiter = a tuple of delimiters for (fields, values, value lists)

    Returns:
    a dictionary of header {field: value} pairs
    """
    # Replace whitespace and delimiter characters
    sub_regex = '[%s\s]+' % re.escape(''.join(delimiter))
    conv = re.sub(sub_regex, '_', desc)
    try:
        # Check if modified header is valid
        header = parseAnnotation(conv, delimiter=delimiter)
    except:
        # Assign header to None if header cannot be converted
        header = None

    return header


def convert454Header(desc):
    """
    Parses 454 headers into the pRESTO format

    Arguments:
    desc = a sequence description string

    Returns:
    a dictionary of header {field: value} pairs

    New header example:
        @GXGJ56Z01AE06X length=222
        @<accession> <length=##>

    Old header example:
        @000034_0199_0169 length=437 uaccno=GNDG01201ARRCR
        @<rank_x_y> <length=##> <uaccno=accession>
    """
    # Split description and assign field names
    try:
        # Build header dictionary
        fields = desc.split(' ')
        header = OrderedDict()
        header['ID'] = fields[0]
        header['LENGTH'] = fields[1].replace('length=', '')

        # Check for old format
        if len(fields) == 3:
            header['UACCNO'] = fields[2].replace('uaccno=', '')
        elif len(fields) != 2:
            raise
    except:
        header = None

    return header


def convertGenbankHeader(desc, delimiter=default_delimiter):
    """
    Converts Genbank and RefSeq headers into the pRESTO format

    Arguments:
    desc = a sequence description string
    delimiter = a tuple of delimiters for (fields, values, value lists)

    Returns:
    a dictionary of header {field: value} pairs

    New header example:
        >CM000663.2 Homo sapiens chromosome 1, GRCh38 reference primary assembly
        <accession>.<version> <description>
    Old header example:
        >gi|568336023|gb|CM000663.2| Homo sapiens chromosome 1, GRCh38 reference primary assembly
        gi|<GI record number>|<dbsrc>|<accession>.<version>|<description>
    """
    # Define special characters to replace
    sub_regex = '[%s\s]+' % re.escape(''.join(delimiter[1:]))

    # Split description and assign field names
    try:
        header = OrderedDict()

        # Try old format and fallback to new format if that fails
        fields = desc.split('|')
        if len(fields) == 5:
            header['ID'] = fields[3]
            header['GI'] = fields[1]
            header['SOURCE'] = fields[2]
            header['DESC'] = re.sub(sub_regex, '_', fields[4].strip())
        else:
            fields = desc.split(' ')
            header['ID'] = fields[0]
            header['DESC'] = re.sub(sub_regex, '_', '_'.join(fields[1:]).strip())
    except:
        header = None

    return header


def convertIlluminaHeader(desc):
    """
    Converts Illumina headers into the pRESTO format

    Arguments:
    desc = a sequence description string

    Returns:
    a dictionary of header {field: value} pairs

    New header example:
        @MISEQ:132:000000000-A2F3U:1:1101:14340:1555 2:N:0:ATCACG
        @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read number>:<is filtered>:<control number>:<index sequence>
    Old header example:
        @HWI-EAS209_0006_FC706VJ:5:58:5894:21141#ATCACG/1
        @<instrument>:<flowcell lane>:<tile>:<x-pos>:<y-pos>#<index sequence>/<read number>
    """
    # Split description and assign field names
    try:
        # Try new format and fallback to old if that fails
        fields = desc.split(' ')
        if len(fields) == 2:
            x = fields[1].split(':')
            index = x[3]
            read_num = x[0]
        else:
            fields = desc.split('#')
            x = fields[1].split('/')
            index = x[0]
            read_num = x[1]

        # Build header dictionary
        header = OrderedDict()
        header['ID'] = fields[0]
        header['INDEX'] = index
        header['READ'] = read_num
    except:
        header = None

    return header


def convertIMGTHeader(desc, simple=False):
    """
    Converts germline headers from IMGT/GENE-DB into the pRESTO format

    Arguments:
    desc = a sequence description string
    simple = if True then the header will be converted to only the allele name

    Returns:
    a dictionary of header {field: value} pairs

    Header specifications from http://imgt.org/genedb
        The FASTA header contains 15 fields separated by '|':
         1. IMGT/LIGM-DB accession number(s)
         2. gene and allele name
         3. species
         4. functionality
         5. exon(s), region name(s), or extracted label(s)
         6. start and end positions in the IMGT/LIGM-DB accession number(s)
         7. number of nucleotides in the IMGT/LIGM-DB accession number(s)
         8. codon start, or 'NR' (not relevant) for non coding labels and
            out-of-frame pseudogenes
         9. +n: number of nucleotides (nt) added in 5' compared to the
            corresponding label extracted from IMGT/LIGM-DB
        10. +n or -n: number of nucleotides (nt) added or removed in 3'
            compared to the corresponding label extracted from IMGT/LIGM-DB
        11. +n, -n, and/or nS: number of added, deleted, and/or substituted
            nucleotides to correct sequencing errors, or 'not corrected' if
            non corrected sequencing errors
        12. number of amino acids (AA): this field indicates that the
            sequence is in amino acids
        13. number of characters in the sequence: nt (or AA)+IMGT gaps=total
        14. partial (if it is)
        15. reverse complementary (if it is)

    Header example:
        >X60503|IGHV1-18*02|Homo sapiens|F|V-REGION|142..417|276 nt|1| | | | |276+24=300|partial in 3'| |
    """
    # Split description and assign field names
    try:
        fields = desc.split('|')

        # Build header dictionary
        header = OrderedDict()
        header['ID'] = fields[1]

        if not simple:
            header['SPECIES'] = re.sub('\s', '_', fields[2])
            header['REGION'] = fields[4]
            header['FUNCTIONALITY'] = re.sub('[\(\)\[\]]', '', fields[3])
            header['PARTIAL'] = 'FALSE' if re.sub('\s', '', fields[13]) == '' else 'TRUE'
            header['ACCESSION'] = fields[0]

        # Position and length data
        #header['NUCLEOTIDES'] = re.sub('[^0-9]', '', fields[6])
        #header['LENGTH'] = fields[12].split('=')[1]
    except:
        header = None

    return header


def convertSRAHeader(desc):
    """
    Parses NCBI SRA headers into the pRESTO format

    Arguments:
    desc = a sequence description string

    Returns:
    a dictionary of header {field: value} pairs

    Header example from fastq-dum --split-files:
        @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
        @SRR1383326.1 1 length=250
        @<accession>.<spot> <original sequence description> <length=#>
    Header example from fastq-dum --split-files -I:
        @SRR1383326.1.1 1 length=250
        @<accession>.<spot>.<read number> <original sequence description> <length=#>
    """
    # Split description and assign field names
    try:
        fields = desc.split(' ')


        # Build header dictionary
        header = OrderedDict()

        # Check for read number if sequence id
        read_id = fields[0].split('.')
        if len(read_id) == 3:
            header['ID'] = '.'.join(read_id[:2])
            header['READ'] = read_id[2]
        else:
            header['ID'] = fields[0]

        header['DESC'] = fields[1]
        header['LENGTH'] = fields[2].replace('length=', '')
    except:
        header = None

    return header


def convertHeaders(seq_file, convert_func, convert_args={}, out_args=default_out_args):
    """
    Converts sequence headers to the pRESTO format

    Arguments:
    seq_file = the sequence file name
    convert_func = the function used to convert sequence headers
    convert_args = a dictionary of arguments to pass to convert_func
    out_args = common output argument dictionary from parseCommonArgs

    Returns:
    the output sequence file name
    """
    # Define subcommand label dictionary
    cmd_dict = {convertGenericHeader:'generic',
                convert454Header:'454',
                convertGenbankHeader:'genbank',
                convertIlluminaHeader:'illumina',
                convertIMGTHeader:'imgt',
                convertSRAHeader:'sra'}

    log = OrderedDict()
    log['START'] = 'ConvertHeaders'
    log['COMMAND'] = cmd_dict[convert_func]
    log['FILE'] = os.path.basename(seq_file)
    printLog(log)

    # Open input file
    in_type = getFileType(seq_file)
    seq_iter = readSeqFile(seq_file)
    if out_args['out_type'] is None:  out_args['out_type'] = in_type

    # Count records
    result_count = countSeqFile(seq_file)

    # Open output file handles
    pass_handle = getOutputHandle(seq_file,
                                  'convert-pass',
                                  out_dir=out_args['out_dir'],
                                  out_name=out_args['out_name'],
                                  out_type=out_args['out_type'])
    if out_args['failed']:
        fail_handle = getOutputHandle(seq_file,
                                      'convert-fail',
                                      out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'],
                                      out_type=out_args['out_type'])
    else:
        fail_handle = None

    # Set additional conversion arguments
    if convert_func in [convertGenericHeader, convertGenbankHeader]:
        convert_args.update({'delimiter':out_args['delimiter']})

    # Iterate over sequences
    start_time = time()
    seq_count = pass_count = fail_count = 0
    for seq in seq_iter:
        # Print progress for previous iteration and update count
        printProgress(seq_count, result_count, 0.05, start_time)
        seq_count += 1

        # Convert header
        header = convert_func(seq.description, **convert_args)

        if header is not None:
            # Write successfully converted sequences
            pass_count += 1
            seq.id = seq.name = flattenAnnotation(header, out_args['delimiter'])
            seq.description = ''
            SeqIO.write(seq, pass_handle, out_args['out_type'])
        else:
            fail_count += 1
            if fail_handle is not None:
                # Write successfully unconverted sequences
                SeqIO.write(seq, fail_handle, out_args['out_type'])

    # Print counts
    printProgress(seq_count, result_count, 0.05, start_time)
    log = OrderedDict()
    log['OUTPUT'] = os.path.basename(pass_handle.name)
    log['SEQUENCES'] = seq_count
    log['PASS'] = pass_count
    log['FAIL'] = fail_count
    log['END'] = 'ConvertHeaders'
    printLog(log)

    # Close file handles
    pass_handle.close()
    if fail_handle is not None:  fail_handle.close()

    return pass_handle.name


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
                 convert-pass
                     reads passing header conversion.
                 convert-fail
                     raw reads failing header conversion.

             output annotation fields:
                 <format defined>
                     the annotation fields added are specific to the header format of the
                     input file.
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            formatter_class=CommonHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    subparsers = parser.add_subparsers(title='subcommands', metavar='',
                                       help='Conversion method')
    # TODO:  This is a temporary fix for Python issue 9253
    subparsers.required = True

    # Parent parser defining universal argument
    parser_parent = getCommonArgParser(log=False)

    # Subparser for generic header conversion
    parser_generic = subparsers.add_parser('generic', parents=[parser_parent],
                                       formatter_class=CommonHelpFormatter,
                                       help='''Converts sequence headers without a known
                                            annotation system.''')
    parser_generic.set_defaults(convert_func=convertGenericHeader)

    # Subparser for conversion of 454 headers
    parser_454 = subparsers.add_parser('454', parents=[parser_parent],
                                       formatter_class=CommonHelpFormatter,
                                       help='''Converts Roche 454 sequence headers.''')
    parser_454.set_defaults(convert_func=convert454Header)

    # Subparser for conversion of GenBank and RefSeq headers
    parser_genbank = subparsers.add_parser('genbank', parents=[parser_parent],
                                           formatter_class=CommonHelpFormatter,
                                           help='''Converts NCBI GenBank and RefSeq
                                                sequence headers.''')
    parser_genbank.set_defaults(convert_func=convertGenbankHeader)

    # Subparser for conversion of Illumina headers
    parser_illumina = subparsers.add_parser('illumina', parents=[parser_parent],
                                            formatter_class=CommonHelpFormatter,
                                            help='''Converts Illumina sequence headers.''')
    parser_illumina.set_defaults(convert_func=convertIlluminaHeader)

    # Subparser for conversion of IMGT germline headers
    parser_imgt = subparsers.add_parser('imgt', parents=[parser_parent],
                                        formatter_class=CommonHelpFormatter,
                                        help='''Converts sequence headers output by
                                             IMGT/GENE-DB.''')
    parser_imgt.add_argument('--simple', action='store_true', dest='simple',
                             help='''If specified, only the allele name, and no other
                                  annotations, will appear in the converted sequence
                                  header.''')
    parser_imgt.set_defaults(convert_func=convertIMGTHeader)

    # Subparser for conversion of SRA headers
    parser_sra = subparsers.add_parser('sra', parents=[parser_parent],
                                       formatter_class=CommonHelpFormatter,
                                       help='''Converts NCBI SRA sequence headers.''')
    parser_sra.set_defaults(convert_func=convertSRAHeader)

    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    # Parse arguments
    parser = getArgParser()
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)

    # Create convert_args
    convert_keys = ['simple']
    args_dict['convert_args'] = dict((k, args_dict[k]) for k in args_dict \
                                     if k in convert_keys)
    for k in args_dict['convert_args']:  del args_dict[k]

    # Calls header conversion function
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        convertHeaders(**args_dict)
