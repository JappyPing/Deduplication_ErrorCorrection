#!/usr/bin/env python3
"""
Removes primers and annotates sequences with primer and barcode identifiers
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
from Bio import pairwise2

# Presto imports
from presto.Defaults import default_delimiter, default_out_args
from presto.Commandline import CommonHelpFormatter, getCommonArgParser, parseCommonArgs
from presto.Annotation import parseAnnotation, flattenAnnotation, mergeAnnotation
from presto.Sequence import compilePrimers, getDNAScoreDict, reverseComplement
from presto.IO import readPrimerFile, printLog
from presto.Multiprocessing import SeqResult, manageProcesses, feedSeqQueue, \
                                   collectSeqQueue

# Defaults
default_gap_penalty = (1, 1)
default_max_error = 0.2
default_max_len = 50
default_start = 0


class PrimerAlignment:
    """
    A class defining a primer alignment result

    Variables:
    seq          = input SeqRecord object
    primer       = string defining primer id
    align_seq    = string defining input sequence alignment
    align_primer = string defining the primer alignment
    start        = start position of the alignment in the input sequence
    end          = end position of the alignment in the input sequence
    gaps         = number of gaps in input sequence alignment
    error        = error rate
    rev_primer   = True if alignment is tail-ended
    valid        = True if alignment is valid

    Methods:
    __len__      = evaluates to the length of the variable align_seq
    __nonzero__  = evaluates to the value of the variable valid
    """
    # Instantiation
    def __init__(self, seq=None):
        self.seq = seq
        self.primer = None
        self.align_seq = None
        self.align_primer = None
        self.start = None
        self.end = None
        self.gaps = 0
        self.error = 1
        self.rev_primer = False
        self.valid = False

    # Set boolean evaluation to valid value
    def __bool__(self):
        return self.valid

    # Set length evaluation to length of alignment
    def __len__(self):
        if self.align_seq is None:
            return 0
        else:
            return len(self.align_seq)


def alignPrimers(seq_record, primers, primers_regex=None, max_error=default_max_error,
                 max_len=default_max_len, rev_primer=False, skip_rc=False,
                 gap_penalty=default_gap_penalty,
                 score_dict=getDNAScoreDict(mask_score=(0, 1), gap_score=(0, 0))):
    """
    Performs pairwise local alignment of a list of short sequences against a long sequence

    Arguments: 
    seq_record = a SeqRecord object to align primers against
    primers = dictionary of {names: short IUPAC ambiguous sequence strings}
    primers_regex = optional dictionary of {names: compiled primer regular expressions}
    max_error = maximum acceptable error rate before aligning reverse complement
    max_len = maximum length of sample sequence to align
    rev_primer = if True align with the tail end of the sequence
    skip_rc = if True do not check reverse complement sequences
    gap_penalty = a tuple of positive (gap open, gap extend) penalties
    score_dict = optional dictionary of alignment scores as {(char1, char2): score}

    Returns:
    A PrimerAlignment object
    """
    # Defined undefined parameters
    if primers_regex is None:  primers_regex = compilePrimers(primers)
    seq_record = seq_record.upper()
    rec_len = len(seq_record)
    max_len = min(rec_len, max_len)

    # Create empty return object
    align = PrimerAlignment(seq_record)
    align.rev_primer = rev_primer
    
    # Define sequences to align and assign orientation tags
    if not skip_rc:
        seq_list = [seq_record, reverseComplement(seq_record)]
        seq_list[0].annotations['seqorient'] = 'F'
        seq_list[1].annotations['seqorient'] = 'RC'
    else:
        seq_list = [seq_record]
        seq_list[0].annotations['seqorient'] = 'F'
    
    # Assign primer orientation tags
    for rec in seq_list:
        rec.annotations['prorient'] = 'F' if not rev_primer else 'RC' 
    
    # Attempt regular expression match first
    for rec in seq_list:
        scan_seq = str(rec.seq)
        scan_seq = scan_seq[:max_len] if not rev_primer else scan_seq[-max_len:]
        for adpt_id, adpt_regex in primers_regex.items():
            adpt_match = adpt_regex.search(scan_seq)
            # Parse matches
            if adpt_match:
                align.seq = rec
                align.seq.annotations['primer'] = adpt_id
                align.primer = adpt_id
                align.align_seq = scan_seq
                align.align_primer = '-' * adpt_match.start(0) + \
                                     primers[adpt_id] + \
                                     '-' * (max_len - adpt_match.end(0))
                align.gaps = 0
                align.error = 0
                align.valid = True

                # Determine start and end positions
                if not rev_primer:
                    align.start = adpt_match.start(0)
                    align.end = adpt_match.end(0)
                else:
                    rev_pos = rec_len - max_len
                    align.start = adpt_match.start(0) + rev_pos
                    align.end = adpt_match.end(0) + rev_pos

                return align
    
    # Perform local alignment if regular expression match fails
    best_align, best_rec, best_adpt, best_error = None, None, None, None
    for rec in seq_list:
        this_align = dict()
        scan_seq = str(rec.seq)
        scan_seq = scan_seq[:max_len] if not rev_primer else scan_seq[-max_len:]
        for adpt_id, adpt_seq in primers.items():
            pw2_align = pairwise2.align.localds(scan_seq, adpt_seq, score_dict,
                                                -gap_penalty[0], -gap_penalty[1],
                                                one_alignment_only=True)
            if pw2_align:
                this_align.update({adpt_id: pw2_align[0]})
        if not this_align:  continue
        
        # Determine alignment with lowest error rate
        for x_adpt, x_align in this_align.items():
            x_error = 1.0 - x_align[2] / len(primers[x_adpt])
            #x_gaps = len(x_align[1]) - max_len
            #x_error = 1.0 - (x_align[2] + x_gaps) / primers[x_adpt])
            if best_error is None or x_error < best_error:
                best_align = this_align
                best_rec = rec
                best_adpt = x_adpt
                best_error = x_error
        
        # Skip rev_primer complement if forward sequence error within defined threshold
        if best_error <= max_error:  break

    # Set return object to lowest error rate alignment
    if best_align:
        # Define input alignment string and gap count
        align_primer = best_align[best_adpt][1]
        align_len = len(align_primer)
        align_gaps = align_len - max_len

        # Populate return object
        align.seq = best_rec
        align.primer = best_adpt
        align.align_seq = str(best_align[best_adpt][0])
        align.align_primer = align_primer
        align.gaps = align_gaps
        align.error = best_error
        align.valid = True

        # Determine start and end positions
        if not rev_primer:
            # TODO:  need to switch to an aligner that outputs start/end for both sequences in alignment
            align.start = align_len - len(align_primer.lstrip('-'))
            align.end = best_align[best_adpt][4] - align_gaps
        else:
            # Count position from tail and end gaps
            rev_pos = rec_len - align_len
            align.start = rev_pos + best_align[best_adpt][3] + align_gaps
            align.end = rev_pos + len(align_primer.rstrip('-'))

    return align


def scorePrimers(seq_record, primers, start=default_start, rev_primer=False, 
                 score_dict=getDNAScoreDict(mask_score=(0, 1), gap_score=(0, 0))):
    """
    Performs simple alignment of primers with a fixed starting position, 
    no reverse complement alignment, and no tail alignment option

    Arguments: 
    seq_record = a SeqRecord object to align primers against
    primers = dictionary of {names: short IUPAC ambiguous sequence strings}
    start = position where primer alignment starts
    rev_primer = if True align with the tail end of the sequence
    score_dict = optional dictionary of {(char1, char2): score} alignment scores
    
    Returns:
    A PrimerAlignment object
    """
    # Create empty return dictionary
    seq_record = seq_record.upper()
    align = PrimerAlignment(seq_record)
    align.rev_primer = rev_primer

    # Define orientation variables
    seq_record.annotations['seqorient'] = 'F'
    seq_record.annotations['prorient'] = 'F' if not rev_primer else 'RC'

    # Score primers
    this_align = {}
    rec_len = len(seq_record)
    if rev_primer:  end = rec_len - start
    for adpt_id, adpt_seq in primers.items():
        if rev_primer:  start = end - len(adpt_seq)
        else:  end = start + len(adpt_seq)
        chars = zip(seq_record[start:end], adpt_seq)
        score = sum([score_dict[(c1, c2)] for c1, c2 in chars])
        this_align.update({adpt_id: (score, start, end)})

    # Determine primer with lowest error rate
    best_align, best_adpt, best_err = None, None, None
    for adpt, algn in this_align.items():
        #adpt_err = 1.0 - float(algn[0]) / weightSeq(primers[adpt])
        err = 1.0 - float(algn[0]) / len(primers[adpt])
        if best_err is None or err < best_err:
            best_align = algn
            best_adpt = adpt
            best_err = err

    # Set return dictionary to lowest error rate alignment
    if best_align:
        # Populate return object
        align.primer = best_adpt if best_err < 1.0 else None
        align.start = best_align[1]
        align.end = best_align[2]
        align.error = best_err
        align.valid = True

        # Determine alignment sequences
        if not rev_primer:
            align.align_seq = str(seq_record.seq[:best_align[2]])
            align.align_primer = '-' * best_align[1] + primers[best_adpt]
        else:
            align.align_seq = str(seq_record.seq[best_align[1]:])
            align.align_primer = primers[best_adpt] + '-' * (rec_len - best_align[2])
    
    return align


def getMaskedSeq(align, mode='mask', barcode=False, delimiter=default_delimiter):
    """
    Create an output sequence with primers masked or cut

    Arguments: 
    align = a PrimerAlignment object returned from alignPrimers or scorePrimers
    mode = defines the action taken; one of ['cut','mask','tag','trim']
    barcode = if True add sequence preceding primer to description
    delimiter = a tuple of delimiters for (annotations, field/values, value lists) 

    Returns:
    output SeqRecord object
    """
    seq = align.seq

    # Build output sequence
    if mode == 'tag' or not align.align_primer:
        # Do not modify sequence
        out_seq = seq
    elif mode == 'trim':
        # Remove region before primer
        if not align.rev_primer:
            out_seq = seq[align.start:]
        else:  
            out_seq = seq[:align.end]
    elif mode == 'cut':
        # Remove primer and preceding region
        if not align.rev_primer:
            out_seq = seq[align.end:]
        else: 
            out_seq = seq[:align.start]
    elif mode == 'mask':
        # Mask primer with Ns and remove preceding region
        if not align.rev_primer:
            mask_len = align.end - align.start + align.gaps
            out_seq = 'N' * mask_len + seq[align.end:]
            if hasattr(seq, 'letter_annotations') and \
                    'phred_quality' in seq.letter_annotations:
                out_seq.letter_annotations['phred_quality'] = \
                    [0] * mask_len + \
                    seq.letter_annotations['phred_quality'][align.end:]
        else:
            mask_len = min(align.end, len(seq)) - align.start + align.gaps
            out_seq = seq[:align.start] + 'N' * mask_len
            if hasattr(seq, 'letter_annotations') and \
                    'phred_quality' in seq.letter_annotations:
                out_seq.letter_annotations['phred_quality'] = \
                    seq.letter_annotations['phred_quality'][:align.start] + \
                    [0] * mask_len
            
    # Add alignment annotations to output SeqRecord
    out_seq.annotations = seq.annotations    
    out_seq.annotations['primer'] = align.primer
    out_seq.annotations['prstart'] = align.start
    out_seq.annotations['error'] = align.error

    # Parse seq annotation and create output annotation
    seq_ann = parseAnnotation(seq.description, delimiter=delimiter)
    out_ann = OrderedDict([('SEQORIENT', seq.annotations['seqorient']),
                           ('PRIMER', align.primer)])
    
    # Add ID sequence to description
    if barcode:
        seq_code = seq[:align.start].seq if not align.rev_primer \
                   else seq[align.end:].seq
        out_seq.annotations['barcode'] = seq_code
        out_ann['BARCODE'] = seq_code
    
    out_ann = mergeAnnotation(seq_ann, out_ann, delimiter=delimiter)
    out_seq.id = flattenAnnotation(out_ann, delimiter=delimiter)
    out_seq.description = ''

    return out_seq


def processMPQueue(alive, data_queue, result_queue, align_func, align_args={}, 
                   mask_args={}, max_error=default_max_error):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments: 
    alive = a multiprocessing.Value boolean controlling whether processing 
            continues; when False function returns
    data_queue = a multiprocessing.Queue holding data to process
    result_queue = a multiprocessing.Queue to hold processed results
    align_func = the function to call for alignment
    align_args = a dictionary of arguments to pass to align_func
    mask_args = a dictionary of arguments to pass to getMaskedSeq
    max_error = maximum acceptable error rate for a valid alignment
    
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
            
            # Define result object for iteration
            in_seq = data.data
            result = SeqResult(in_seq.id, in_seq)
     
            # Align primers
            align = align_func(in_seq, **align_args)
            
            # Process alignment results
            if not align:
                # Update log if no alignment
                result.log['ALIGN'] = None
            else:
                # Create output sequence
                out_seq = getMaskedSeq(align, **mask_args)        
                result.results = out_seq
                result.valid = bool(align.error <= max_error) if len(out_seq) > 0 else False
                
                # Update log with successful alignment results
                result.log['SEQORIENT'] = out_seq.annotations['seqorient']
                result.log['PRIMER'] = out_seq.annotations['primer']
                result.log['PRORIENT'] = out_seq.annotations['prorient']
                result.log['PRSTART'] = out_seq.annotations['prstart']
                if 'barcode' in out_seq.annotations:  
                    result.log['BARCODE'] = out_seq.annotations['barcode']
                if not align.rev_primer:
                    align_cut = len(align.align_seq) - align.gaps
                    result.log['INSEQ'] = align.align_seq + \
                                          str(align.seq.seq[align_cut:])
                    result.log['ALIGN'] = align.align_primer
                    result.log['OUTSEQ'] = str(out_seq.seq).rjust(len(in_seq) + align.gaps)
                else:
                    align_cut = len(align.seq) - len(align.align_seq) + align.gaps
                    result.log['INSEQ'] = str(align.seq.seq[:align_cut]) + align.align_seq
                    result.log['ALIGN'] = align.align_primer.rjust(len(in_seq) + align.gaps)
                    result.log['OUTSEQ'] = str(out_seq.seq)
                result.log['ERROR'] = align.error
            
            # Feed results to result queue
            result_queue.put(result)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        sys.stderr.write('Error processing sequence with ID: %s.\n' % data.id)
        raise
    
    return None


def maskPrimers(seq_file, primer_file, mode, align_func, align_args={}, 
                max_error=default_max_error, barcode=False,
                out_args=default_out_args, nproc=None, queue_size=None):
    """
    Masks or cuts primers from sample sequences using local alignment

    Arguments: 
    seq_file = name of file containing sample sequences
    primer_file = name of the file containing primer sequences
    mode = defines the action taken; one of 'cut','mask','tag'
    align_func = the function to call for alignment
    align_arcs = a dictionary of arguments to pass to align_func
    max_error = maximum acceptable error rate for a valid alignment
    barcode = if True add sequence preceding primer to description
    out_args = common output argument dictionary from parseCommonArgs
    nproc = the number of processQueue processes;
            if None defaults to the number of CPUs
    queue_size = maximum size of the argument queue;
                 if None defaults to 2*nproc
                 
    Returns:
    a list of successful output file names
    """
    # Define subcommand label dictionary
    cmd_dict = {alignPrimers:'align', scorePrimers:'score'}
    
    # Print parameter info
    log = OrderedDict()
    log['START'] = 'MaskPrimers'
    log['COMMAND'] = cmd_dict.get(align_func, align_func.__name__)
    log['SEQ_FILE'] = os.path.basename(seq_file)
    log['PRIMER_FILE'] = os.path.basename(primer_file)
    log['MODE'] = mode
    log['BARCODE'] = barcode
    log['MAX_ERROR'] = max_error
    if 'start' in align_args: log['START_POS'] = align_args['start']
    if 'max_len' in align_args: log['MAX_LEN'] = align_args['max_len']
    if 'rev_primer' in align_args: log['REV_PRIMER'] = align_args['rev_primer']
    if 'skip_rc' in align_args: log['SKIP_RC'] = align_args['skip_rc']
    if 'gap_penalty' in align_args:
        log['GAP_PENALTY'] = ', '.join([str(x) for x in align_args['gap_penalty']])
    log['NPROC'] = nproc
    printLog(log)

    # Create dictionary of primer sequences to pass to maskPrimers
    primers = readPrimerFile(primer_file)
    if 'rev_primer' in align_args and align_args['rev_primer']:
        primers = {k: reverseComplement(v) for k, v in primers.items()}

    # Define alignment arguments and compile primers for align mode
    align_args['primers'] = primers 
    align_args['score_dict'] = getDNAScoreDict(mask_score=(0, 1), gap_score=(0, 0))
    if align_func is alignPrimers:
        align_args['max_error'] = max_error
        align_args['primers_regex'] = compilePrimers(primers)
    
    # Define sequence masking arguments
    mask_args = {'mode': mode, 
                 'barcode': barcode, 
                 'delimiter': out_args['delimiter']}

    # Define feeder function and arguments
    feed_func = feedSeqQueue
    feed_args = {'seq_file': seq_file}
    # Define worker function and arguments
    work_func = processMPQueue
    work_args = {'align_func': align_func, 
                 'align_args': align_args,
                 'mask_args': mask_args,
                 'max_error': max_error}
    
    # Define collector function and arguments
    collect_func = collectSeqQueue
    collect_args = {'seq_file': seq_file,
                    'task_label': 'primers',
                    'out_args': out_args}
    
    # Call process manager
    result = manageProcesses(feed_func, work_func, collect_func, 
                             feed_args, work_args, collect_args, 
                             nproc, queue_size)

    # Print log
    result['log']['END'] = 'MaskPrimers'
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
                 mask-pass
                     processed reads with successful primer matches.
                 mask-fail
                     raw reads failing primer identification.

             output annotation fields:
                 SEQORIENT
                     the orientation of the output sequence. Either F (input) or RC
                     (reverse complement of input).
                 PRIMER
                     name of the best primer match.
                 BARCODE
                     the sequence preceding the primer match. Only output when the
                     --barcode flag is specified.
             ''')

    # Define ArgumentParser
    parser = ArgumentParser(description=__doc__, epilog=fields,
                            formatter_class=CommonHelpFormatter)
    parser.add_argument('--version', action='version',
                        version='%(prog)s:' + ' %s-%s' %(__version__, __date__))
    subparsers = parser.add_subparsers(title='subcommands', metavar='',
                                       help='Alignment method')
    # TODO:  This is a temporary fix for Python issue 9253
    subparsers.required = True
    
    # Parent parser
    parser_parent = getCommonArgParser(multiproc=True)
    parser_parent.add_argument('-p', action='store', dest='primer_file', required=True, 
                               help='A FASTA or REGEX file containing primer sequences.')
    parser_parent.add_argument('--mode', action='store', dest='mode',
                               choices=('cut', 'mask', 'trim', 'tag'), default='mask',
                               help='''Specifies the action to take with the primer sequence.
                                    The "cut" mode will remove both the primer region and
                                    the preceding sequence. The "mask" mode will replace the
                                    primer region with Ns and remove the preceding sequence.
                                    The "trim" mode will remove the region preceding the primer,
                                    but leave the primer region intact. The "tag" mode will
                                    leave the input sequence unmodified.''')
    parser_parent.add_argument('--maxerror', action='store', dest='max_error', type=float,
                               default=default_max_error, help='Maximum allowable error rate.')
    parser_parent.add_argument('--revpr', action='store_true', dest='rev_primer', 
                              help='Specify to match the tail-end of the sequence against the \
                                    reverse complement of the primers.')
    parser_parent.add_argument('--barcode', action='store_true', dest='barcode', 
                               help='''Specify to encode sequences with barcode sequences
                                    (unique molecular identifiers) found preceding the primer
                                    region.''')
    
    # Align mode argument parser
    parser_align = subparsers.add_parser('align', parents=[parser_parent],
                                         formatter_class=CommonHelpFormatter,
                                         help='Find primer matches using pairwise local alignment')
    parser_align.add_argument('--maxlen', action='store', dest='max_len', type=int,
                              default=default_max_len, help='Maximum sequence length to scan for primers.')
    parser_align.add_argument('--skiprc', action='store_true', dest='skip_rc', 
                              help='Specify to prevent checking of sample reverse complement sequences.')
    parser_align.add_argument('--gap', nargs=2, action='store', dest='gap_penalty',
                              type=float, default=default_gap_penalty,
                              help='''A list of two positive values defining the gap open
                                   and gap extension penalties for aligning the primers.
                                   Note: the error rate is calculated as the percentage
                                   of mismatches from the primer sequence with gap
                                   penalties reducing the match count accordingly; this may
                                   lead to error rates that differ from strict mismatch
                                   percentage when gaps are present in the alignment.''')
    #parser_align.set_defaults(start=None)
    parser_align.set_defaults(align_func=alignPrimers)
    

    # Score mode argument parser
    parser_score = subparsers.add_parser('score', parents=[parser_parent], 
                                         formatter_class=CommonHelpFormatter,
                                         help='Find primer matches by scoring primers at a fixed position')
    parser_score.add_argument('--start', action='store', dest='start', type=int, default=default_start, 
                              help='The starting position of the primer')

    parser_score.set_defaults(align_func=scorePrimers)
    
    return parser


if __name__ == '__main__':
    """
    Parses command line arguments and calls main function
    """
    # Parse arguments
    parser = getArgParser()
    args = parser.parse_args()
    args_dict = parseCommonArgs(args)
    
    # Define align_args dictionary to pass to maskPrimers
    if args_dict['align_func'] is alignPrimers:
        args_dict['align_args'] = {'max_len':args_dict['max_len'],
                                   'rev_primer':args_dict['rev_primer'],
                                   'skip_rc':args_dict['skip_rc'],
                                   'gap_penalty':args_dict['gap_penalty']}
        del args_dict['max_len']
        del args_dict['rev_primer']
        del args_dict['skip_rc']
        del args_dict['gap_penalty']
    elif args_dict['align_func'] is scorePrimers:
        args_dict['align_args'] = {'start':args_dict['start'],
                                   'rev_primer':args_dict['rev_primer']}
        del args_dict['start']
        del args_dict['rev_primer']
    
    # Call maskPrimers for each sample file
    del args_dict['seq_files']
    for f in args.__dict__['seq_files']:
        args_dict['seq_file'] = f
        maskPrimers(**args_dict)
    