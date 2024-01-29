"""
External application wrappers
"""
# Info
__author__    = 'Jason Anthony Vander Heiden, Namita Gupta'
from presto import __version__, __date__

# Imports
import csv
import tempfile
import pandas as pd
from io import StringIO
from subprocess import CalledProcessError, check_output, PIPE, Popen, STDOUT
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

# Presto imports
from presto.Defaults import default_muscle_exec, default_usearch_exec, default_blastn_exec

# Defaults
default_ident = 0.9
default_evalue = 1e-5
default_max_hits = 100


def runMuscle(seq_list, muscle_exec=default_muscle_exec):
    """
    Multiple aligns a set of sequences using MUSCLE

    Arguments:
    seq_list = a list of SeqRecord objects to align
    muscle_exec = the MUSCLE executable

    Returns:
    a MultipleSeqAlignment object containing the alignment
    """
    # Return sequence if only one sequence in seq_list
    if len(seq_list) < 2:
        align = MultipleSeqAlignment(seq_list)
        return align

    # Set MUSCLE command
    cmd = [muscle_exec, '-diags', '-maxiters', '2']

    # Convert sequences to FASTA and write to string
    stdin_handle = StringIO()
    SeqIO.write(seq_list, stdin_handle, 'fasta')
    stdin_str = stdin_handle.getvalue()
    stdin_handle.close()

    # Open MUSCLE process
    child = Popen(cmd, bufsize=-1, stdin=PIPE, stdout=PIPE, stderr=PIPE,
                  universal_newlines=True)

    # Send sequences to MUSCLE stdin and retrieve stdout, stderr
    stdout_str, __ = child.communicate(stdin_str)

    # Capture sequences from MUSCLE stdout
    stdout_handle = StringIO(stdout_str)
    align = AlignIO.read(stdout_handle, 'fasta')
    stdout_handle.close()

    return align


def runUClust(seq_list, ident=default_ident, seq_start=0, seq_end=None,
              usearch_exec=default_usearch_exec):
    """
    Cluster a set of sequences using the UCLUST algorithm from USEARCH

    Arguments:
    seq_list = a list of SeqRecord objects to align
    ident = the sequence identity cutoff to be passed to usearch
    seq_start = the start position to trim sequences at before clustering
    seq_end = the end position to trim sequences at before clustering
    usearch_exec = the path to the usearch executable

    Returns:
    a dictionary object containing {sequence id: cluster id}
    """
    # Return sequence if only one sequence in seq_list
    if len(seq_list) < 2:
        #return {seq_list[0].id:0}
        return {1:[seq_list[0].id]}

    # Format sequences and make a copy so we don't mess up original sequences
    short_list = list()
    for rec in seq_list:
        seq = rec.seq[seq_start:seq_end]
        seq = seq.ungap('-')
        seq = seq.ungap('.')
        short_list.append(SeqRecord(seq, id=rec.id, name=rec.name,
                                    description=rec.description))

    # Open temporary files
    in_handle = tempfile.NamedTemporaryFile(mode='w+t')
    out_handle = tempfile.NamedTemporaryFile(mode='w+t')

    # Define usearch command
    cmd = [usearch_exec,
           '-cluster_fast', in_handle.name,
           '-uc', out_handle.name,
           '-id', str(ident),
           '-minseqlength', '1',
           '-threads', '1']

    # Write usearch input fasta file
    SeqIO.write(short_list, in_handle, 'fasta')
    in_handle.seek(0)

    # Run usearch uclust algorithm
    try:
        stdout_str = check_output(cmd, stderr=STDOUT, shell=False,
                                  universal_newlines=True)
        #check_call(cmd, stderr=STDOUT, shell=False)
    except CalledProcessError:
        group_dict = None
    else:
        # TODO:  unsure about this return object.
        # Parse the results of usearch
        # Output columns for the usearch 'uc' output format
        #   0 = entry type -- S: centroid seq, H: hit, C: cluster record (redundant with S)
        #   1 = group the sequence is assigned to
        #   8 = the id of the sequence
        group_dict = {}
        for row in csv.reader(out_handle, delimiter='\t'):
            if row[0] in ('S', 'H'):
                key = int(row[1]) + 1
                group = group_dict.setdefault(key, [])
                group.append(row[8])
        #out_list = [r for r in csv.reader(out_handle, delimiter='\t')]
        #group_dict = {r[8]: int(r[1]) + 1 for r in out_list if r[0] in ('S', 'H')}

    return group_dict if group_dict else None


def makeUSearchDb(ref_file, usearch_exec=default_usearch_exec):
    """
    Makes a usearch database file for ublast

    Arguments:
    ref_file = the path to the reference database file
    usearch_exec = the path to the usearch executable

    Returns:
    a handle to the named temporary file containing the database file
    """
    # Open temporary file
    db_handle = tempfile.NamedTemporaryFile(suffix='.udb')

    # Define usearch command
    cmd = ' '.join([usearch_exec,
               '-makeudb_ublast', ref_file,
               '-output', db_handle.name])

    child = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=(sys.platform != 'win32'))
    stdout_str, stderr_str = child.communicate()

    return db_handle


def runUBlastAlignment(seq, ref_file, evalue=default_evalue, max_hits=default_max_hits,
                       usearch_exec=default_usearch_exec):
    """
    Aligns a sequence against a reference database using the UBLAST algorithm of USEARCH

    Arguments:
    seq = a SeqRecord object to align
    ref_file = the path to the reference database file
    evalue = the E-value cut-off for ublast
    max_hits = the maxhits output limit for ublast
    usearch_exec = the path to the usearch executable

    Returns:
    a DataFrame of alignment results
    """
    # Open temporary files
    in_handle = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')
    out_handle = tempfile.NamedTemporaryFile(mode='w+t', encoding='utf-8')

    # Define usearch command
    cmd = [usearch_exec,
           '-ublast', in_handle.name,
           '-db', ref_file,
           '-strand', 'plus',
           '-evalue', str(evalue),
           '-maxhits', str(max_hits),
           '-userout', out_handle.name,
           '-userfields', 'query+target+qlo+qhi+tlo+thi+alnlen+evalue+id',
           '-threads', '1']


    # Write usearch input fasta file
    SeqIO.write(seq, in_handle, 'fasta')

    # Run usearch ublast algorithm
    in_handle.seek(0)
    stdout_str = check_output(cmd, stderr=STDOUT, shell=False, universal_newlines=True)

    # Parse usearch output
    field_names = ['query', 'target', 'query_start', 'query_end',
                   'target_start', 'target_end',
                   'length', 'evalue', 'identity']
    align_df = pd.read_table(out_handle, header=None, names=field_names, encoding='utf-8')
    # Convert to base-zero indices
    align_df[['query_start', 'query_end', 'target_start', 'target_end']] -= 1

    # Close temp file handles
    in_handle.close()
    out_handle.close()

    return align_df



def runBlastnAlignment(seq, ref_file, evalue=default_evalue, max_hits=default_max_hits,
                       blastn_exec=default_blastn_exec):
    """
    Aligns a sequence against a reference database using BLASTN

    Arguments:
    seq = a SeqRecord objects to align
    ref_dict = a dictionary of reference SeqRecord objects
    evalue = the E-value cut-off for ublast
    maxhits = the maxhits output limit for ublast
    blastn_exec = the path to the usearch executable

    Returns:
    a DataFrame of alignment results
    """
    seq_fasta = seq.format('fasta')

    # Define blastn command
    cmd = ' '.join([blastn_exec,
                    '-query -',
                    #'-query', str(seq.seq),
                    '-subject', ref_file,
                    '-strand plus',
                    '-evalue', str(evalue),
                    '-max_target_seqs', str(max_hits),
                    '-outfmt "6 qseqid sseqid qstart qend sstart send length evalue pident"',
                    '-num_threads 1'])

    # Run blastn
    child = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE,
                  shell=(sys.platform != 'win32'), universal_newlines=True)
    stdout_str, stderr_str = child.communicate(seq_fasta)
    out_handle = StringIO(stdout_str)

    # Parse blastn output
    field_names = ['query', 'target', 'query_start', 'query_end', 'target_start', 'target_end',
                   'length', 'evalue', 'identity']
    align_df = pd.read_table(out_handle, header=None, names=field_names)
    # Convert to base-zero indices
    align_df[['query_start', 'query_end', 'target_start', 'target_end']] -= 1

    return align_df