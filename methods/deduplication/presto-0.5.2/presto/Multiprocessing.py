"""
Multiprocessing functions
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from presto import __version__, __date__

# Imports
import ctypes
import os
import signal
import sys
import multiprocessing as mp
from collections import OrderedDict
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from presto.IO import getFileType, readSeqFile, countSeqFile, countSeqSets, \
                      getOutputHandle, printLog, printProgress

# Constants
TERMINATION_SENTINEL = None
EXCEPTION_SENTINEL = None


class SeqData:
    """
    A class defining sequence data objects for worker processes
    """
    # Instantiation
    def __init__(self, key, records):
        self.id = key
        self.data = records
        self.valid = (key is not None and records is not None)

    # Set boolean evaluation to valid value
    def __bool__(self):
        return self.valid

    # Set length evaluation to number of data records
    def __len__(self):
        if isinstance(self.data, SeqRecord) or isinstance(self.data, Seq):
            return 1
        elif self.data is None:
            return 0
        else:
            return len(self.data)


class SeqResult:
    """
    A class defining sequence result objects for collector processes
    """
    # Instantiation
    def __init__(self, key, records):
        self.id = key
        self.data = records
        self.results = None
        self.valid = False
        self.log = OrderedDict([('ID', key)])

    # Set boolean evaluation to valid value
    def __bool__(self):
        return self.valid

    # Set length evaluation to number of results
    def __len__(self):
        if isinstance(self.results, SeqRecord) or isinstance(self.results, Seq):
            return 1
        elif self.results is None:
            return 0
        else:
            return len(self.results)

    # Set data_count to number of data records
    @property
    def data_count(self):
        if isinstance(self.data, SeqRecord) or isinstance(self.data, Seq):
            return 1
        elif self.data is None:
            return 0
        else:
            return len(self.data)


def manageProcesses(feed_func, work_func, collect_func,
                    feed_args={}, work_args={}, collect_args={},
                    nproc=None, queue_size=None):
    """
    Manages feeder, worker and collector processes

    Arguments:
      feed_func : Data Queue feeder function
      work_func : Worker function
      collect_func : Result Queue collector function
      feed_args : Dictionary of arguments to pass to feed_func
      work_args : Dictionary of arguments to pass to work_func
      collect_args : Dictionary of arguments to pass to collect_func
      nproc : Number of processQueue processes;
              if None defaults to the number of CPUs
      queue_size : Maximum size of the argument queue;
                   if None defaults to 2*nproc

    Returns:
      dict : Dictionary of collector results
    """
    # Define signal handler that raises KeyboardInterrupt
    def _signalHandler(s, f):
        raise SystemExit

    # Define function to terminate child processes
    def _terminate():
        sys.stderr.write('Terminating child processes...')
        # Terminate feeders
        feeder.terminate()
        feeder.join()
        # Terminate workers
        for w in workers:
            w.terminate()
            w.join()
        # Terminate collector
        collector.terminate()
        collector.join
        sys.stderr.write('  Done.\n')

    # Raise SystemExit upon termination signal
    signal.signal(signal.SIGTERM, _signalHandler)

    # Define number of processes and queue size
    if nproc is None:  nproc = mp.cpu_count()
    if queue_size is None:  queue_size = nproc * 2

    # Define shared child process keep alive flag
    alive = mp.Value(ctypes.c_bool, True)

    # Define shared data queues
    data_queue = mp.Queue(queue_size)
    result_queue = mp.Queue(queue_size)
    # TODO:  find out what's up with this context shenanigans
    ctx = mp.get_context()
    collect_queue = ctx.SimpleQueue()
    # Initiate manager and define shared data objects

    try:
        # Initiate feeder process
        feeder = mp.Process(target=feed_func,
                            args=(alive, data_queue),
                            kwargs=feed_args)
        feeder.start()

        # Initiate worker processes
        workers = []
        for __ in range(nproc):
            w = mp.Process(target=work_func,
                           args=(alive, data_queue, result_queue),
                           kwargs=work_args)
            w.start()
            workers.append(w)

        # Initiate collector process
        collector = mp.Process(target=collect_func,
                               args=(alive, result_queue, collect_queue),
                               kwargs=collect_args)
        collector.start()

        # Wait for feeder to finish and add sentinel objects to data_queue
        feeder.join()
        for __ in range(nproc):  data_queue.put(None)

        # Wait for worker processes to finish and add sentinel to result_queue
        for w in workers:  w.join()
        result_queue.put(None)

        # Wait for collector process to finish and add sentinel to collect_queue
        collector.join()
        collect_queue.put(None)

        # Get collector return values
        collected = collect_queue.get()
    except (KeyboardInterrupt, SystemExit):
        sys.stderr.write('Exit signal received\n')
        _terminate()
        sys.exit()
    except Exception as e:
        sys.stderr.write('ERROR:  %s\n' % e)
        _terminate()
        sys.exit()
    except:
        sys.stderr.write('ERROR:  Exiting with unknown exception\n')
        _terminate()
        sys.exit()
    else:
        if not alive.value:
            sys.stderr.write('ERROR:  Exiting due to child process error\n')
            _terminate()
            sys.exit()

    return collected


def feedSeqQueue(alive, data_queue, seq_file, index_func=None, index_args={}):
    """
    Feeds the data queue with SeqRecord objects

    Arguments:
      alive : multiprocessing.Value boolean controlling whether processing
              continues; when False function returns
      data_queue : multiprocessing.Queue to hold data for processing
      seq_file : Sequence file to read input from
      index_func : Function to use to define sequence sets
                   if None do not index sets and feed individual records
      index_args : Dictionary of arguments to pass to index_func

    Returns:
      None
    """
    try:
        # Read input file and index sequence sets if required
        if index_func is None:
            seq_iter = readSeqFile(seq_file)
            data_iter = ((s.id, s) for s in seq_iter)
        else:
            seq_dict = readSeqFile(seq_file, index=True)
            index_dict = index_func(seq_dict, **index_args)
            data_iter = ((k, [seq_dict[i] for i in v]) \
                         for k, v in index_dict.items())
    except:
        alive.value = False
        raise

    try:
        # Iterate over data_iter and feed data queue
        while alive.value:
            # Get data from queue
            if data_queue.full():  continue
            else:  data = next(data_iter, None)
            # Exit upon reaching end of iterator
            if data is None:  break

            # Feed queue
            data_queue.put(SeqData(*data))
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None
    except:
        alive.value = False
        raise

    return None


def processSeqQueue(alive, data_queue, result_queue, process_func, process_args={}):
    """
    Pulls from data queue, performs calculations, and feeds results queue

    Arguments:
      alive : multiprocessing.Value boolean controlling whether processing
              continues; when False function returns
      data_queue : multiprocessing.Queue holding data to process
      result_queue : multiprocessing.Queue to hold processed results
      process_func : function to use for filtering sequences
      process_args : Dictionary of arguments to pass to process_func

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

            # Perform work
            result = process_func(data, **process_args)

            #import cProfile
            #prof = cProfile.Profile()
            #result = prof.runcall(process_func, data, **process_args)
            #prof.dump_stats('worker-%d.prof' % os.getpid())

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


def collectSeqQueue(alive, result_queue, collect_queue, seq_file,
                    task_label, out_args, index_field=None):
    """
    Pulls from results queue, assembles results and manages log and file IO

    Arguments:
      alive : a multiprocessing.Value boolean controlling whether processing
              continues; when False function returns
      result_queue : Multiprocessing.Queue holding worker results
      collect_queue : Multiprocessing.Queue to store collector return values
      seq_file : Sample sequence file name
      task_label : Task label used to tag the output files
      out_args : Common output argument dictionary from parseCommonArgs
      index_field : Field defining set membership for sequence sets
                    if None data queue contained individual records

    Returns:
      None : Adds a dictionary with key value pairs to collect_queue containing
            'log' defining a log object,
            'out_files' defining the output file names
    """
    try:
        # Count records
        if index_field is None:
            result_count = countSeqFile(seq_file)
        else:
            result_count = countSeqSets(seq_file, index_field, out_args['delimiter'])

        # Define output format
        out_type = getFileType(seq_file) if out_args['out_type'] is None \
                   else out_args['out_type']

        # Defined valid alignment output handle
        pass_handle = getOutputHandle(seq_file,
                                      '%s-pass' % task_label,
                                      out_dir=out_args['out_dir'],
                                      out_name=out_args['out_name'],
                                      out_type=out_type)
        # Defined failed alignment output handle
        if out_args['failed']:
            fail_handle = getOutputHandle(seq_file,
                                          '%s-fail'  % task_label,
                                          out_dir=out_args['out_dir'],
                                          out_name=out_args['out_name'],
                                          out_type=out_type)
        else:
            fail_handle = None

        # Define log handle
        if out_args['log_file'] is None:
            log_handle = None
        else:
            log_handle = open(out_args['log_file'], 'w')
    except:
        alive.value = False
        raise

    try:
        # Iterator over results queue until sentinel object reached
        start_time = time()
        set_count = seq_count = pass_count = fail_count = 0
        while alive.value:
            # Get result from queue
            if result_queue.empty():  continue
            else:  result = result_queue.get()
            # Exit upon reaching sentinel
            if result is None:  break

            # Print progress for previous iteration
            printProgress(set_count, result_count, 0.05, start_time)

            # Update counts for current iteration
            set_count += 1
            seq_count += result.data_count

            # Write log
            printLog(result.log, handle=log_handle)

            # Write alignments
            if result:
                pass_count += 1
                SeqIO.write(result.results, pass_handle, out_type)
            else:
                fail_count += 1
                if fail_handle is not None:
                    SeqIO.write(result.data, fail_handle, out_type)
        else:
            sys.stderr.write('PID %s:  Error in sibling process detected. Cleaning up.\n' \
                             % os.getpid())
            return None

        # Print total counts
        printProgress(set_count, result_count, 0.05, start_time)

        # Update return values
        log = OrderedDict()
        log['OUTPUT'] = os.path.basename(pass_handle.name)
        log['SEQUENCES'] = seq_count
        if index_field is not None:
            log['SETS'] = set_count
        log['PASS'] = pass_count
        log['FAIL'] = fail_count
        collect_dict = {'log':log, 'out_files': [pass_handle.name]}
        collect_queue.put(collect_dict)

        # Close file handles
        pass_handle.close()
        if fail_handle is not None:  fail_handle.close()
        if log_handle is not None:  log_handle.close()
    except:
        alive.value = False
        raise

    return None
