import psutil
from Bio import SeqIO
import gzip
from tqdm import tqdm
from mpire import WorkerPool
from Bio.SeqRecord import SeqRecord

def monitor_memory(pid):
    process = psutil.Process(pid)
    memory_bytes = process.memory_info().rss
    memory_mb = memory_bytes / (1024 * 1024)  # Convert bytes to megabytes
    return memory_mb

def parse_file_type(data_set):
    items = data_set.split(".")
    ext = items[-1]
    if ext == 'fa' or ext == 'fasta':
        f_type = 'fasta'
    elif ext == 'fq' or ext == 'fastq':
        f_type = 'fastq'
    elif ext == 'gz':
        f_type = items[-2] + '.' + ext
    return f_type
    
def parse_data_index(data_set):
    file_type = parse_file_type(data_set)
    if file_type == 'fastq.gz' or file_type == 'fq.gz' or file_type == 'fa.gz' or file_type == 'fasta.gz':
        ff_type = file_type.split('.')[0]
        handle = gzip.open(data_set, 'rt')
        records = SeqIO.index(handle, ff_type)
        return records, ff_type
    else:
        records = SeqIO.index(data_set, file_type)     
        return records, file_type
    
def parse_data_dict(data_set):
    file_type = parse_file_type(data_set)
    if file_type == 'fastq.gz' or file_type == 'fq.gz' or file_type == 'fa.gz' or file_type == 'fasta.gz':
        ff_type = file_type.split('.')[0]
        handle = gzip.open(data_set, 'rt')
        records_dict = SeqIO.to_dict(SeqIO.parse(handle, ff_type))
        return records_dict, ff_type
    else:
        records_dict = SeqIO.to_dict(SeqIO.parse(data_set, file_type))   
        return records_dict, file_type
    
def covert_fq2fa(fin, fout):
    records_dict, _ = parse_data_dict(fin)
    records, _ = parse_data_index(fin)
    name_lst = list(records)
    fa_records = []
    with WorkerPool(60, shared_objects=records_dict, start_method='fork') as pool:
        with tqdm(total=len(name_lst), desc="Convert fq to fa") as pbar:   
            for tmp_rec in pool.imap(fq2fa, name_lst):
                fa_records.append(tmp_rec)
                pbar.update() 
    with open(fout, "w") as handle:
        SeqIO.write(fa_records, handle, 'fasta')

def fq2fa(shared_objects, name):
    return SeqRecord(shared_objects[name].seq, id=shared_objects[name].id, description=shared_objects[name].description)