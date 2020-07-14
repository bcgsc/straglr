from pathos.multiprocessing import ProcessingPool as Pool
import tempfile
import sys
import os

def split_tasks(args, n):
    k, m = divmod(len(args), n)
    return (args[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))
                
def parallel_process(func, args, nprocs, bam=None, fasta=None):
    p = Pool(nprocs)
    if bam is not None and fasta is not None:
        results = p.map(func, args, [bam]*len(args), [fasta]*len(args))
    elif bam is not None and fasta is None:
        results = p.map(func, args, [bam]*len(args))
    elif bam is None and fasta is not None:
        results = p.map(func, args, [fasta]*len(args))
    else:
        results = p.map(func, args)
    return results

def combine_batch_results(batch_results, data_type):
    all_results = None
    if data_type is list or data_type is tuple:
        all_results = []
        for batch_result in batch_results:
            all_results.extend(batch_result)
    elif data_type is dict:
        all_results = {}
        for batch_result in batch_results:
            for key, val in batch_result.items():
                all_results[key] = val

    return all_results

def create_tmp_file(content):
    fd, path = tempfile.mkstemp()
    try:
        with os.fdopen(fd, 'w') as out:
            out.write(content)
    except:
        sys.exit("can't generate temp file: {}".format(path))

    return path

def reverse_complement(seq):
    """Reverse complements sequence string"""
    complement = str.maketrans("agtcAGTC", "tcagTCAG")
    return seq[::-1].translate(complement)
