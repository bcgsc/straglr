from pathos.multiprocessing import ProcessingPool as Pool
import tempfile
import sys
import os
from operator import itemgetter

def split_tasks(args, n):
    k, m = divmod(len(args), n)
    batches = (args[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))
    return [b for b in batches if b]

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

def create_tmp_file(content, ext=None):
    fd, path = tempfile.mkstemp(suffix=ext)
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

def merge_spans(spans):
    spans_checked = [span for span in spans if len(span) == 2 and type(span[0]) is int and type(span[1]) is int and span[0]<=span[1]]
    spans_sorted = sorted(spans_checked, key=itemgetter(0,1))
    spans_merged = [list(spans_sorted[0])]
    for i in range(1, len(spans_sorted)):
        if spans_sorted[i][0] <= spans_merged[-1][1] + 1:
            if spans_sorted[i][1] > spans_merged[-1][1]:
                spans_merged[-1][1] = spans_sorted[i][1]
        else:
            spans_merged.append(list(spans_sorted[i]))

    return spans_merged

def complement_spans(spans):
    spans_merged = merge_spans(spans)
    complement = []
    for i in range(1, len(spans_merged)):
        complement.append([spans[i-1][1] + 1, spans[i][0] - 1])
    return complement

def is_same_repeat(reps, same_pats=None, min_fraction=0.5):
        def check_same_pats(rep1, rep2):
            if rep1 in same_pats:
                if rep2 in same_pats[rep1]:
                    return True

                for rep in same_pats[rep1]:
                    if rep2 in rep:
                        if float(rep.count(rep2) * len(rep2)) / len(rep) >= min_fraction:
                            return True
                    elif rep in rep2:
                        if float(rep2.count(rep) * len(rep)) / len(rep2) >= min_fraction:
                            return True

            return False

        if len(reps[0]) <= len(reps[1]):
            rep1, rep2 = reps[0], reps[1]
        else:
            rep2, rep1 = reps[0], reps[1]

        if rep1 == rep2:
            return True

        perms1 = []
        for i in range(len(rep1)):
            pat = rep1[i:] + rep1[:i]
            perms1.append(pat)

        for p1 in perms1:
            if p1 in rep2:
                if float(rep2.count(p1) * len(p1)) / len(rep2) >= min_fraction:
                    return True

        if same_pats:
            if check_same_pats(reps[0], reps[1]) or check_same_pats(reps[1], reps[0]):
                return True

            for i in range(0, len(reps[0]), 5):
                pat = reps[0][i:] + reps[0][:i]
                if check_same_pats(pat, reps[1]) or check_same_pats(reps[1], pat):
                    return True

            for i in range(0, len(reps[1]), 5):
                pat = reps[1][i:] + reps[1][:i]
                if check_same_pats(pat, reps[0]) or check_same_pats(reps[0], pat):
                    return True

        return False

