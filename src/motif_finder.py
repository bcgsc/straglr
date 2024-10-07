import itertools
import pysam
import re
from operator import itemgetter
from collections import defaultdict
import sys
from itertools import chain, combinations, combinations_with_replacement, permutations

def permute(kmer):
    perms = []
    for i in range(len(kmer)):
            pat = kmer[i:] + kmer[:i]
            perms.append(pat)
    return perms[1:]

def get_kmers(k):
    bases = ['A', 'G', 'T', 'C']
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
    uniq_kmers = set(kmers)
    for kmer in kmers:
        if not kmer in uniq_kmers:
            continue
        perms = permute(kmer)
        uniq_kmers = uniq_kmers - set(perms)
    return list(uniq_kmers)

def find_matches(seq, kmer):
    return [m.start() for m in re.finditer(kmer, seq)]

'''
def subtract_spans(l, spans):
    left = []

    for i in range(len(spans)):
        if i == 0 and spans[i][0] > 0:
            left.append([0, spans[i][0] - 1])

        if i < len(spans) - 1:
            left.append([spans[i][1] + 1, spans[i + 1][0] -1])

        if i == len(spans) - 1 and spans[i][1] < l - 1:
            left.append([spans[i][1] + 1, l - 1])

    return left
'''
def merge_matches(starts, k):
    merged = []
    merged_set = set()
    end = None
    i = 0
    while i < len(starts):
        j = i + 1
        while j < len(starts):
            if starts[j] - starts[j - 1] <= k:
                merged_set.add(starts[j])
                merged_set.add(starts[j - 1])
                end = j
                j += 1
            else:
                break
        if end is not None:
            copies = int((starts[end] + k - starts[i]) / k)
            merged.append([starts[i], starts[end] + k, copies])
            end = None
        i = j

    unmerged = [[i, i + k, 1] for i in sorted(list(set(starts) - merged_set))]
    return sorted(merged + unmerged, key=itemgetter(0))

def get_partials(seq):
    max_del_size = len(seq) - 2
    subseqs = set()
    for d in range(1, max_del_size + 1, 1):
        for i in range(len(seq)):
            if i + d <= len(seq):
                subseq = seq[:i] + seq[i+d:]
                subseqs.add(subseq)
    return sorted(list(subseqs), key=lambda s: (-len(s), s))

def trim_ins(seq):
    def all_sublists(l):
        return chain(*(combinations(l, i) for i in range(len(l) + 1)))
    motifs = set()
    indices = []
    for i in range(len(seq) - 1):
        if seq[i].upper() == seq[i + 1].upper():
            indices.append(i)

    for skips in all_sublists(indices):
        if not skips:
            continue

        motif = ''
        for i in range(len(seq)):
            if i not in skips:
                motif += seq[i]
        motifs.add(motif)

    return sorted(list(motifs), key=lambda s: (-len(s), s))

def merge_blocks(blocks, pat, seq):
    ''' indels '''
    merged = []
    merged_set = set()
    i = 0
    end = None

    partials = None
    trimmed = None

    while i < len(blocks):
        copies = 0
        j = i + 1
        while j < len(blocks):
            to_merge = False
            gap_seq = seq[blocks[j - 1][1]:blocks[j][0]]
            gap_len = blocks[j][0] - blocks[j - 1][1]
            filled = False
            # one base gaps (ins)
            if gap_len < 2:
                filled = True
                to_merge = True
            # small indels
            elif gap_len < len(pat):
                if partials is None:
                    partials = get_partials(pat)
                if gap_seq.upper() in partials:
                    filled = True
                    copies += 1
                    to_merge = True
            
            elif gap_len > len(pat) and gap_len < 2 * len(pat):
                if trimmed is None:
                    trimmed = trim_ins(gap_seq)
                if pat.upper() in trimmed:
                    filled = True
                    copies += 1
                    to_merge = True

            # larger gaps, try fill with combinations of partials
            if not filled and len(pat) > 3 and gap_len > len(pat) and gap_len <= len(pat) * 3:
                if partials is None:
                    partials = get_partials(pat)
                filled = fill_gaps(gap_seq, pat)
                print('gg', gap_len, gap_seq, pat, blocks[j - 1][1], blocks[j][0], filled)
                if filled and filled[1] == 0:
                    filled_start = blocks[j - 1][1] + filled[-2]
                    filled_end = blocks[j - 1][1] + filled[-1]
                    print('gg2', gap_seq, pat, blocks[j - 1], blocks[j], filled_start, filled_end, filled)
                    if filled_start == blocks[j - 1][1] and filled_end == blocks[j][0]:
                        copies += filled[-3]
                        to_merge = True
                    elif blocks[j - 1][1] == filled_start:
                        blocks[j - 1][1] = filled_end
                        blocks[j - 1][2] += filled[-3]
                    elif blocks[j][0] == filled_end:
                        blocks[j][0] = filled_start
                        blocks[j][2] += filled[-3]
                    else:
                        merged.append([filled_start, filled_end, filled[-3]])
            
            if to_merge:
                merged_set.add(tuple(blocks[j]))
                merged_set.add(tuple(blocks[j - 1]))
                end = j
                j += 1
            else:
                break
        if end is not None:
            for n in range(i, end + 1, 1):
                copies += blocks[n][2]
            merged.append([blocks[i][0], blocks[end][1], copies])
            end = None
        i = j

    unmerged = [list(b) for b in sorted(list(set([tuple(b) for b in blocks]) - merged_set))]
    return sorted(merged + unmerged, key=itemgetter(0))

def find_seeds(matches, k, n=2):
    return [m for m in matches if (m[1] - m[0]) / k >= n]

def combine_pats(matches, seq):
    matches_sorted = sorted(matches, key=itemgetter(1,2))
    combined_pats = []
    for i in range(len(matches_sorted)):
        pat, start, end, gap_size = matches_sorted[i]
        #copies = (end - start) / len(pat)
        copies = seq[start:end].upper().count(pat.upper())
        match = pat, start, end, copies, gap_size
        combined_pats.append(match)
    return combined_pats

def combine_pats_old(matches):
    matches_sorted = sorted(matches, key=itemgetter(1,2))
    combined = []
    for i in range(len(matches_sorted)):
        check_before = False
        check_after = False
        pat, start, end, gap_size = matches_sorted[i]
        copies = (end - start) / len(pat)
        #match = pat, start, end, copies

        if i == 0:
            check_before = True
        elif matches_sorted[i][1] >= matches_sorted[i - 1][2]:
            check_before = True
        else:
            diff = matches_sorted[i - 1][2] - matches_sorted[i][1]
            motif_before = matches_sorted[i - 1][0]
            motif_current =  matches_sorted[i][0]
            print('cc1', matches_sorted[i], matches_sorted[i - 1], diff, motif_current, motif_before)
            # current motif is an extension of the previous motif
            if len(motif_current) > len(motif_before):
                print('rr', motif_current[:diff], motif_current[:diff][-1 * len(motif_before):])
                if motif_current[:diff][-1 * len(motif_before):].upper() == motif_before.upper():
                    print('rra')
                    check_before = True
            # previous motif is an extension of the current motif
            else:
                #print('vv1', diff, matches_sorted[i - 1][0][-1 * diff:], matches_sorted[i - 1][0], matches_sorted[i][0])
                if matches_sorted[i - 1][0][-1 * diff:].upper() == matches_sorted[i][0].upper():
                    check_before = True

        if i == len(matches) - 1:
            check_after = True
        elif matches_sorted[i + 1][1] >= matches_sorted[i][2]:
            check_after = True
        else:
            diff = matches_sorted[i][2] - matches_sorted[i + 1][1]
            motif_after = matches_sorted[i + 1][0]
            motif_current =  matches_sorted[i][0]
            print('cc2', matches_sorted[i], matches_sorted[i + 1], diff, motif_current, motif_after)
            if len(motif_current) > len(motif_after):
                print('rr2', motif)
                if matches_sorted[i + 1][0][-1 * diff][:len(matches_sorted[i][0])].upper() == matches_sorted[i][0].upper():
                    print('rr2a')
                    #end -= len(pat)
                    #copies -= 1
                    check_after = True
            elif len(matches_sorted[i + 1][0]) < len(matches_sorted[i][0]):
                #print('vv2', diff, matches_sorted[i][0][-1 * diff:], matches_sorted[i][0], matches_sorted[i + 1][0])
                if matches_sorted[i][0][-1 * diff:].upper() == matches_sorted[i + 1][0].upper():
                    check_after = True

        match = pat, start, end, copies
        if check_before and check_after:
            combined.append(match)

    #print('ww', combined)
    return combined

def fill_gaps(seq, pat):
    #print(seq, pat)
    partials = get_partials(pat)
    n = int(len(seq)/2)
    #print(seq, pat, n, partials)

    choices = []
    #for i in range(4, 5, 1):
    for i in range(2, int(len(seq)/2) + 1, 1):
        for subset in combinations_with_replacement(partials, i):
            if len(''.join(subset)) > len(seq):
                continue
            # homopolymer
            if len(set(subset)) == 1 and len(set(subset[0])) == 1:
                continue

            # start base must be the same
            same_starts = [s for s in subset if s[0].upper() == pat[0].upper()]
            if len(same_starts) != len(subset):
                continue

            for order in permutations(subset):
                #print('ee', subset, order)
                test_seq = ''.join(order)
                if test_seq not in seq:
                    continue
                diff = len(seq) - len(test_seq)
                m = re.search(test_seq, seq)
                choices.append((','.join(order), diff, len(subset), m.start(), m.end()))

    if choices:
        choices.sort(key=itemgetter(1,2))
        #for c in choices:
        #   print(c)
        return choices[0]
    else:
        return None

def merge_blocks_multiple(blocks_start, pat, seq, name=None):
    blocks = blocks_start
    blocks_merged = []
    n = 0
    while blocks_merged != blocks:
        if blocks_merged:
            blocks = blocks_merged 
        blocks_merged = merge_blocks(blocks, pat, seq)
        n += 1
    print('multi', name, pat, n)
    return blocks_merged

def extract_pats(seq, kmers, name=None):
    pats = []
    #print('ss', seq)
    all_seeds = []
    all_matches = {}
    for pat in kmers:
        starts = find_matches(seq, pat)
        if not starts:
            continue
        blocks = merge_matches(starts, len(pat))
        seeds = find_seeds(blocks, len(pat), n=3)
        if not seeds:
            continue
        
        print('uu', name, pat, len(seeds), seeds)

        if len(seeds) > 1:
            start_index = blocks.index(seeds[0])
            end_index = blocks.index(seeds[-1])
            
            blocks_between = blocks[start_index:end_index + 1]
            print('ww', name, start_index, end_index, blocks_between)

            #seeds_merged = merge_blocks(blocks_between, pat, seq)
            #print('qq', name, pat, seeds_merged)

            seeds_merged = merge_blocks_multiple(blocks_between, pat, seq, name=name)
            print('qq', name, pat, seeds_merged)
            seeds = seeds_merged
        
        copies = sum([s[2] for s in seeds])
        start = min([s[0] for s in seeds])
        end = max([s[1] for s in seeds])
        print('zz', name, pat, copies, start, end, seeds)
        
        for i in range(len(seeds)-1):
            ss = seq[seeds[i][1]:seeds[i+1][0]]
            print('ww', name, pat, seeds[i], seeds[i+1], ss)

kmers = []
for k in range(2,7,1):
    kmers.extend(get_kmers(k))
    #for pat in get_kmers(k):
    #    print(pat)
#kmers = get_kmers(4)
#kmers = ['AGGC']

fa = pysam.FastaFile('/projects/btl/rchiu/jira/BTL-2134/test.fa')
name = 'ATXN10.HG01961.p1'
seq = fa.fetch(name)

fa_file = '/projects/btl/rchiu/jira/BTL-2147/tmp/tmptohjf5_5.fa'
fa_file = '/projects/btl/rchiu/jira/BTL-2147/tmp2/tmpdlcik4ss.fa'

fa_file = sys.argv[1]
fa = pysam.FastaFile(fa_file)

#kmers = sorted(get_kmers(5))
for seq_name in fa.references:
    #if not 'chr3:129172577:129172656' in seq_name:
    #if not 'chr3:129172577:129172656' in seq_name or not '9202ed12-74ba-47b6-a7bd-e66a4e3fb901' in seq_name:
    #   continue

    #LRP12
    #if not '012732fc-97e2-463a-b74d-729183b053e1' in seq_name:
    #    continue

    seq = fa.fetch(seq_name)
    combined_pats = extract_pats(seq, kmers, name=seq_name)
    #print(seq_name, combined_pats)
    #break


