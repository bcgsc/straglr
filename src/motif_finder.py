import itertools
import pysam
import re
from operator import itemgetter
from collections import defaultdict, Counter
import sys
from itertools import chain, combinations, combinations_with_replacement, permutations, product
from pybedtools import BedTool

def permute(kmer):
    perms = []
    for i in range(len(kmer)):
            pat = kmer[i:] + kmer[:i]
            perms.append(pat)
    return perms[1:]

def get_kmers(k):
    bases = ['C', 'A', 'G', 'T']
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
    uniq_kmers = set(kmers)
    '''
    for kmer in kmers:
        if not kmer in uniq_kmers:
            continue
        perms = permute(kmer)
        uniq_kmers = uniq_kmers - set(perms)
    '''
    return [k for k in uniq_kmers if len(set(k)) > 1]
    #return list(uniq_kmers)

def find_matches(seq, kmer):
    return [m.start() for m in re.finditer(kmer, seq)]

def get_holes(spans):
    ''' used by fill_repeat() '''
    holes = []

    for i in range(len(spans) - 1):
        j = i + 1
        if spans[i][1] != spans[j][0]:
            holes.append([spans[i][1], spans[j][0]])

    return holes

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

def fill_repeat(repeat, other_matches, seq, name):
    def make_others_bed():
        bed_str = ''
        for pat in other_matches.keys():
            for match in other_matches[pat]:
                bed_str += '{}\n'.format('\t'.join([name] + list(map(str, match[:3])) + [pat]))
        return BedTool(bed_str, from_string=True).sort().saveas('oo.bed')

    def make_holes_bed(holes):
        bed_str = ''
        for hole in holes:
            hole_span = tuple(hole)
            hole_seq = seq[hole[0]:hole[1]].upper()
            bed_str += '{}\n'.format('\t'.join([name, str(hole[0]), str(hole[1]), hole_seq]))
        return BedTool(bed_str, from_string=True).sort().saveas('hh.bed')

    filled_blocks = []
    holes = get_holes(repeat)
    print('hh', name, holes)
    if not holes:
        return filled_blocks

    others_bed = make_others_bed()
    holes_bed = make_holes_bed(holes)

    print('hh0', name, holes)
    candidates = dict([(tuple(h), []) for h in holes])
    filled = dict([(tuple(h), []) for h in holes])
    for cols in holes_bed.intersect(others_bed, wo=True):
        hole_cols = cols[:4]
        candidate_cols = cols[4:9]
        # hole length < candidate motif length
        if len(cols[3]) < len(candidate_cols[-1]):
            #print('aa', cols)
            continue

        hole_span = int(hole_cols[2]) - int(hole_cols[1])
        candidate_span = int(candidate_cols[2]) - int(candidate_cols[1])

        hole = (int(hole_cols[1]), int(hole_cols[2]))
        if filled[hole]:
            continue

        candidate = list(map(int, candidate_cols[1:4])) + [candidate_cols[-1]]
        # exact span = filled
        if hole_cols[:3] == candidate_cols[:3]:
            #filled[hole].append(candidate)
            filled[hole] = [candidate]
            #print(hole, candidate, filled)
            continue

        # candidate pattern subseq of neighbor "seed" pattern
        if not filled[hole] and candidate_span - hole_span == len(candidate_cols[-1]) and (hole_cols[1] == candidate_cols[1] or hole_cols[2] == candidate_cols[2]):
            candidate_modified = [hole[0], hole[1], candidate[2] - 1, candidate[3]]
            print('hv', hole, hole_cols[-1], candidate, candidate_modified)
            filled[hole] = [candidate_modified]
            continue

        # candidate spap subsume hole = candidate
        candidates[hole].append(candidate)
        '''
        elif int(hole_cols[1]) >= int(candidate_cols[1]) and int(hole_cols[2]) <= int(candidate_cols[2]):
            candidates[hole].append(candidate)
        '''

    print('hh1', name, filled)
    for hole in filled:
        if not filled[hole]:
            hole_seq = seq[hole[0]:hole[1]].upper()
            #print('ccc', name, hole, hole_seq, candidates[hole])

        filled_blocks.extend(filled[hole])

    return filled_blocks
    
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

            # all dimers
            lens = list(set([len(s) for s in subset]))
            if len(lens) == 1 and lens[0] <= 2:
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
    #print('multi', name, pat, n)
    return blocks_merged

def extend_blocks(seeds, blocks, pat, seq, before=True, name=None):
    blocks_iter = blocks if not before else blocks[::-1]
    block_bound = seeds[-1] if not before else seeds[0]

    blocks_extended = []
    for block_iter in blocks_iter:
        blocks_test = sorted([block_iter, block_bound], key=itemgetter(0))
        blocks_merged = merge_blocks(blocks_test, pat, seq)
        if blocks_merged != blocks_test:
            #print('yy2', name, pat, seeds, blocks, blocks_test, seq[blocks_test[0][1]:blocks_test[1][0]], pat, blocks_merged)
            blocks_extended.append(block_iter)
            block_bound = block_iter
        else:
            break
    return sorted(blocks_extended, key=itemgetter(0))

def get_seeds(seq, pat, name=None):
    seeds = []
    blocks = []
    starts = find_matches(seq, pat)
    if starts:
        blocks = merge_matches(starts, len(pat))
        seeds = find_seeds(blocks, len(pat), n=3)

        if seeds:
            print('uu', name, pat, len(seeds), seeds)
            start_index = blocks.index(seeds[0])
            end_index = blocks.index(seeds[-1])

            # extensions
            blocks_before = [blocks[i] for i in range(len(blocks)) if i < start_index]
            blocks_after = [blocks[i] for i in range(len(blocks)) if i > end_index]
            blocks_extended_before = []
            blocks_extended_after = []
            if blocks_before:
                blocks_extended_before = extend_blocks(seeds, blocks_before, pat, seq, before=True, name=name)
                #if blocks_extended_before:
                #    print('yy before', name, pat, seeds, blocks_extended_before, start_index, blocks.index(blocks_extended_before[0]))
            if blocks_after:
                blocks_extended_after = extend_blocks(seeds, blocks_after, pat, seq, before=False, name=name)
                #if blocks_extended_after:
                #    print('yy after', name, pat, seeds, blocks_extended_after, end_index, blocks.index(blocks_extended_after[-1]))
       
            #if blocks_extended_before or blocks_extended_after:
            #    print('ee', name, pat, seeds, 'bf', blocks_extended_before, 'af', blocks_extended_after, 'new', blocks_extended_before + seeds + blocks_extended_after)
            seeds = blocks_extended_before + seeds + blocks_extended_after
            start_index = blocks.index(seeds[0])
            end_index = blocks.index(seeds[-1])

            if len(seeds) > 1:
                blocks_between = blocks[start_index:end_index + 1]
                #print('ww', name, start_index, end_index, blocks_between)

                seeds_merged = merge_blocks_multiple(blocks_between, pat, seq, name=name)
                #print('qq', name, pat, seeds_merged)
                seeds = seeds_merged

            copies = sum([s[2] for s in seeds])
            start = min([s[0] for s in seeds])
            end = max([s[1] for s in seeds])
            print('zz', name, pat, copies, start, end, seeds)

    return seeds, blocks

def split_seeds(seeds, f=0.1):
    groups = []
    starts = [0]
    seeds_span = seeds[-1][1] - seeds[0][0]
    max_gap_size = f * seeds_span
    for i in range(len(seeds)-1):
        gap_size = seeds[i+1][0] - seeds[i][1]
        if gap_size > max_gap_size:
            starts.append(i+1)
    starts.append(len(seeds))
    for i in range(len(starts)):
        if i == len(starts) - 1:
             break
        group = seeds[starts[i]:starts[i+1]]
        groups.append(group)
    return groups

def find_perms(seeds, name=None):
    seeds_copies = []
    for seed in seeds:
        copies = sum([m[2] for m in seed])
        seeds_copies.append((copies, seed))
    seeds_copies.sort(key=itemgetter(0), reverse=True)
    
    i = 0
    groups = []
    group = []
    while i < len(seeds_copies):
        seed1 = seeds_copies[i][1]
        pat1 = seed1[0][-1]
        group = [seed1]
        perms = permute(pat1)
        if i == len(seeds_copies) - 1:
            break
        for j in range(i+1, len(seeds_copies), 1):
            seed2 = seeds_copies[j][1]
            pat2 = seed2[0][-1]

            if pat2 in perms:
                group.append(seed2)
                if j == len(seeds_copies) - 1:
                    j += 1
            else:
                break
        i = j
        groups.append(group)
        group = []
    if group:
        groups.append(group)

    return groups

def combine_blocks(blocks):
    blocks_sorted = sorted(blocks, key=itemgetter(0,1))
    combined = [blocks_sorted[0]]
    olaps = 0
    for i in range(1, len(blocks_sorted), 1):
        b1 = combined[-1]
        b2 = blocks_sorted[i]
        if b2[0] >= b1[1]:
            combined.append(b2)
            #print('qqq normal', b1, b2)
            
        # b1 subsume b2, replace b1 by b2
        elif b2[0] == b1[0]:
            #print('qqq subsume', b1, b2)
            combined[-1] = b2
            olaps += 1
        # b1 overlap b2, adjust coord
        elif b2[0] > b1[0]:
            #print('qqq overlap', b1, b2)
            olaps += 1

            if b1[2] < b2[2]:
                b1[2] -= 1
                b1[1] -= len(b1[3])
                #print('qqq1a', b1, b2)
                if b1[2] == 0:
                    combined[-1] = b2
                else:
                    combined.append(b2)

            elif b2[2] < b1[2]:
                b2[2] -= 1
                b2[0] += len(b2[3])
                #print('qqq1b', b1, b2)
                if b2[2] > 0:
                    combined.append(b2)

            elif b1[2] == 1:
                combined[-1] = b2

            elif b2[2] == 1:
                continue

            else:
                continue

        else:
            print('qqq ??', b1, b2)
            combined.append(b2)

    coverage = sum([b[1] - b[0] for b in combined])
    return olaps, coverage, combined

def config_repeat(seeds, name=None):
    groups = find_perms(seeds, name=name)
    
    group_indexes = []
    for group in groups:
        group_indexes.append(list(range(len(group))))
    combos = list(product(*group_indexes))
    olaps_combos = []
    for combo in combos:
        blocks = []
        for index, group in zip(combo, groups):
            for block in group[index]:
                if block[2] <= 0:
                    print('mm', block)
                # need to copy block because coordinate will get changed in combining
                blocks.append(block[:])

        olaps, coverage, combined = combine_blocks(blocks)
        olaps_combos.append((olaps, -1 * coverage, combined))

        for block in combined:
            if block[2] <= 0:
                print('nn', block)
    
    olaps_combos.sort(key=itemgetter(0, 1))
    for olaps, coverage, combo in olaps_combos:
        pat_counts = Counter([b[3] for b in combo])
        print('jj', name, olaps, -1 * coverage, pat_counts, len(combo), combo)

    return olaps_combos[0][-1]

def extract_pats(seq, kmers, name=None):
    def add_pat(matches, pat):
        [m.append(pat) for m in matches]
    
    all_matches = {}
    all_repeats = []
    all_seeds = []
    for pat in kmers:
        seeds, matches = get_seeds(seq, pat, name=name)
        all_matches[pat] = matches
        if not seeds:
            continue
        
        [add_pat(seed, pat) for seed in [seeds]]
        all_seeds.append(seeds)
        continue
        # break up seeds that are too far apart (not necessary if only no flanking sequences given)
        # skip dimer seeds
        if len(pat) > 2:
            seeds_split = split_seeds(seeds)
            [add_pat(seed, pat) for seed in seeds_split]
            all_repeats.extend(seeds_split)

    combined_seeds = config_repeat(all_seeds, name=name)
    print('mm', name, combined_seeds)
    #print('nn', name, [seq[s[0]:s[1]] for s in combined_seeds])
    filled_holes = fill_repeat(combined_seeds, all_matches, seq, name)
    repeat = sorted(combined_seeds + filled_holes, key=itemgetter(0))
    print('ff', name, repeat)
    return repeat

kmers = []
max_k = 6
exclude = []
for k in range(2, max_k + 1,1):
    multiples = []
    for m in range(k + 1, max_k + 1, 1):
        if m % k == 0 and k != max_k:
            multiples.append(int(m/k))
    ks = get_kmers(k)
    if multiples:
        for m in multiples:
            exclude.extend([k * m for k in ks])
    kmers.extend(ks)
    #for pat in get_kmers(k):
    #    print(pat)
kmers = sorted(list(set(kmers) - set(exclude)))
#kmers = get_kmers(3)
#kmers = ['AGGC']
'''
fa = pysam.FastaFile('/projects/btl/rchiu/jira/BTL-2134/test.fa')
name = 'ATXN10.HG01961.p1'
seq = fa.fetch(name)

fa_file = '/projects/btl/rchiu/jira/BTL-2147/tmp/tmptohjf5_5.fa'
fa_file = '/projects/btl/rchiu/jira/BTL-2147/tmp2/tmpdlcik4ss.fa'
'''
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
    if not '61692e52-0988-4fa6-a588-' in seq_name:
        continue

    seq = fa.fetch(seq_name)
    combined_pats = extract_pats(seq, kmers, name=seq_name)
    #print(seq_name, combined_pats)
    #break


