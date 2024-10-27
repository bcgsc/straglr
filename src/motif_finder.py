import itertools
import pysam
import re
from operator import itemgetter
from collections import defaultdict, Counter
import sys
from itertools import chain, combinations, combinations_with_replacement, permutations, product
from pybedtools import BedTool
import math

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

def find_seeds(matches, k, n=2):
    return [m for m in matches if (m[1] - m[0]) / k >= n]

def fill_seeds(seeds, matches):
    holes = get_holes(seeds)

    for hole in holes:
        for match in matches:
            if match[0] >= hole[0] and match[1] <= hole[1]:
                seeds.append(match)

    seeds.sort(key=itemgetter(0,1))

def fill_partials(seq, pat):
    partials = get_partials(pat)
    choices = []
    for i in range(2, int(len(seq)/2) + 1, 1):
        for subset in combinations_with_replacement(partials, i):
            if len(''.join(subset)) > len(seq):
                continue
            # homopolymer
            if len(set(subset)) == 1 and len(set(subset[0])) == 1:
                continue
            for order in permutations(subset):
                test_seq = ''.join(order)
                if test_seq not in seq:
                    continue
                diff = len(seq) - len(test_seq)
                m = re.search(test_seq, seq)
                choices.append((order, diff, len(subset), m.start(), m.end()))

    if choices:
        choices.sort(key=itemgetter(1,2))
        return choices[0]
    else:
        return None

def extend_seeds(seeds, matches, pat_len, seq_len, n=2):
    gap_size = n * pat_len

    extended = True
    while extended and seeds[0][0] > 0:
        matches_beyond = [m for m in matches if m[1] < seeds[0][0]]
        if not matches_beyond:
            break

        if matches_beyond[-1][1] <= seeds[0][0] - gap_size:
            seeds.insert(0, matches_beyond[-1])
            extended = True
        else:
            extended = False

    extended = True
    while extended and seeds[-1][1] < seq_len:
        matches_beyond = [m for m in matches if m[0] > seeds[-1][1]]
        if not matches_beyond:
            break

        if matches_beyond[0][0] <= seeds[-1][1] + gap_size:
            seeds.append(matches_beyond[0])
            extended = True
        else:
            extended = False

def blocks_olapped(b1, b2):
    if (b1[0] >= b2[0] and b1[0] <= b2[1]) or (b2[0] >= b1[0] and b2[0] <= b1[1]) or\
        (b1[0] <= b2[0] and b1[1] >= b2[1]) or (b2[0] <= b1[0] and b2[1] >= b1[1]):
        return True
    return False

def fill_repeat(repeat, other_matches, seq, name):
    def make_others_bed():
        bed_str = ''
        for pat in other_matches.keys():
            for match in other_matches[pat]:
                bed_str += '{}\n'.format('\t'.join([name] + list(map(str, match[:3])) + [pat]))
        bed = BedTool(bed_str, from_string=True)
        return bed.sort()

    def make_holes_bed(holes):
        bed_str = ''
        for hole in holes:
            hole_span = tuple(hole)
            hole_seq = seq[hole[0]:hole[1]].upper()
            bed_str += '{}\n'.format('\t'.join([name, str(hole[0]), str(hole[1]), hole_seq]))
        bed = BedTool(bed_str, from_string=True)
        return bed.sort()

    seed_pat_counts = Counter([b[3] for b in repeat if len(b[3]) > 2])
    filled_blocks = []
    holes = get_holes(repeat)
    if not holes:
        return repeat

    others_bed = make_others_bed()
    holes_bed = make_holes_bed(holes)

    candidates = dict([(tuple(h), []) for h in holes])
    filled = dict([(tuple(h), []) for h in holes])
    for cols in holes_bed.intersect(others_bed, wo=True):
        hole_cols = cols[:4]
        candidate_cols = cols[4:9]
        # hole length < candidate motif length
        if len(cols[3]) < len(candidate_cols[-1]):
            continue

        hole_span = int(hole_cols[2]) - int(hole_cols[1])
        candidate_span = int(candidate_cols[2]) - int(candidate_cols[1])

        hole = (int(hole_cols[1]), int(hole_cols[2]))
        if filled[hole]:
            continue

        candidate = list(map(int, candidate_cols[1:4])) + [candidate_cols[-1]]
        # exact span = filled
        if hole_cols[:3] == candidate_cols[:3]:
            filled[hole] = [candidate]
            continue

        # candidate pattern subseq of neighbor "seed" pattern
        if not filled[hole] and candidate_span - hole_span == len(candidate_cols[-1]) and (hole_cols[1] == candidate_cols[1] or hole_cols[2] == candidate_cols[2]):
            candidate_modified = [hole[0], hole[1], candidate[2] - 1, candidate[3]]
            filled[hole] = [candidate_modified]
            continue

        # candidate spap subsume hole = candidate
        candidates[hole].append(candidate)

    new_holes = []
    for hole in filled:
        if not filled[hole] and hole in candidates:
            seed_blocks = []
            for pat, count in seed_pat_counts.most_common():
                blocks = sorted([b for b in candidates[hole] if b[3].upper() == pat.upper() and b[0] >= hole[0] and b[1] <= hole[1]], key=itemgetter(2), reverse=True)
                #print('uu0', hole, pat, count, blocks)
                for block in blocks:
                    if not seed_blocks:
                        seed_blocks.append(block)
                        break
                    olapped = False
                    for s in seed_blocks:
                        #print('uu2', hole, seed_blocks, s, block, blocks_olapped(s, block))
                        if blocks_olapped(s, block):
                        #if (s[0] >= block[0] and s[0] <= block[1]) or (s[1] >= block[0] and s[1] <= block[1]):
                            olapped = True
                            break
                    if not olapped:
                        seed_blocks.append(block)
            #print('uu', hole, seed_blocks)

            if seed_blocks:
                seed_blocks.sort(key=itemgetter(0,1))
                if seed_blocks[0][0] != hole[0]:
                    new_holes.append((hole[0], seed_blocks[0][0]))
                new_holes.extend([tuple(h) for h in get_holes(seed_blocks)])
                if seed_blocks[-1][1] != hole[1]:
                    new_holes.append((seed_blocks[-1][1], hole[1]))
                #print('uu2', hole, seed_blocks, new_holes)
                
                filled[hole] = seed_blocks

    for hole in new_holes:
        filled[hole] = []

    # fill with concatenations of partials
    for hole in filled:
        if not filled[hole]:
            hole_seq = seq[hole[0]:hole[1]].upper()
            if len(hole_seq) <= 2:
                continue
            flanks = [b for b in repeat if b[1] == hole[0] or b[0] == hole[1]]
            if len(flanks) == 2 and flanks[0][3] == flanks[1][3]:
                filled_partials = fill_partials(hole_seq, flanks[0][3])
                if filled_partials and filled_partials[1] == 0 and filled_partials[3] == 0 and filled_partials[4] == len(hole_seq):
                    ppats = filled_partials[0]
                    start = filled_partials[3] + hole[0]
                    new_blocks = []
                    for i in range(len(ppats)):
                        end = start + len(ppats[i])
                        new_blocks.append([start, end, 1, ppats[i]])
                        start = end
                    filled[hole] = new_blocks

    for hole in filled:
        if not filled[hole]:
            hole_seq = seq[hole[0]:hole[1]].upper()
            filled_blocks.append([hole[0], hole[1], 1, hole_seq])
        filled_blocks.extend(filled[hole])

    filled_repeat = sorted(repeat + filled_blocks, key=itemgetter(0))

    fix_fillings(filled_repeat, [s[0] for s in seed_pat_counts.most_common()], name=name)

    #print('nn', name, filled_repeat)
    #print('nn2', name, check_blocks_order(filled_repeat))

    return filled_repeat

def fix_fillings(blocks, seed_pats, name=None):
    partials = {}
    for seed_pat in seed_pats:
        partials[seed_pat] = get_partials(seed_pat)

    changed = False
    last_index = 1
    num_blocks = len(blocks)
    while last_index == 1 or (last_index < num_blocks and changed):
        changed = False
        for i in range(last_index, len(blocks), 1):
            last_index = i
            block1 = blocks[i - 1]
            block = blocks[i]
            if block[3] not in seed_pats:
                #print('qqq', name, block, seed_pats)
                for seed_pat in seed_pats:
                    if block[3] in partials[seed_pat]:
                        if len(seed_pat) - len(block[3]) >= 2:
                            continue
                        changed = True
                        block[3] = seed_pat
                        last_index = i + 1
                if not changed:
                    if block[2] == 1 and len(block[3]) < 10:
                        trimmed = trim_ins(block[3])
                        for pat in trimmed:
                            if pat in seed_pats:
                                block[3] = pat
                                last_index = i + 1
                                changed = True
                                break

            if changed:
                num_blocks = len(blocks)
                break

    # merge
    i = 0
    remove_index = []
    while i < len(blocks) - 1:
        copies = blocks[i][2]
        end = blocks[i][1]
        for j in range(i + 1, len(blocks), 1):
            if blocks[j][3] == blocks[i][3]:
                copies += blocks[j][2]
                end = blocks[j][1]
                remove_index.append(j)
            elif blocks[j][1] - blocks[j][0] == 1 and (blocks[j][3].upper() == blocks[i][3][-1].upper() or blocks[j][3].upper() == blocks[j+1][3][0].upper()):
                #print('bbb', name, blocks[j])
                end = blocks[j][1]
                remove_index.append(j)
            else:
                break
        if copies != blocks[i][2] or end != blocks[i][1]:
            blocks[i][2] = copies
            blocks[i][1] = end
        i = j
    for i in remove_index[::-1]:
        del blocks[i]

def get_seeds(seq, pat, name=None):
    seeds = []
    blocks = []
    starts = find_matches(seq, pat)
    if starts:
        blocks = merge_matches(starts, len(pat))
        seeds = find_seeds(blocks, len(pat), n=3)

    return seeds, blocks

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

def combine_seeds(blocks, seq):
    def remove_olaps(olaps):
        pat_counts = Counter([b[3] for b in blocks])
        olaps.sort(key=itemgetter(2), reverse=True)

        for i in range(len(olaps)):
            for j in range(len(olaps)):
                if i == j or olaps[j][2] == 0:
                    continue
                b1, b2 = sorted([olaps[i], olaps[j]], key=itemgetter(0,1))
                #print('ooo', b1, b2)
                if b2[0] >= b1[1]:
                    continue
                elif b1[0] <= b2[0] and b1[1] >= b2[1]:
                    b2[2] = 0
                elif b1[0] >= b2[0] and b1[1] <= b2[1]:
                    b1[2] = 0
                else:
                    olap = b1[1] - b2[0]
                    if pat_counts[b1[3]] < pat_counts[b2[3]]:
                        change_first = True
                    elif pat_counts[b1[3]] > pat_counts[b2[3]]:
                        change_first = False
                    elif b1[2] < b2[2]:
                        change_first = True
                    elif b1[2] > b2[2]:
                        change_first = False
                    elif pat_counts[b1[3]] > pat_counts[b2[3]]:
                        change_first = False
                    elif pat_counts[b1[3]] < pat_counts[b2[3]]:
                        change_first = True
                    else:
                        change_first = False

                    if change_first:
                        copies_change = math.ceil(olap/len(b1[3]))
                        b1[2] -= copies_change
                        b1[1] -= copies_change * len(b1[3])
                    else:
                        copies_change = math.ceil(olap/len(b2[3]))
                        b2[2] -= copies_change
                        b2[0] += copies_change * len(b2[3])
                    
    num_olaps = 0
    blocks.sort(key=itemgetter(0,1))

    # remove subsumes first
    subsumes = []
    for i in range(len(blocks)):
        for j in range(len(blocks)):
            if i != j:
                if blocks[i][0] > blocks[j][0] and blocks[i][1] < blocks[j][1]:
                    subsumes.append(i)
                    break

    num_olaps += len(subsumes)
    for i in subsumes[::-1]:
        del blocks[i]

    # count overlaps and coverage for picking
    i = 0
    if len(blocks) > 1:
        while i <= len(blocks):
            if i == len(blocks) - 1:
                break
            olaps = [blocks[i]]
            for j in range(i + 1, len(blocks), 1):
                if blocks[j][0] >= max([b[1] for b in olaps]):
                    break
                olaps.append(blocks[j])
                
            if len(olaps) > 1:
                num_olaps += len(olaps) - 1
                remove_olaps(olaps)
            i = j

    empties = [i for i in range(len(blocks)) if blocks[i][2] == 0]
    for i in empties[::-1]:
        del blocks[i]

    coverage = sum([b[1] - b[0] for b in blocks])

    highest_count = Counter([b[3] for b in blocks]).most_common(1)[0][1]

    return num_olaps, coverage, highest_count

def check_blocks_order(blocks):
    problems = []
    for i in range(len(blocks) - 1):
        if blocks[i+1][0] < blocks[i][1]:
            problems.append((blocks[i], blocks[i+1]))
    return problems

def check_blocks_seq(blocks, seq, name=None):
    for block in blocks:
        seq1 = block[3] * block[2]
        seq2 = seq[block[0]:block[1]]
        print('yy', name, block, seq1, seq2, seq1==seq2)

def config_repeat(seeds, seq, name=None):
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
                # need to copy block because coordinate will get changed in combining
                blocks.append(block[:])

        olaps, coverage, highest_count = combine_seeds(blocks, seq)
        olaps_combos.append((olaps, -1 * coverage, -1 * highest_count, blocks))
    
    olaps_combos.sort(key=itemgetter(0, 1, 2))
    '''
    for olaps, coverage, highest_count, combo in olaps_combos:
        pat_counts = Counter([b[3] for b in combo])
        print(olaps, coverage, pat_counts, highest_count, combo)
    '''
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
        
        fill_seeds(seeds, matches)
        extend_seeds(seeds, matches, len(pat), len(seq))
        [add_pat(seed, pat) for seed in [seeds]]
        all_seeds.append(seeds)

    if not all_seeds:
        return None
    
    combined_seeds = config_repeat(all_seeds, seq, name=name)
    print('mm', name, combined_seeds)
    print('mm2', name, check_blocks_order(combined_seeds))

    repeat = fill_repeat(combined_seeds, all_matches, seq, name)

    print('ff', name, repeat)
    print('ff2', name, check_blocks_order(repeat))
    #check_blocks_seq(repeat, seq, name=seq_name)

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
    #if not '61692e52-0988-4fa6-a588-' in seq_name:
    #    continue
    #if not 'd7ac326c-54dc' in seq_name:
    #    continue
    #if not '08db7117' in seq_name:
    #    continue

    seq = fa.fetch(seq_name)
    combined_pats = extract_pats(seq, kmers, name=seq_name)
    #print(seq_name, combined_pats)
    #break


