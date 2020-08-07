import pysam
from distutils import spawn
import sys
import subprocess
import os
import re
from intspan import intspan
from .variant import Variant, Allele
from collections import defaultdict, Counter
from operator import itemgetter, attrgetter
import itertools
from .utils import split_tasks, parallel_process, combine_batch_results, create_tmp_file, reverse_complement
from .ins import INSFinder, INS
import math
import random
from pybedtools import BedTool

class TREFinder:
    def __init__(self, bam, genome_fasta, reads_fasta=None, check_split_alignments=True,
                 max_str_len=50, min_str_len=2, flank_size=100, min_support=2, nprocs=1,
                 clustering='gmm', max_num_clusters=3, genotype_in_size=False, eps=None,
                 remove_tmps=False):
        self.bam = bam
        self.genome_fasta = genome_fasta
        trf_path = spawn.find_executable("trf")
        if not trf_path:
            sys.exit('ABORT: {}'.format("can't find trf in PATH"))
            
        #self.trf_args = '2 5 5 80 10 50 500 -d -h'
        self.trf_args = '2 5 5 80 10 10 500 -d -h'
        self.flank_len = 2000

        self.reads_fasta = reads_fasta

        # for checking sequences flanking repeat
        self.trf_flank_min_mapped_fraction = 0.7
        self.trf_flank_size = flank_size
        self.nprocs = nprocs

        self.check_split_alignments = check_split_alignments

        self.update_loci = True

        self.max_str_len = max_str_len
        self.min_str_len = min_str_len

        self.min_support = min_support

        self.tmp_files = set()

        # Cluster object for genotyping
        self.clustering = clustering
        self.max_num_clusters = max_num_clusters

        # report genotype in size instead of copy numbers (default)
        self.genotype_in_size = genotype_in_size

        # epsilon(distance between points) parameter for dbscan clustering
        self.eps = eps

        # True when running in genotyping mode - strictly genotyping within given coordinates
        self.strict = False

        self.remove_tmps = remove_tmps

    def construct_trf_output(self, input_fasta):
        m = re.search('(\d[\d\s]*\d)', self.trf_args)
        if m is not None:
            return '{}/{}.{}.dat'.format(os.getcwd(), os.path.basename(input_fasta), m.group(1).replace(' ', '.'))
    
    def run_trf(self, input_fasta):
        cmd = ' '.join(['trf', input_fasta, self.trf_args])
        # redirect stdout and stderr to devnull
        FNULL = open(os.devnull, 'w')
        returncode = subprocess.call(cmd, shell=True, stdout=FNULL, stderr=FNULL)
        
        output = self.construct_trf_output(input_fasta)
        if os.path.exists(output):
            return output
        else:
            sys.exit('cannot run {}'.format(cmd))

    def type_trf_cols(self, cols):
        return list(map(int, cols[:3])) + [float(cols[3])] + list(map(int, cols[4:12])) + [float(cols[12])] + cols[13:]
    
    def parse_trf(self, trf_output):
        """
            columns:
            start_index, end_index, periold_size, copy_number, consensus_size,
            percent_matches, percent_index, score, A, C, G, T, entropy,
            consensus_pattern

            start and end index are 1-based
        """
        self.tmp_files.add(trf_output)
        results = defaultdict(list)
        with open(trf_output, 'r') as ff:
            for line in ff:
                cols = line.rstrip().split()
                if not cols:
                    continue
                if cols[0] == 'Sequence:':
                    seq = cols[1]
                elif len(cols) == 15:
                    if len(cols[13]) >= self.min_str_len and len(cols[13]) <= self.max_str_len:
                        results[seq].append(self.type_trf_cols(cols))

        return results

    def extract_tres(self, ins_list, target_flank=2000):
        genome_fasta = pysam.Fastafile(self.genome_fasta)

        # prepare input for trf
        trf_input = ''
        for ins in ins_list:
            target_start, target_end = ins[1] - target_flank + 1, ins[1] + target_flank
            prefix = '.'.join(map(str, [INS.eid(ins), len(ins[5]), target_start, target_end]))
            trf_input += '>{}.q\n{}\n'.format(prefix, ins[7])
            trf_input += '>{}.t\n{}\n'.format(prefix, self.extract_genome_neighbour(ins[0],
                                                                                    ins[1],
                                                                                    target_flank,
                                                                                    genome_fasta))
            trf_input += '>{}.i\n{}\n'.format(prefix, ins[5])

        results = self.perform_trf(trf_input)

        grouped_results = {}
        for seq in results.keys():
            read, read_type = seq.rsplit('.', 1)
            if not read in grouped_results:
                grouped_results[read] = {}
            grouped_results[read][read_type] = results[seq]

        expansions = self.analyze_trf(grouped_results, target_flank)
        print('uue', expansions)

        if self.remove_tmps:
            self.cleanup()

        return expansions

    def extract_genome_neighbour(self, chrom, tpos, w, genome_fasta):
        return genome_fasta.fetch(chrom, tpos - w, tpos + w)
    
    def analyze_trf(self, results, target_flank, full_cov=0.7):
        same_pats = self.find_similar_long_patterns_ins(results)

        expansions = {}
        for seq_id in sorted(results.keys()):
            eid, ins_len, gstart, gend = seq_id.rsplit('.', 3)

            if 't' in results[seq_id]:
                print('uu1', seq_id)
                pat, pgstart, pgend = self.analyze_trf_per_seq(results[seq_id], int(ins_len), int(gstart), int(gend), same_pats, target_flank, full_cov=full_cov)
                if pat is not None:
                    expansions[eid] = pat, pgstart, pgend

        return expansions
    
    def find_similar_long_patterns_ins(self, results, min_len=4):
        all_pats = {}
        for seq_id in results.keys():
            all_pats[seq_id] = {'query': set(), 'target': set()}
            for seq_type in results[seq_id].keys():
                for result in results[seq_id][seq_type]:
                    if len(result[13]) >= min_len:
                        if seq_type == 'i':
                            all_pats[seq_id]['query'].add(result[13])
                        else:
                            all_pats[seq_id]['target'].add(result[13])

        same_pats = {}
        for seq_id in all_pats.keys():
            queries = set()
            targets = set()
            if all_pats[seq_id]['query'] and all_pats[seq_id]['target']:
                queries |= all_pats[seq_id]['query']
                targets |= all_pats[seq_id]['target']
                blastn_out = self.align_patterns(queries, targets, word_size=4)
                if os.path.exists(blastn_out):
                    results = self.parse_pat_blastn(blastn_out)
                    if results:
                        for query in results.keys():
                            same_pats[query] = results[query]

        return same_pats

    def is_same_repeat(self, reps, same_pats=None, min_fraction=0.5):
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
            if reps[0] in same_pats and reps[1] in same_pats[reps[0]]:
                return True

            if reps[1] in same_pats and reps[0] in same_pats[reps[1]]:
                return True

            for i in range(0, len(reps[0]), 5):
                pat = reps[0][i:] + reps[0][:i]
                if pat in same_pats and reps[1] in same_pats[pat]:
                    return True

            for i in range(0, len(reps[1]), 5):
                pat = reps[1][i:] + reps[1][:i]
                if pat in same_pats and reps[0] in same_pats[pat]:
                    return True

        return False

    def combine_trf_coords(self, coords, bounds, buf=20, max_sep=50):
        # screen out repeat completely in the flanks
        coords = [c for c in coords if not (c[1] < bounds[0] - buf or c[0] > bounds[1] + buf)]
        if not coords:
            return []

        merged = intspan.from_ranges(coords)
        gaps = merged.complement()

        merged_list = merged.ranges()
        gap_list = gaps.ranges()
        gaps_filled = []
        for i in range(len(merged_list)-1):
            if gap_list[i]:
                gap_size = gap_list[i][1] - gap_list[i][0] + 1
                if gap_size <= max_sep:
                    gaps_filled.append(gap_list[i])

        return intspan.from_ranges(merged_list + gaps_filled).ranges()

    def analyze_trf_per_seq(self, result, ins_len, gstart, gend, same_pats, target_flank, full_cov=0.8, mid_pt_buf=200):
        pattern_matched = None
        pgstart, pgend = None, None

        # filter by locations
        filtered_results = {'i':[], 'q':[], 't':[]}
        filtered_patterns = {'i':[], 'q':[], 't':[]}
        if 'i' in result:
            filtered_results['i'] = [r for r in result['i'] if float(r[1] - r[0] + 1) / ins_len >= full_cov]
            filtered_patterns['i'] = [r[13] for r in filtered_results['i']]

        mid_pts = self.flank_len + mid_pt_buf, self.flank_len + 1 - mid_pt_buf
        for pt in ('q', 't'):
            if pt in result:
                filtered_results[pt] = [r for r in result[pt] if r[0] <= mid_pts[0] and r[1] >= mid_pts[1]]
                filtered_patterns[pt] = [r[13] for r in filtered_results[pt]]

        for i_result in sorted(filtered_results['i'], key=lambda r: len(r[13])):
            i_pat = i_result[13]
            if len(i_pat) == 1:
                continue
            i_pat_matches = {'q': False, 't': False}
            for pt in ('t', 'q'):
                for r in sorted(filtered_results[pt], key=lambda r: len(r[13]), reverse=True):
                    pat = r[13]
                    if i_pat in pat or pat in i_pat or self.is_same_repeat((i_pat, pat), same_pats):
                        i_pat_matches[pt] = True
                        break

            if i_pat_matches['q'] or i_pat_matches['t']:
                pattern_matched = i_pat

                if i_pat_matches['t']:
                    dist_from_mid = {}
                    for i in range(len(filtered_results['t'])):
                        result = filtered_results['t'][i]
                        dist_from_mid[i] = min(abs(float(result[0]) - target_flank), abs(float(result[1]) - target_flank))

                    dist_from_mid_picked = None
                    # check for identical repeats first
                    for i in range(len(filtered_results['t'])):
                        r = filtered_results['t'][i]
                        print('uui', r, i_pat, self.is_same_repeat((i_pat, r[13]), min_fraction=1), dist_from_mid[i], dist_from_mid_picked)
                        if self.is_same_repeat((i_pat, r[13]), min_fraction=1):
                            if dist_from_mid_picked is None or dist_from_mid[i] < dist_from_mid_picked:
                               pgstart, pgend = gstart + r[0] - 1, gstart + r[1] - 1
                               dist_from_mid_picked = dist_from_mid[i]
                               print('uu', r, i_pat, pgstart, pgend, pattern_matched)
            
                    if pgstart is None:
                        for i in range(len(filtered_results['t'])):
                            r = filtered_results['t'][i]
                            print('uub', r, i_pat, self.is_same_repeat((i_pat, r[13]), same_pats), dist_from_mid[i], dist_from_mid_picked)
                            if i_pat in r[13] or r[13] in i_pat or self.is_same_repeat((i_pat, r[13]), same_pats):
                                if dist_from_mid_picked is None or dist_from_mid[i] < dist_from_mid_picked:
                                    pgstart, pgend = gstart + r[0] - 1, gstart + r[1] - 1
                                    dist_from_mid_picked = dist_from_mid[i]
                                    print('uu', r, i_pat, pgstart, pgend, pattern_matched)

                break

        # reference and query(flanking) may not have the repeat long enough for trf to detect
        if not pattern_matched and filtered_patterns['i']:
            # if target patterns don't match, but there is, just use longest one
            if filtered_results['t']:
                tpats = [r for r in filtered_results['t'] if int(r[0]) < target_flank and int(r[1]) > target_flank]
                if tpats:
                    tpats_sorted = sorted(tpats, key=lambda p:len(p[-1]), reverse=True)
                    pgstart, pgend = gstart + tpats_sorted[0][0] - 1, gstart + tpats_sorted[0][1] - 1
            if not pgstart and not pgend:
                # just take the pattern with the longest repeat
                candidates = sorted([r for r in filtered_results['i'] if r[13] in filtered_patterns['i']], key=lambda r:len(r[-1]), reverse=True)
                pattern_matched = candidates[0][13]
                # deduce the insertion pos from gstart and gend, and use it as the repeat insertion point
                mid = float(gstart + gend) / 2
                pgstart, pgend = int(math.floor(mid)), int(math.ceil(mid))

        print('uu3', pattern_matched, pgstart, pgend)
        return pattern_matched, pgstart, pgend

    def annotate(self, ins_list, expansions):
        #ins_dict = dict((ins.eid, ins) for ins in ins_list)
        ins_dict = dict((INS.eid(ins), ins) for ins in ins_list)
        for eid in expansions.keys():
            if eid in ins_dict:
                ins_dict[eid][6] = 'tre'
                ins_dict[eid].append(expansions[eid][0])

                if expansions[eid][1] is not None and expansions[eid][2] is not None:
                    ins_dict[eid][1] = expansions[eid][1]
                    ins_dict[eid][2] = expansions[eid][2]

    def merge_loci(self, tres, d=50):
        loci_sorted = sorted([[tre[0], tre[1], tre[2], tre[8], len(tre[8])] for tre in tres], key=itemgetter(0, 1, 4))

        bed_lines = []
        for tre in tres:
            bed_lines.append('\t'.join(map(str, [tre[0], tre[1], tre[2], tre[8]])))
        bed_file = create_tmp_file('\n'.join(bed_lines))
        print('tre loci: {}'.format(bed_file))

        tres_bed = BedTool(bed_file)
        tres_merged = tres_bed.sort().merge(d=d, c='4', o='distinct')

        merged = []
        for tre in tres_merged:
            repeats = sorted(tre[-1].split(','), key=len)
            if len(repeats[0]) < self.min_str_len:
                continue
            merged.append([tre[0], int(tre[1]), int(tre[2]), tre[-1]])

        return merged

    def extract_aln_tuple(self, aln, tcoord, search_direction, max_extend=200):
        found = []
        tpos = tcoord
        if search_direction == 'left':
            while not found and tpos >= tcoord - max_extend:
                found = [p for p in aln.aligned_pairs if p[1] == tpos and p[0] is not None]
                tpos -= 1
        else:
            while not found and tpos <= tcoord + max_extend:
                found = [p for p in aln.aligned_pairs if p[1] == tpos and p[0] is not None]
                tpos += 1

        if found:
            return found[0]
        else:
            return found

    def extract_aln_tuple2(self, aln, tcoord, search_direction, max_extend=200):
        tuples_in_range = []
        if search_direction == 'left':
            tstart = tcoord - max_extend
            tend = tcoord

        else:
            tstart = tcoord
            tend = tcoord + max_extend

        tuples_in_range = [p for p in aln.aligned_pairs if p[1] is not None and p[1] >= tstart and p[1] <= tend and p[0] is not None]


    def extract_subseq(self, aln, target_start, target_end, reads_fasta=None, max_extend=50):
        tstart = None
        tend = None
        qstart = None
        qend = None
        seq = None

        self.extract_aln_tuple2(aln, target_start, 'left')
        self.extract_aln_tuple2(aln, target_end, 'right')

        if target_start >= aln.reference_start and target_start <= aln.reference_end:
            aln_tuple = self.extract_aln_tuple(aln, target_start, 'left', max_extend=max_extend)
            if aln_tuple:
                qstart, tstart = aln_tuple

        if target_end >= aln.reference_start and target_end <= aln.reference_end:
            aln_tuple = self.extract_aln_tuple(aln, target_end, 'right', max_extend=max_extend)
            if aln_tuple:
                qend, tend = aln_tuple

        if qstart is not None and qend is not None:
            if not reads_fasta:
                seq = aln.query_sequence[qstart:qend]
            else:
                seq = INSFinder.get_seq(reads_fasta, aln.query_name, aln.is_reverse, [qstart, qend])

        return seq, tstart, tend, qstart

    def create_trf_fasta(self, locus, read, tstart, tend, qstart, seq, read_len):
        """ for genotyping """
        fields = list(locus) + [read, tstart, tend, len(seq), qstart, read_len]
        header = '{}'.format(':'.join(map(str, fields)))

        return header, '>{}\n{}\n'.format(header, seq)

    def examine_repeats(self, seq, repeat, max_sep=100, min_cov=0.8):
        """ for regex extracting, return the most common motif """
        pat = re.compile(repeat.upper().replace('*', '[AGCT]'))
        pat_counts = Counter()
        coords = []
        for m in pat.finditer(seq):
            if not coords or m.start() - coords[-1] - 1 <= max_sep:
                coords.extend([m.start(), m.start() + len(repeat) - 1])
                pattern = seq[m.start():m.start() + len(repeat)]
                pat_counts.update([pattern])
        sorted_coords = sorted(coords)
        return float(sorted_coords[-1] - sorted_coords[0] + 1) / len(seq) >= min_cov,\
               (sorted_coords[0], sorted_coords[-1]),\
               set([pat_counts.most_common()[0][0]])

    def perform_trf(self, seqs):
        trf_fasta = create_tmp_file(seqs)
        print('trf input {}'.format(trf_fasta))
        self.tmp_files.add(trf_fasta)
        output = self.run_trf(trf_fasta)
        results = self.parse_trf(output)
        return results

    def find_similar_long_patterns_gt(self, results, patterns, min_len=20):
        same_pats = {}
        queries = defaultdict(set)
        targets = defaultdict(set)
        for seq in results.keys():
            cols = seq.split(':')
            if len(cols) < 7:
                continue
            locus = tuple(cols[:3])
            same_pats[locus] = None
            targets[locus] |= set(patterns[seq].split(','))
            for result in results[seq]:
                if len(result[13]) >= min_len:
                    queries[locus].add(result[13])

        for locus in queries.keys():
            if locus in queries and locus in targets:
                blastn_out = self.align_patterns(queries[locus], targets[locus], locus=locus)
                if os.path.exists(blastn_out):
                    same_pats[locus] = self.parse_pat_blastn(blastn_out)

        return same_pats

    def align_patterns(self, queries, targets, locus=None, word_size=9):
        query_fa = ''
        for seq in queries:
            for i in range(0, len(seq), 5):
                pat = seq[i:] + seq[:i]
                query_fa += '>{}\n{}\n'.format(pat, pat)

        target_fa = ''
        for seq in targets:
            target_fa += '>{}\n{}\n'.format(seq, seq)

        query_file = create_tmp_file(query_fa)
        target_file = create_tmp_file(target_fa)
        blastn_out = create_tmp_file('')
        print('blastn {} {} {} {}'.format(locus, query_file, target_file, blastn_out))
        self.tmp_files.add(query_file)
        self.tmp_files.add(target_file)
        self.tmp_files.add(blastn_out)

        cmd = ' '.join(['blastn',
                        '-query',
                        query_file,
                        '-subject',
                        target_file,
                        '-task blastn -word_size {} -outfmt 6 -out'.format(word_size),
                        blastn_out])
        print(cmd)
        # redirect stdout and stderr to devnull
        FNULL = open(os.devnull, 'w')
        returncode = subprocess.call(cmd, shell=True, stdout=FNULL, stderr=FNULL)

        if os.path.exists(blastn_out):
            return blastn_out
        else:
            sys.exit('cannot run {}'.format(cmd))

    def parse_pat_blastn(self, blastn_out, min_pid=0.8, min_alen=0.8):
        matches = defaultdict(set)
        with open(blastn_out, 'r') as ff:
            for line in ff:
                cols = line.rstrip().split('\t')
                query = cols[0]
                subject = cols[1]
                qlen = len(query)
                slen = len(subject)
                pid = float(cols[2])
                alen = int(cols[3])

                if pid >= min_pid and float(alen)/min(qlen, slen) >= 0.8:
                    matches[query].add(subject)

        return matches

    def extract_alleles_trf(self, trf_input, repeat_seqs, flank, clipped, bam, strands, patterns, min_span=0.6, too_close_to_read_end=200):
        # min_span used to be 0.8
        results = self.perform_trf(trf_input)
        same_pats = self.find_similar_long_patterns_gt(results, patterns)

        # group by locus
        alleles = defaultdict(dict)

        for seq in results.keys():
            cols = seq.split(':')
            if len(cols) < 7:
                print('problematic seq id: {}'.format(seq))
                continue
            locus = tuple(cols[:3])

            expected_pats = patterns[seq].split(',')
            expected_pat_sizes = [len(p) for p in expected_pats]

            read = ':'.join(cols[3:-5])
            gstart, gend = cols[-5], cols[-4]
            rstart = int(cols[-2])
            read_len = int(cols[-1])

            results_matched = []
            for result in results[seq]:
                if len(result[13]) >= self.min_str_len and len(result[13]) <= self.max_str_len:
                    for pat in expected_pats:
                        if self.is_same_repeat((result[13], pat), same_pats=same_pats[locus]):
                            results_matched.append(result)
                            continue

            if results_matched:
                #seq_len = int(cols[-1])
                seq_len = int(cols[-3])
                bounds = (flank, seq_len - flank)
                combined_coords = self.combine_trf_coords([(r[0], r[1]) for r in results_matched], bounds)

                # if coordinates can't be merged pick largest span
                if len(combined_coords) > 1:
                    combined_coords.sort(key=lambda c:c[1]-c[0], reverse=True)

                if combined_coords:
                    #covered = True
                    check_seq_len = abs(len(repeat_seqs[seq]) - 2 * flank)
                    span = float(combined_coords[0][1] - combined_coords[0][0] + 1)
                    #if check_seq_len < 300:
                        #min_span = 0.3
                    #covered = (span / check_seq_len) >= min_span

                    if check_seq_len == 0 or (span / check_seq_len) < min_span:
                        continue
                    #if abs(len(repeat_seqs[seq]) - 2 * flank) > 0 and span >= 100:
                        #covered = (span / abs(len(repeat_seqs[seq]) - 2 * flank)) >= min_span
                    #if not covered:
                        #continue

                #if combined_coords and\
                   #(abs(len(repeat_seqs[seq]) - 2 * flank) > 0 and\
                    #(combined_coords[0][1] - combined_coords[0][0] + 1) >= 100 and\
                    #float(combined_coords[0][1] - combined_coords[0][0] + 1) / abs(len(repeat_seqs[seq]) - 2 * flank) > min_span):
                    coords = combined_coords[0]
                    repeat_seq = repeat_seqs[seq][coords[0]-1:coords[-1]]
                    size = coords[-1] - coords[0] + 1
                    rpos = rstart + coords[0] - 1
                    pats = set([r[-2] for r in results_matched])
                    genome_start = int(gstart) + coords[0]
                    genome_end = int(gend) - (seq_len - coords[-1])

                    print('ff {} {} {} {} {} {} {} {} {} {}'.format(read, locus, size, rpos, gstart, gend, seq_len, coords, genome_start, genome_end))
                    if not read in alleles[locus]:
                        if genome_start < genome_end:
                            alleles[tuple(locus)][read] = (rpos, pats, size, genome_start, genome_end)
                        else:
                            alleles[tuple(locus)][read] = (rpos, pats, size, int(gstart), int(gend))

        return self.alleles_to_variants(alleles)

    def get_read_seqs(self, headers, bam):
        grouped_reads = defaultdict(list)
        for header in headers:
            chrom, start, end, read = header.split(':')[:4]
            grouped_reads[(chrom, start, end)].append(read)

        read_seqs = {}
        for region, reads in grouped_reads.iteritems():
            for aln in bam.fetch(region[0], int(region[1]), int(region[2])):
                if aln.query_name in reads:
                    read_seqs[aln.query_name] = aln.query_sequence

        return read_seqs

    def extract_alleles_regex(self, headers, repeat_seqs, flank, clipped, bam, strands, patterns):
        alleles = defaultdict(dict)

        read_seqs = self.get_read_seqs(headers, bam)

        for header in headers:
            cols = header.split(':')
            read = cols[3]
            if not header in repeat_seqs.keys():
                continue
            if not read in read_seqs:
                print('cannot get sequence:{}'.format(read))
                continue
            repeat_seq = repeat_seqs[header][flank:-1*flank]
            repeat = patterns[header]
            has_repeat, rspan, pats = self.examine_repeats(repeat_seq, repeat)
            matched_seq = repeat_seq[rspan[0] : rspan[1] + 1]
            if has_repeat and pats:
                locus = tuple(cols[:3])
                read_seq = read_seqs[read]
                rpos = read_seq.find(matched_seq)
                size = len(matched_seq)
                print('ff {} {} {} {}'.format(read, locus, size, rpos))

                if not read in alleles[locus]:
                    alleles[tuple(locus)][read] = (rpos, pats, size)

        return self.alleles_to_variants(alleles)

    def alleles_to_variants(self, alleles):
        variants = []
        for locus in alleles.keys():
            variant = [locus[0], locus[1], locus[2], [], None, [], '-']

            # check for minimum number of supporting reads
            if len(alleles[locus]) < self.min_support:
                continue

            pat_counts = Counter()
            for read in alleles[locus]:
                variant[3].append([read,
                                   alleles[locus][read][0], # rstart
                                   None, # repeat
                                   None, # copy number
                                   alleles[locus][read][2], # size
                                   alleles[locus][read][3], # genome_start
                                   alleles[locus][read][4], # genome_end
                                   ])

                # update pattern counts
                pat_counts.update(alleles[locus][read][1])

            pat_counts_sorted = pat_counts.most_common()
            top_pats = [pat_counts_sorted[0][0]]
            for i in range(1, len(pat_counts_sorted)):
                if pat_counts_sorted[i][1] == pat_counts_sorted[0][1]:
                    top_pats.append(pat_counts_sorted[i][0])
                else:
                    break
            variant[4] = (sorted(top_pats, key=len)[0])
            for allele in variant[3]:
                allele[2] = variant[4]
                allele[3] = round(float(allele[4]) / len(variant[4]), 1)

            variants.append(variant)

        return variants

    def get_alleles(self, loci, reads_fasta=None, closeness_to_end=200):
        bam = pysam.Samfile(self.bam, 'rb')
        if self.reads_fasta:
            reads_fasta = pysam.Fastafile(self.reads_fasta)

        trf_input = ''
        strands = {}
        repeat_seqs = {}
        generic = set()
        patterns = {}
        single_neighbour_size = 100
        split_neighbour_size = 500
        min_mapped = 0.5

        if self.strict:
            self.trf_flank_size = 50

        for locus in loci:
            clipped = defaultdict(dict)
            alns = []
            for aln in bam.fetch(locus[0], locus[1] - split_neighbour_size, locus[2] + split_neighbour_size):
                alns.append(aln)
                locus_size = locus[2] - locus[1] + 1

                # check split alignments first
                if self.check_split_alignments and\
                   ((aln.reference_start >= locus[1] - split_neighbour_size and aln.reference_start <= locus[2] + split_neighbour_size) or\
                    (aln.reference_end >= locus[1] - split_neighbour_size and aln.reference_end <= locus[2] + split_neighbour_size)):
                    check_end = None
                    start_olap = False
                    end_olap = False
                    if aln.reference_start >= locus[1] - split_neighbour_size and aln.reference_start <= locus[2] + split_neighbour_size:
                        start_olap = True
                    if aln.reference_end >= locus[1] - split_neighbour_size and aln.reference_end <= locus[2] + split_neighbour_size:
                        end_olap = True
                    if start_olap and not end_olap:
                        check_end = 'start'
                    elif end_olap and not start_olap:
                        check_end = 'end'
                    clipped_end, partner_start = INSFinder.is_split_aln_potential_ins(aln, min_split_size=400, closeness_to_end=10000, check_end=check_end)
                    if clipped_end is not None:
                        clipped[aln.query_name][clipped_end] = (aln, partner_start)

            # clipped alignment
            remove = set()
            for read in clipped.keys():
                if len(clipped[read].keys()) == 2:
                    aln1 = clipped[read]['end'][0]
                    aln2 = clipped[read]['start'][0]
                    if reads_fasta or aln1.query_alignment_end < aln2.query_alignment_start:
                        aln1_tuple = self.extract_aln_tuple(aln1, locus[1] - self.trf_flank_size, 'left')
                        aln2_tuple = self.extract_aln_tuple(aln2, locus[2] + self.trf_flank_size, 'right')
                        qstart = None
                        qend = None
                        tstart = None
                        tend = None
                        if aln1_tuple and aln2_tuple:
                            qstart, tstart = aln1_tuple
                            qend, tend = aln2_tuple
                            if not reads_fasta:
                                seq = aln1.query_sequence[qstart:qend]
                            else:
                                if aln1.cigartuples[0][0] == 5:
                                    qstart = aln1.cigartuples[0][1] + qstart
                                if aln2.cigartuples[0][0] == 5:
                                    qend = aln2.cigartuples[0][1] + qend
                                seq = INSFinder.get_seq(reads_fasta, aln1.query_name, aln1.is_reverse, [qstart, qend])
                            alns.remove(aln1)
                            alns.remove(aln2)
                        else:
                            print('problem getting seq2 {} {}'.format(aln.query_name, locus))
                            continue

                        # leave patterns out, some too long for trf header
                        header, fa_entry = self.create_trf_fasta(locus[:3], aln1.query_name, tstart, tend, qstart, seq, aln1.infer_read_length())
                        patterns[header] = locus[-1]
                        repeat_seqs[header] = seq

                        if not '*' in locus[-1]:
                            trf_input += fa_entry
                        else:
                            generic.add(header)
                            repeat_seqs[header] = seq

                        strands[read] = aln1.is_reverse

                    else:
                        remove.add(read)
                else:
                    clipped_end = list(clipped[read].keys())[0]
                    aln = clipped[read][clipped_end][0]
                    alns.remove(aln)

            for read in remove:
                del clipped[read]

            for aln in alns:
                if aln.reference_start <= locus[1] - single_neighbour_size and\
                   aln.reference_end >= locus[2] + single_neighbour_size:
                    seq, tstart, tend, qstart = self.extract_subseq(aln, locus[1] - self.trf_flank_size, locus[2] + self.trf_flank_size, reads_fasta=reads_fasta)
                    if seq is None:
                        print('problem getting seq1 {} {} {} {} {}'.format(aln.query_name, locus, tstart, tend, qstart))
                        continue
                    # leave patterns out, some too long for trf header
                    header, fa_entry = self.create_trf_fasta(locus[:3], aln.query_name, tstart, tend, qstart, seq, aln.infer_read_length())
                    patterns[header] = locus[-1]
                    repeat_seqs[header] = seq

                    if not '*' in locus[-1]:
                        trf_input += fa_entry
                    else:
                        generic.add(header)

                    # for finding position of repeat sequence inside read
                    strands[aln.query_name] = aln.is_reverse

        variants = []
        if trf_input:
            variants.extend(self.extract_alleles_trf(trf_input,
                                                     repeat_seqs,
                                                     self.trf_flank_size,
                                                     clipped,
                                                     bam,
                                                     strands,
                                                     patterns))

        if generic:
            variants.extend(self.extract_alleles_regex(generic,
                                                       repeat_seqs,
                                                       self.trf_flank_size,
                                                       clipped,
                                                       bam,
                                                       strands,
                                                       patterns))

        # set genotyping configuration (class variable)
        Variant.set_genotype_config(method=self.clustering,
                                    min_reads=self.min_support,
                                    max_num_clusters=self.max_num_clusters,
                                    eps=self.eps)

        for variant in variants:
            # genotype
            Variant.genotype(variant, report_in_size=self.genotype_in_size)
            Variant.summarize_genotype(variant)

        if self.remove_tmps:
            self.cleanup()

        return variants

    def examine_ins(self, ins_list, min_expansion=0):
        if self.nprocs > 1:
            print('gg {}'.format(len(ins_list)))
            random.shuffle(ins_list)
            batches = list(split_tasks(ins_list, self.nprocs))
            batched_results = parallel_process(self.extract_tres, batches, self.nprocs)
            expansions = combine_batch_results(batched_results, type(batched_results[0]))
        else:
            expansions = self.extract_tres(ins_list)

        # label ins
        self.annotate(ins_list, expansions)

        tre_events = [ins for ins in ins_list if ins[6] == 'tre' and\
                      len(ins[8]) >= self.min_str_len and\
                      len(ins[8]) <= self.max_str_len]

        if tre_events:
            merged_loci = self.merge_loci(tre_events)
            print('mm {}'.format(len(merged_loci)))
            self.make_trf_sensitive()
            variants = self.collect_alleles(merged_loci)

            # remove redundant variants
            if variants:
                self.remove_redundants(variants)

                # udpate coordinates
                for variant in variants:
                    Variant.update_coords(variant)

                return [v for v in variants if Variant.above_min_expansion(v, min_expansion)]

        return []

    def remove_redundants(self, variants):
        groups = defaultdict(list)
        for i in range(len(variants)):
            variant = variants[i]
            key = '-'.join(map(str, [variant[0], variant[1], variant[2]]))
            groups[key].append(i)

        remove = []
        for key in groups.keys():
            if len(groups[key]) > 1:
                remove.extend(sorted(groups[key])[1:])

        for i in sorted(remove, reverse=True):
            del variants[i]

    def make_trf_sensitive(self):
        self.trf_args = '2 5 5 80 10 10 500 -d -h'

    def collect_alleles(self, loci):
        tre_variants = []
        if self.nprocs > 1:
            random.shuffle(loci)
            batches = list(split_tasks(loci, self.nprocs))
            batched_results = parallel_process(self.get_alleles, batches, self.nprocs)
            tre_variants = combine_batch_results(batched_results, type(batched_results[0]))
        else:
            tre_variants = self.get_alleles(loci)

        return tre_variants
    
    def genotype(self, loci_bed):
        # extract loci (return fake insertions)
        loci = []
        with open(loci_bed, 'r') as ff:
            for line in ff:
                cols = line.rstrip().split()
                if len(cols) == 4:
                    if len(cols[3]) <= self.max_str_len:
                        loci.append((cols[0], int(cols[1]), int(cols[2]), cols[3]))

        # use more senstitive trf params to capture smaller alleles
        self.make_trf_sensitive()

        # use give loci coordinates for reporting
        self.update_loci = False

        self.strict = True

        return self.collect_alleles(loci)
    
    def output(self, variants, out_file):
        with open(out_file, 'w') as out:
            out.write('#{}\n'.format('\t'.join(Variant.tsv_headers + Allele.tsv_headers)))
            for variant in sorted(variants, key=itemgetter(0, 1, 2)):
                if not variant[5]:
                    continue
                variant_cols = Variant.to_tsv(variant)
                for allele in sorted(variant[3], key=itemgetter(3), reverse=True):
                    allele_cols = Allele.to_tsv(allele)
                    out.write('{}\n'.format('\t'.join(variant_cols + allele_cols)))

    def cleanup(self):
        if self.tmp_files:
            for ff in self.tmp_files:
                if os.path.exists(ff):
                    os.remove(ff)
