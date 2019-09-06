import pysam
from distutils import spawn
import sys
import subprocess
import os
import re
from variant import Variant, Allele
from collections import defaultdict
from operator import itemgetter
import itertools
from utils import split_tasks, parallel_process, combine_batch_results, create_tmp_file, reverse_complement

class ExpandedRepeatFinder:
    def __init__(self, bam, genome_fasta):
        self.bam = bam
        self.genome_fasta = genome_fasta
        trf_path = spawn.find_executable("trf")
        if not trf_path:
            sys.exit('ABORT: %s' % "can't find trf in PATH")
            
        self.trf_args = '2 5 5 80 10 50 500 -d -h'
        self.flank_len = 50

        # for checking sequences flanking repeat
        self.trf_flank_min_mapped_fraction = 0.9
        self.trf_flank_size = 1000

    def construct_trf_output(self, input_fasta):
        m = re.search('(\d[\d\s]*\d)', self.trf_args)
        if m is not None:
            return '%s/%s.%s.dat' % (os.getcwd(), os.path.basename(input_fasta), m.group(1).replace(' ', '.'))
    
    def run_trf(self, input_fasta):
        cmd = ' '.join(['trf', input_fasta, self.trf_args])
        # redirect stdout and stderr to devnull
        FNULL = open(os.devnull, 'w')
        returncode = subprocess.call(cmd, shell=True, stdout=FNULL, stderr=FNULL)
        
        output = self.construct_trf_output(input_fasta)
        if os.path.exists(output):
            return output
        else:
            sys.exit('cannot run %s' % cmd)

    #def chain_trf(self, results):
        #all_to_chain = []
        #skip = set()
        #for i in range(len(results) - 1):
            #if i in skip:
                #continue
            #to_chain = set()
            #j = i + 1
            #while j < len(results):
                #if not j in skip and \
                   #len(results[i][13]) == len(results[j][13]) and \
                   #self.is_same_repeat((results[i][13], results[j][13])):
                    #to_chain.add(i)
                    #to_chain.add(j)
                #j += 1
            #if to_chain:
                #skip |= to_chain
                #all_to_chain.append(to_chain)

        #chained_results = []
        #for indexes in all_to_chain:
            #chained_results.append(self.combine_trf([results[i] for i in indexes]))
        #for i in range(len(results)):
            #if not i in skip:
                #chained_results.append(results[i])

        #return chained_results

    #def merge_trf(self, results):
        #if len(results) == 1:
            #return results

        #results_merged = []
        #i = 0
        #j = None
        #while i < len(results) - 1:
            #j = i + 1
            #to_combine = [results[i]]
            #pat = results[i][13]
            #while j < len(results):
                #if len(results[i][13]) == len(results[j][13]) and self.is_same_repeat((results[i][13], results[j][13])):
                    #to_combine.append(results[j])
                    #j += 1
                #else:
                    #break

            ## combine
            #combined = self.combine_trf(to_combine)
            #results_merged.append(combined)

            #i = j

        #if j is not None and j == len(results) - 1:
            #results_merged.append(self.type_trf_cols(results[j]))

        #return results_merged

    #def combine_trf(self, results):
        #"""
            #0 = start index
            #1 = end index
            #2 = period size
            #3 = copy number
            #4 = size of consensus pattern
            #5 = percent match between adjacent copies
            #6 = percent indels between adjacent copies
            #7 = alignment score
            #8 = percentage composition A
            #9 = percentage composition C
            #10 = percentage composition G
            #11 = percentage composition T
            #12 = entropy
            #13 = pattern sequence
            #14 = repeat sequence
        #"""
        #if len(results) == 1:
            #return results[0]

        #combined = []
        #for i in range(15):
            #if i == 0:
                #col = min([r[i] for r in results])
            #elif i == 1 or (i >= 5 and i <= 7) or i == 12:
                #col = max([r[i] for r in results])
            #elif i == 2 or i == 4 or i == 13:
                #col = results[0][i]
            #elif i == 3 or (i >= 8 and i <=11):
                #col = sum([r[i] for r in results])
            #elif i == 14:
                #col = ','.join([r[i] for r in results])
            #combined.append(col)

        #return combined

    def type_trf_cols(self, cols):
        return map(int, cols[:3]) + [float(cols[3])] + map(int, cols[4:12]) + [float(cols[12])] + cols[13:]
    
    def parse_trf(self, trf_output, merge=False):
        """
            columns:
            start_index, end_index, periold_size, copy_number, consensus_size,
            percent_matches, percent_index, score, A, C, G, T, entropy,
            consensus_pattern

            start and end index are 1-based
        """
        results = defaultdict(list)
        with open(trf_output, 'r') as ff:
            for line in ff:
                cols = line.rstrip().split()
                if not cols:
                    continue
                if cols[0] == 'Sequence:':
                    seq = cols[1]
                elif len(cols) == 15:
                    results[seq].append(self.type_trf_cols(cols))
                    
        if not merge:
            return results

        results_merged = {}
        for seq in results.keys():
            #results_merged[seq] = self.merge_trf(results[seq])
            results_merged[seq] = self.chain_trf(results[seq])

        return results_merged

    def extract_tres(self, ins_list, target_flank=2000):
        genome_fasta = pysam.Fastafile(self.genome_fasta)

        # prepare input for trf
        trf_input = ''
        for ins in ins_list:
            target_start, target_end = ins.gpos - target_flank + 1, ins.gpos + target_flank
            prefix = '.'.join(map(str, [ins.eid, len(ins.seq), target_start, target_end]))
            trf_input += '>%s.q\n%s\n' % (prefix, ins.neighbour_seq)
            #trf_input += '>%s.t\n%s\n' % (prefix, self.extract_genome_neighbour(ins.chrom,
                                                                                #ins.gpos,
                                                                                #len(ins.neighbour_seq)/2,
                                                                                #genome_fasta))
            trf_input += '>%s.t\n%s\n' % (prefix, self.extract_genome_neighbour(ins.chrom,
                                                                                ins.gpos,
                                                                                target_flank,
                                                                                genome_fasta))
            trf_input += '>%s.i\n%s\n' % (prefix, ins.seq)

        trf_fasta = create_tmp_file(trf_input)

        output = self.run_trf(trf_fasta)

        results = self.parse_trf(output)

        grouped_results = {}
        for seq in results.keys():
            read, read_type = seq.rsplit('.', 1)
            if not grouped_results.has_key(read):
                grouped_results[read] = {}
            grouped_results[read][read_type] = results[seq]

        expansions = self.analyze_trf(grouped_results)

        return expansions

    def extract_genome_neighbour(self, chrom, tpos, w, genome_fasta):
        return genome_fasta.fetch(chrom, tpos - w, tpos + w)
    
    def analyze_trf(self, results, full_cov=0.8):
        expansions = {}
        for seq_id in sorted(results.keys()):
            #eid, ins_len = seq_id.rsplit('.', 1)
            eid, ins_len, gstart, gend = seq_id.rsplit('.', 3)
            
            if (results[seq_id].has_key('i') or results[seq_id].has_key('q')) and results[seq_id].has_key('t'):
                pat, pgstart, pgend = self.analyze_trf_per_seq(results[seq_id], int(ins_len), int(gstart), int(gend))
                if pat is not None:
                    expansions[eid] = pat, pgstart, pgend

        return expansions
    
    def is_same_repeat(self, reps, min_fraction=0.6):
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
                if float(rep2.count(p1) * len(p1)) / len(rep2) > min_fraction:
                    return True

        return False

    def is_same_repeat_old(self, reps):
        if len(reps[0]) <= len(reps[1]):
            rep1, rep2 = reps[0], reps[1]
        else:
            rep2, rep1 = reps[0], reps[1]
        perms2 = []
        for i in range(len(rep2)):
            perms2.append(rep2[i:] + rep2[:i])

        if rep1 in perms2 or [p2 for p2 in perms2 if rep1 in p2]:
            return True
        else:
            return False

    def analyze_trf_per_seq(self, result, ins_len, gstart, gend, full_cov=0.8):
        pgstart, pgend = None, None

        filtered_patterns = {'i':[], 'q':[], 't':[]}
        if result.has_key('i'):
            filtered_patterns['i'] = [r[13] for r in result['i'] if float(r[1] - r[0] + 1) / ins_len >= full_cov]

        mid_pts = self.flank_len, self.flank_len + 1
        for pt in ('q', 't'):
            if result.has_key(pt):
                filtered_patterns[pt] = [r[13] for r in result[pt] if r[0] <= mid_pts[0] or r[1] >= mid_pts[1]]

        pattern_matched = None
        for i_pat in sorted(filtered_patterns['i'], key=len):
            if len(i_pat) == 1:
                continue
            i_pat_matches = {'q': False, 't': False}
            for pt in ('t', 'q'):
                for pat in sorted(filtered_patterns[pt], key=len, reverse=True):
                    if self.is_same_repeat((i_pat, pat)):
                        i_pat_matches[pt] = True
                        break

            if i_pat_matches['q'] or i_pat_matches['t']:
                pattern_matched = i_pat

                if i_pat_matches['t']:
                    for r in result['t']:
                        if self.is_same_repeat((r[13], i_pat)):
                            pgstart, pgend = gstart + r[0] - 1, gstart + r[1] - 1

                break

        return pattern_matched, pgstart, pgend

    def analyze_trf_per_seq_old(self, result, ins_len, full_cov=0.8):
        filtered_patterns = {'i':[], 'q':[], 't':[]}
        filtered_patterns['i'] = [p for p in result['i'].keys() if float(result['i'][p][1] - result['i'][p][0] + 1) / ins_len >= full_cov]
        
        mid_pts = self.flank_len, self.flank_len + 1
        for pt in ('q', 't'):
            #mid_pts = lengths[pt] / 2, -(-1 * lengths[pt] / 2)
            filtered_patterns[pt] = [p for p in result[pt].keys() if result[pt][p][0] <= mid_pts[0] or result[pt][p][1] >= mid_pts[1]]
            
        #print 'kk', filtered_patterns
        pattern_matched = None
        for i_pat in sorted(filtered_patterns['i'], key=len):
            if len(i_pat) == 1:
                continue
            i_pat_matches = {'q': False, 't': False}
            for pt in ('q', 't'):
                for pat in sorted(filtered_patterns[pt], key=len, reverse=True):
                    if self.is_same_repeat((i_pat, pat)):
                        i_pat_matches[pt] = True
                        break

            if i_pat_matches['q'] or i_pat_matches['t']:
                pattern_matched = i_pat
                break
            
        return pattern_matched

    def annotate(self, ins_list, expansions):
        ins_dict = dict((ins.eid, ins) for ins in ins_list)
        for eid in expansions.keys():
            if ins_dict.has_key(eid):
                ins_dict[eid].etype = 'tre'
                ins_dict[eid].repeat_unit = expansions[eid][0]
                if expansions[eid][1] is not None and expansions[eid][2] is not None:
                    ins_dict[eid].gpos = expansions[eid][1]
                    ins_dict[eid].gend = expansions[eid][2]

    def merge_loci(self, tres, w=100, s=0.2):
        #loci_sorted = sorted([[tre.chrom, tre.gpos, tre.gpos+1, tre.repeat_unit] for tre in tres], key=itemgetter(0, 1))
        loci_sorted = sorted([[tre.chrom, tre.gpos, tre.gend, tre.repeat_unit, len(tre.repeat_unit)] for tre in tres], key=itemgetter(0, 1, 4))

        # list of new tuples representing merges
        merges = []
        # mapping of index of merges to index of sorted_tuples
        merge_indices = defaultdict(set)
        # list of singleton tuples
        singles = [loci_sorted[0]]
        # corresponding indices of singleton tuple
        single_indices = [0]

        for i in range(1, len(loci_sorted)):
            locus = loci_sorted[i]
            merged = False

            j = len(merges) - 1
            while j >= 0:
                before = merges[j]
                if before[0] == locus[0] and before[2] + w >= locus[1]:
                    if self.is_same_repeat((before[3], locus[3])):
                        #before = (before[0], before[1], locus[2], min(before[3], locus[3], key=len))
                        before[2] = max(before[2], locus[2])
                        before[3] = min(before[3], locus[3], key=len)
                        merge_indices[j].add(i)
                        merged = True
                        break
                    j -= 1
                else:
                    break

            if not merged and singles:
                before = singles[-1]
                if before[0] == locus[0] and before[2] + w >= locus[1] and self.is_same_repeat((before[3], locus[3])):
                    new_locus = [before[0], before[1], locus[2], min(before[3], locus[3], key=len)]
                    merges.append(new_locus)
                    merge_indices[len(merges) - 1].add(single_indices[-1])
                    merge_indices[len(merges) - 1].add(i)

                    singles.pop()
                    single_indices.pop()
                    merged = True

            if not merged:
                singles.append(locus)
                single_indices.append(i)

        return merges

    #def merge_loci_old(self, tres, w=100, s=0.2):
        #loci_sorted = sorted([(tre.chrom[0], tre.gpos[0], tre.gpos[1], tre.repeat_unit) for tre in tres], key=itemgetter(1))

        #merged = []
        #for i in range(len(loci_sorted)):
            #locus = loci_sorted[i]
            #print 'zz', locus
            #if not merged:
                #merged.append(locus)
            #else:
                #before = merged.pop()
                #if before[2] + w >= locus[1] and self.is_same_repeat((before[3], locus[3])):
                    #new_locus = (before[0], before[1], locus[2], min(before[3], locus[3], key=len))
                    #merged.append(new_locus)
                #else:
                    #merged.append(before)
                    #merged.append(locus)

        #for locus in merged:
            #print 'll', locus

        #return [merged[0]]

    def extract_subseq(self, aln, target_start, target_end, max_extend=50):
        aligned_pairs = aln.get_aligned_pairs()
        query_start = None
        query_end = None
        # actual target start and end used
        tstart, tend = None, None

        # start
        found = []
        tpos = target_start
        while (not found or found[0] is None) and tpos >= target_start - max_extend:
            found = [p[0] for p in aligned_pairs if p[1] == tpos]
            tstart = tpos
            tpos -= 1
        if found:
            query_start = found[0]

        # end
        found = []
        tpos = target_end
        while (not found or found[0] is None) and tpos <= target_end + max_extend:
            found = [p[0] for p in aligned_pairs if p[1] == tpos]
            tend = tpos
            tpos += 1
        if found:
            query_end = found[0]

        if query_start and query_end:
            return aln.query_sequence[query_start:query_end], tstart, tend

        return None, tstart, tend
    
    def check_trf_prediction_fullness(self, start, end, aligned_pairs, read=None):
        left_flank = sorted([p[1] for p in aligned_pairs if (p[0] >= start - self.trf_flank_size and p[0] < start and p[1] is not None)])
        right_flank = sorted([p[1] for p in aligned_pairs if (p[0] <= end + self.trf_flank_size and p[0] > end and p[1] is not None)])

        chrom_start = None
        chrom_end = None
        if left_flank and left_flank[0] is not None and right_flank and right_flank[0] is not None:
            chrom_start = max(left_flank)
            chrom_end = min(right_flank)

        #if read:
            #print 'ff', read, start, end, len(left_flank), float(len(left_flank)) / self.trf_flank_size, len(right_flank), float(len(right_flank)) / self.trf_flank_size
        if float(len(left_flank)) / self.trf_flank_size >= self.trf_flank_min_mapped_fraction and \
           float(len(right_flank)) / self.trf_flank_size >= self.trf_flank_min_mapped_fraction:
            return chrom_start, chrom_end
        else:
            return False

    #def sort_trf_results(self, results, pat):
        #results_matched = [r for r in results if self.is_same_repeat((r[13], pat))]
        #counts = []
        #for result in results_matched:
            #counts.append(max(result[-1].count(pat), result[-1].count(reverse_complement(pat))))

        #sorted_indices = sorted(range(len(counts)), key = lambda i: counts[i], reverse=True)
        #results_matched_sorted = [results_matched[i] for i in sorted_indices]
        #counts_sorted = [counts[i] for i in sorted_indices]

        #return results_matched_sorted, counts_sorted

    def get_alleles(self, loci, flank=50):
        def get_distance(repeat_start, repeat_end, pat_start, pat_end):
            d = 0
            # pat subsumed in repeat
            if pat_start <= repeat_start and pat_end >= repeat_end:
                d = -1 * (repeat_end - repeat_start + 1)
            # no overlap on left
            elif pat_end < repeat_start:
                d = repeat_start - pat_end
            # no overlap on right
            elif pat_start > repeat_end:
                d = pat_start - repeat_end
            # overlap on left
            elif pat_end >= repeat_start:
                d = -1 * (pat_end - repeat_start)
            # overlap on right
            elif pat_start <= repeat_end:
                d = -1 * (repeat_end - pat_start)

            return d

        bam = pysam.Samfile(self.bam, 'rb')

        trf_input = ''
        read_seqs = {}
        repeat_seqs = {}
        aligned_pairs = {}
        for locus in loci:
            for aln in bam.fetch(locus[0], locus[1], locus[2]):
                #if aln.query != 'b8d719b4-0f83-4ff7-8108-4451e58d2b13':
                    #continue
                if aln.reference_start <= locus[1] - flank and\
                   aln.reference_end >= locus[2] + flank:
                    seq, tstart, tend = self.extract_subseq(aln, locus[1] - flank, locus[2] + flank)
                    if seq is None:
                        print 'problem getting seq', aln.query_name, locus, len(aln.query_sequence)
                        continue
                    header = '>%s' % ':'.join(map(str, [locus[0],
                                                        locus[1],
                                                        locus[2],
                                                        locus[3],
                                                        aln.query_name,
                                                        tstart,
                                                        tend,
                                                        len(seq)]))
                    if seq is not None:
                        trf_input += '%s\n%s\n' % (header, seq)
                    # for finding position of repeat sequence inside read
                    read_seqs[aln.query_name] = aln.query_sequence
                    repeat_seqs[aln.query_name] = seq
                    aligned_pairs[aln.query_name] = aln.get_aligned_pairs()

        trf_fasta = create_tmp_file(trf_input)
        print trf_fasta
        output = self.run_trf(trf_fasta)

        results = self.parse_trf(output, merge=False)

        # group by locus
        alleles = defaultdict(dict)
        for seq in results.keys():
            cols = seq.split(':')
            locus = tuple(cols[:4])
            expected_pat = cols[3]
            read = cols[4]
            gstart, gend = int(cols[5]), int(cols[6])

            results_matched = [result for result in results[seq] if self.is_same_repeat((expected_pat, result[13]))]
            #for result in results[seq]:
                #print 'hh2', seq, result, self.tmp_same_repeat((expected_pat, result[13]))
            #for result in results_matched:
                #print 'hh', seq, result
            if results_matched:
                coords = sorted(list(itertools.chain.from_iterable([[int(result[0]), int(result[1])] for result in results_matched])))
                repeat_seq = repeat_seqs[read][coords[0]-1:coords[-1]]
                size = len(repeat_seq)
                cn = round(float(len(repeat_seq)) / len(expected_pat), 1)
                rpos = read_seqs[read].find(repeat_seq)
                rstart = rpos
                rend = rstart + len(repeat_seq)

                flanks_aligned = self.check_trf_prediction_fullness(rstart, rend, aligned_pairs[read], read=read)
                #print 'hhh', seq, coords, coords[-1] - coords[0], float((coords[-1] - coords[0] + 1))/len(expected_pat), size, cn, repeat_seq, flanks_aligned
                if flanks_aligned:
                    d  = get_distance(int(locus[1]), int(locus[2]), flanks_aligned[0], flanks_aligned[1])
                    if d < 20:
                        if not alleles[locus].has_key(read) or \
                           len(alleles[locus][read][3]) < len(result[-1]):
                            alleles[locus][read] = (rpos, expected_pat, cn, repeat_seq, flanks_aligned)

        variants = []
        for locus in alleles.keys():
            variant = Variant(locus[0], locus[1], locus[2], locus[3])
            for read in alleles[locus]:
                variant.alleles.append(Allele(read,
                                              alleles[locus][read][0],
                                              alleles[locus][read][1],
                                              alleles[locus][read][2],
                                              alleles[locus][read][3],
                                              alleles[locus][read][4][0],
                                              alleles[locus][read][4][1]))
            variants.append(variant)

        for variant in variants:
            # update variant coordinates
            variant.update_coords()

            # genotype
            variant.genotype()

        return variants

    def examine_ins(self, ins_list, nprocs=1):
        if nprocs > 1:
            batches = list(split_tasks(ins_list, nprocs))
            batched_results = parallel_process(self.extract_tres, batches, nprocs)
            expansions = combine_batch_results(batched_results, type(batched_results[0]))
        else:
            expansions = self.extract_tres(ins_list)

        # label ins
        self.annotate(ins_list, expansions)

        # triage tre and non-tre events
        tre_events = []
        non_tre_events = []
        for ins in ins_list:
            if ins.etype == 'tre':
                tre_events.append(ins)
            else:
                non_tre_events.append(ins)

        merged_loci = self.merge_loci(tre_events)

        tre_variants = self.collect_alleles(merged_loci, nprocs)

        # remove redundant variants
        self.remove_redundants(tre_variants)

        return tre_variants

    def remove_redundants(self, variants):
        groups = defaultdict(list)
        for i in range(len(variants)):
            variant = variants[i]
            key = '-'.join(map(str, [variant.chrom, variant.start, variant.end]))
            groups[key].append(i)

        remove = []
        for key in groups.keys():
            if len(groups[key]) > 1:
                remove.extend(sorted(groups[key])[1:])

        for i in sorted(remove, reverse=True):
            del variants[i]

    def make_trf_sensitive(self):
        self.trf_args = '2 5 5 80 10 10 500 -d -h'

    def collect_alleles(self, loci, nprocs):
        tre_variants = []
        if nprocs > 1:
            batches = list(split_tasks(loci, nprocs))
            batched_results = parallel_process(self.get_alleles, batches, nprocs)
            tre_variants = combine_batch_results(batched_results, type(batched_results[0]))
        else:
            tre_variants = self.get_alleles(loci)

        return tre_variants