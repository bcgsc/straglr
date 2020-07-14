import pysam
from distutils import spawn
import sys
import subprocess
import os
import re
from intspan import intspan
from variant import Variant, Allele
from collections import defaultdict
from operator import itemgetter, attrgetter
import itertools
from utils import split_tasks, parallel_process, combine_batch_results, create_tmp_file, reverse_complement
from ins import INSFinder

class TREFinder:
    def __init__(self, bam, genome_fasta, check_split_alignments=False, max_str_len=50, flank_size=100, nprocs=1):
        self.bam = bam
        self.genome_fasta = genome_fasta
        trf_path = spawn.find_executable("trf")
        if not trf_path:
            sys.exit('ABORT: %s' % "can't find trf in PATH")
            
        self.trf_args = '2 5 5 80 10 50 500 -d -h'
        #self.flank_len = 100
        self.flank_len = flank_size

        # for checking sequences flanking repeat
        self.trf_flank_min_mapped_fraction = 0.8
        #self.trf_flank_size = 100
        self.trf_flank_size = flank_size
        self.nprocs = nprocs

        self.check_split_alignments = check_split_alignments

        self.update_loci = True

        self.max_str_len = max_str_len

        self.tmp_files = set()

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

    def type_trf_cols(self, cols):
        return map(int, cols[:3]) + [float(cols[3])] + map(int, cols[4:12]) + [float(cols[12])] + cols[13:]
    
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
                    results[seq].append(self.type_trf_cols(cols))

        return results

    def extract_tres(self, ins_list, target_flank=2000):
        genome_fasta = pysam.Fastafile(self.genome_fasta)

        # prepare input for trf
        trf_input = ''
        for ins in ins_list:
            target_start, target_end = ins.gpos - target_flank + 1, ins.gpos + target_flank
            prefix = '.'.join(map(str, [ins.eid, len(ins.seq), target_start, target_end]))
            trf_input += '>%s.q\n%s\n' % (prefix, ins.neighbour_seq)
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

    def combine_trf_coords(self, coords, max_sep=20):
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

    def analyze_trf_per_seq(self, result, ins_len, gstart, gend, full_cov=0.8):
        pgstart, pgend = None, None

        filtered_patterns = {'i':[], 'q':[], 't':[]}
        if result.has_key('i'):
            filtered_patterns['i'] = [r[13] for r in result['i'] if float(r[1] - r[0] + 1) / ins_len >= full_cov]

        # make sure detected repeat pattern crosses mid-point, i.e. the intended location
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
        """
            minimum flank size = 100, one of them has to minimum 100,
            each side at least 80% mapped
        """
        flank_start = max(0, start - self.trf_flank_size)
        if len(aligned_pairs) == 2:
            flank_end = min(end + self.trf_flank_size, aligned_pairs[0][-1][0], aligned_pairs[1][-1][0])
        else:
            flank_end = min(end + self.trf_flank_size, aligned_pairs[-1][0])
        left_flank_size = start - flank_start
        right_flank_size = flank_end - end

        pairs = aligned_pairs
        if len(aligned_pairs) == 2:
            pairs = aligned_pairs[0]
        left_flank = sorted([p[1] for p in pairs if (p[0] >= flank_start and p[0] < start and p[1] is not None)])

        if len(aligned_pairs) == 2:
            pairs = aligned_pairs[1]
        right_flank = sorted([p[1] for p in pairs if (p[0] <= flank_end and p[0] > end and p[1] is not None)])

        chrom_start = None
        chrom_end = None
        if left_flank and left_flank[0] is not None and right_flank and right_flank[0] is not None:
            chrom_start = max(left_flank)
            chrom_end = min(right_flank)
            
        if left_flank_size > 0 and right_flank_size > 0 and\
           (left_flank_size >= self.trf_flank_size or right_flank_size >= self.trf_flank_size) and\
           float(len(left_flank)) / left_flank_size >= self.trf_flank_min_mapped_fraction and \
           float(len(right_flank)) / right_flank_size >= self.trf_flank_min_mapped_fraction:
            return sorted([chrom_start, chrom_end])
        else:
            return False

    def create_trf_fasta(self, locus, read, tstart, tend, seq):
        """ for genotyping """
        header = '%s' % ':'.join(map(str, [locus[0],
                                           locus[1],
                                           locus[2],
                                           locus[3],
                                           read,
                                           tstart,
                                           tend,
                                           len(seq)]))

        return header, '>%s\n%s\n' % (header, seq)

    def get_alleles(self, loci, flank=100, closeness_to_end=100):
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
        clipped = defaultdict(dict)
        for locus in loci:
            for aln in bam.fetch(locus[0], locus[1], locus[2]):
                locus_size = locus[2] - locus[1] + 1

                if INSFinder.is_uniquely_mapped(aln) and\
                   aln.reference_start <= locus[1] - flank and\
                   aln.reference_end >= locus[2] + flank:
                    seq, tstart, tend = self.extract_subseq(aln, locus[1] - flank, locus[2] + flank)
                    if seq is None:
                        print 'problem getting seq', aln.query_name, locus, len(aln.query_sequence)
                        continue
                    header, fa_entry = self.create_trf_fasta(locus, aln.query_name, tstart, tend, seq)
                    trf_input += fa_entry
                    repeat_seqs[header] = seq

                    # for finding position of repeat sequence inside read
                    read_seqs[aln.query_name] = aln.query_sequence
                    aligned_pairs[aln.query_name] = aln.get_aligned_pairs()

                elif self.check_split_alignments and aln.has_tag('SA') and\
                     ((aln.reference_start >= locus[1] - flank and aln.reference_start <= locus[2] + flank) or\
                      (aln.reference_end >= locus[1] - flank and aln.reference_end <= locus[2] + flank)):
                    clipped_end, partner_start = INSFinder.is_split_aln_potential_ins(aln)
                    if clipped_end is not None:
                        clipped[aln.query_name][clipped_end] = (aln, partner_start)

            INSFinder.find_missing_clipped_pairs(clipped, bam)

            # clipped alignment
            for read in clipped.keys():
                if len(clipped[read].keys()) == 2:
                    aln1 = clipped[read]['end'][0]
                    aln2 = clipped[read]['start'][0]

                    if aln1.query_alignment_end < aln2.query_alignment_start:
                        seq = aln1.query_sequence[aln1.query_alignment_end - locus_size - flank:aln2.query_alignment_start + locus_size + flank]
                        header, fa_entry = self.create_trf_fasta(locus, aln1.query_name, aln1.reference_end, aln2.reference_start, seq)
                        trf_input += fa_entry
                        repeat_seqs[header] = seq
                        read_seqs[read] = aln1.query_sequence

                    else:
                        del clipped[read]

                else:
                    del clipped[read]

        trf_fasta = create_tmp_file(trf_input)
        print trf_fasta
        output = self.run_trf(trf_fasta)

        results = self.parse_trf(output)

        # group by locus
        alleles = defaultdict(dict)
        for seq in results.keys():
            cols = seq.split(':')
            if len(cols) < 8:
                print 'problematic seq id: %s' % seq
                continue
            locus = tuple(cols[:4])

            expected_pats = cols[3].split(',')
            expected_pat_sizes = [len(p) for p in expected_pats]

            read = cols[4]
            gstart, gend = int(cols[5]), int(cols[6])

            results_matched = []
            for result in results[seq]:
                if len(result[13]) > 1:
                    for pat in expected_pats:
                        if self.is_same_repeat((pat, result[13])):
                            results_matched.append(result)
            # don't care if results' patterns are same as "expected"
            #results_matched = results[seq]

            if results_matched:
                combined_coords = self.combine_trf_coords([(r[0], r[1]) for r in results_matched])
                if len(combined_coords) == 1 and\
                   len(repeat_seqs[seq]) - 2 * flank > 0 and\
                   float(combined_coords[0][1] - combined_coords[0][0] + 1) / (len(repeat_seqs[seq]) - 2 * flank) > 0.8:
                    coords = sorted(list(itertools.chain.from_iterable([[int(result[0]), int(result[1])] for result in results_matched])))
                    repeat_seq = repeat_seqs[seq][coords[0]-1:coords[-1]]

                    size = len(repeat_seq)
                    cn = round(float(len(repeat_seq)) / min(expected_pat_sizes), 1)
                    rpos = read_seqs[read].find(repeat_seq)
                    rstart = rpos
                    rend = rstart + len(repeat_seq)

                    if clipped.has_key(read):
                        pairs = [clipped[read]['end'][0].get_aligned_pairs(), clipped[read]['start'][0].get_aligned_pairs()]
                    else:
                        pairs = aligned_pairs[read]
                    flanks_aligned = self.check_trf_prediction_fullness(rstart, rend, pairs, read=read)

                    if flanks_aligned:
                        d  = get_distance(int(locus[1]), int(locus[2]), flanks_aligned[0], flanks_aligned[1])
                        if d < 20:
                            if not alleles[locus].has_key(read) or \
                               len(alleles[locus][read][3]) < len(result[-1]):
                                alleles[locus][read] = (rpos, cols[3], cn, repeat_seq, flanks_aligned)

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
            if self.update_loci:
                variant.update_coords()

            # genotype
            variant.genotype()

        return variants

    def examine_ins(self, ins_list):
        if self.nprocs > 1:
            batches = list(split_tasks(ins_list, self.nprocs))
            batched_results = parallel_process(self.extract_tres, batches, self.nprocs)
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
                if len(ins.repeat_unit) <= self.max_str_len:
                    tre_events.append(ins)
            else:
                non_tre_events.append(ins)

        if tre_events:
            merged_loci = self.merge_loci(tre_events)

            variants = self.collect_alleles(merged_loci)

            # remove redundant variants
            if variants:
                self.remove_redundants(variants)
        
                return variants

        return []

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

    def collect_alleles(self, loci):
        tre_variants = []
        if self.nprocs > 1:
            batches = list(split_tasks(loci, self.nprocs))
            batched_results = parallel_process(self.get_alleles, batches, self.nprocs)
            if batched_results:
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

        return self.collect_alleles(loci)
    
    def output(self, variants, out_file):
        with open(out_file, 'w') as out:
            out.write('%s\n' % '\t'.join(Variant.tsv_headers + Allele.tsv_headers))
            for variant in sorted(variants, key=attrgetter('chrom', 'start', 'end')):
                if not variant.genotypes:
                    continue
                variant_cols = variant.as_tsv()
                for allele in sorted(variant.alleles, key=attrgetter('copy_number'), reverse=True):
                    allele_cols = allele.as_tsv()

                    out.write('%s\n' % '\t'.join(variant_cols + allele_cols))

    def cleanup(self):
        if self.tmp_files:
            for ff in self.tmp_files:
                if os.path.exists(ff):
                    os.remove(ff)
