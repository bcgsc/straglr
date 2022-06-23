import pysam
from pybedtools import BedTool
import sys
import os
from .utils import create_tmp_file, parallel_process, combine_batch_results, reverse_complement, split_tasks
from collections import defaultdict
import re
import random
import numpy as np

class INS:
    """
    chrom, gpos, gpos + 1, aln.query_name, rpos, ins_seq, 'ins', None
    instance is a list:
    0: chrom
    1: gstart
    2: gend
    3: read
    4: rpos
    5: seq
    6: etype ('ins' or 'tre')
    7: neighbour_seq
    8: repeat_unit (for tre)
    """
    @classmethod
    def extract_neighbour_seqs(cls, read_seq, rpos, ins_size, flank_size):
        return read_seq[rpos - flank_size : rpos] + \
               read_seq[rpos + ins_size : rpos + ins_size + flank_size]

    @classmethod
    def eid(cls, ins):
        return '_'.join(map(str, [ins[3], ins[4], len(ins[5])]))

    @classmethod
    def to_bed(cls, ins, include_read=False):
        cols = [ins[0], ins[1], ins[1] + 1, cls.eid(ins)]
        if include_read:
            cols.append(ins[3])
        return '\t'.join(map(str, cols))

class INSFinder:
    def __init__(self, bam, genome_fasta, min_size, w=100, reads_fasta=None, exclude=None, chroms=None, nprocs=None, min_support=None, max_cov=100, debug=False):
        self.bam = bam
        self.genome_fasta = genome_fasta

        self.min_size = min_size

        self.reads_fasta = reads_fasta
        
        # for extracting neighbor sequence - for tre search
        self.w = w

        # bed file of exclude regions
        self.exclude = exclude
        
        self.chroms = chroms
        self.nprocs = nprocs

        self.check_split_alignments = True

        self.min_support = min_support

        self.max_cov = max_cov

        self.debug = debug

        self.tmp_files = set()

    def find_ins(self, regions_bed_file=None):
        bam = pysam.Samfile(self.bam, 'rb')

        if regions_bed_file is None:
            chroms = [chrom for chrom in bam.references if not '_' in chrom and not 'M' in chrom]
            if self.chroms:
                chroms = self.chroms

            regions = self.get_examine_regions(bam, chroms)
        else:
            regions = self.get_examine_regions(None, None, regions_bed_file=regions_bed_file)

        all_ins = []
        if regions:
            if self.nprocs > 1:
                random.shuffle(regions)
                batches = list(split_tasks(regions, self.nprocs))
                batched_results = parallel_process(self.examine_regions, batches, self.nprocs)
                all_ins = combine_batch_results(batched_results, type(batched_results[0]))
            else:
                reads_fasta = []
                if self.reads_fasta:
                    for fa in self.reads_fasta:
                        reads_fasta.append(pysam.Fastafile(fa))
                all_ins = []
                for region in regions:
                    all_ins.extend(self.examine_region(region, bam=bam, reads_fasta=reads_fasta))

        ins_filtered = []
        if all_ins:
            ins_cleared = all_ins
            if self.exclude:
                ins_cleared = self.screen(all_ins)

            if ins_cleared:
                eids_merged = self.merge(ins_cleared)

                if eids_merged:
                    ins_filtered = [ins for ins in ins_cleared if INS.eid(ins) in eids_merged]

        if not self.debug:
            self.cleanup()

        return ins_filtered

    def cleanup(self):
        if self.tmp_files:
            for ff in self.tmp_files:
                if os.path.exists(ff):
                    os.remove(ff)

    def has_secondary_alignment(self, aln):
        try:
            sa = aln.get_tag('SA')
            return sa
        except:
            return False

    def get_mapped_fractions(self, aln):
        mapped_size = sum([t[1] for t in aln.cigartuples if t[0] == 0])
        fraction_query = float(mapped_size) / (aln.query_alignment_end - aln.query_alignment_start + 1)
        fraction_reference = float(mapped_size) / (aln.reference_end - aln.reference_start + 1)
        return fraction_query, fraction_reference

    @classmethod
    def is_uniquely_mapped(cls, aln, max_clipping_allowed=0.05):
        closeness_to_end = aln.infer_read_length() * max_clipping_allowed
        if float(aln.query_alignment_length) / aln.infer_read_length() > 0.5 and not aln.has_tag('SA'):
            return True
        if aln.query_alignment_start <= closeness_to_end and\
           aln.query_alignment_end >= aln.query_length - closeness_to_end and\
           not aln.has_tag('SA'):
            return True

    @classmethod
    def get_seq(cls, fastas, name, reverse, coords=None):
        try:
            for fasta in fastas:
                seq = cls.get_seq_from_fasta(fasta, name, reverse, coords)
                if seq is not None:
                    return seq
        except:
            return None

    @classmethod
    def get_seq_from_fasta(cls, fasta, name, reverse, coords=None):
        try:
            seq = fasta.fetch(name)
            if reverse:
                seq = reverse_complement(seq)
            if coords:
                return seq[coords[0]:coords[1]]
            else:
                return seq
        except:
            return None

    def examine_regions(self, regions):
        """worker wrapper for examine_region"""
        bam = pysam.Samfile(self.bam, 'rb')
        reads_fasta = []
        if self.reads_fasta:
            for fa in self.reads_fasta:
                reads_fasta.append(pysam.FastaFile(fa))

        ins_list = []
        for region in regions:
            ins_list.extend(self.examine_region(region, bam=bam, reads_fasta=reads_fasta))

        return ins_list

    def get_coverage(self, bam, region):
        covs = []
        for pileup_col in bam.pileup(region[0], int(region[1]), int(region[2])):
            covs.append(pileup_col.n)

        if covs:
            return np.mean(covs)
        else:
            return None

    def examine_region(self, region, bam=None, reads_fasta=None):
        if bam is None:
            bam = pysam.Samfile(self.bam, 'rb')

        ins_list = []
        mean_cov = self.get_coverage(bam, region)
        if mean_cov is None or mean_cov > self.max_cov:
            return ins_list

        clipped_pairs = defaultdict(dict)
        read_spans = {}
        for aln in bam.fetch(region[0], int(region[1]), int(region[2])):
            if not reads_fasta and not aln.query_sequence:
                continue

            if not aln.query_name in read_spans:
                read_spans[aln.query_name] = {'starts':defaultdict(int), 'ends':defaultdict(int)}
            read_spans[aln.query_name]['starts'][aln.reference_start] += 1
            read_spans[aln.query_name]['ends'][aln.reference_end] += 1

            if self.check_split_alignments:
                clipped_end, partner_start = self.is_split_aln_potential_ins(aln)
                if clipped_end is not None:
                    clipped_pairs[aln.query_name][clipped_end] = (aln, partner_start)

            if not aln.query_name in clipped_pairs and self.is_uniquely_mapped(aln):
                ins_list.extend(self.extract_ins(aln, region, reads_fasta=reads_fasta))

        if clipped_pairs:
            ins_list.extend(self.extract_ins_from_clipped(clipped_pairs, bam, read_spans, reads_fasta=reads_fasta))

        return ins_list

    def extract_ins_from_clipped(self, clipped_pairs, bam, read_spans, reads_fasta=None):
        ins_from_clipped = []
        for read in clipped_pairs.keys():
            if 'start' in clipped_pairs[read] and 'end' in clipped_pairs[read]:
                aln1, aln2 = clipped_pairs[read]['end'][0], clipped_pairs[read]['start'][0]

                if aln2.reference_end <= aln1.reference_end or aln1.reference_start >= aln2.reference_start:
                    continue
                if read_spans[read]['starts'][aln1.reference_start] > 1 or read_spans[read]['ends'][aln2.reference_end] > 1:
                    continue

                if not reads_fasta:
                    ins_seq = aln1.query_sequence[aln1.query_alignment_end:aln2.query_alignment_start]
                else:
                    # need to adjust coordinates for hard-clips
                    qstart = aln1.query_alignment_end
                    qend = aln2.query_alignment_start
                    if aln1.cigartuples[0][0] == 5:
                        qstart += aln1.cigartuples[0][1]
                    if aln2.cigartuples[0][0] == 5:
                        qend += aln2.cigartuples[0][1]
                    ins_seq = self.get_seq(reads_fasta, read, aln1.is_reverse, [qstart, qend])

                if len(ins_seq) >= self.min_size:
                    if aln1 != aln2:
                        mid = int((aln1.reference_end + aln2.reference_start) / 2)
                        ins_from_clipped.append([aln1.reference_name,
                                                 mid,
                                                 mid + 1,
                                                 aln1.query_name,
                                                 aln1.query_alignment_end,
                                                 ins_seq,
                                                 'ins',
                                                 None],
                                                )

        return ins_from_clipped

    @classmethod
    def is_split_aln_potential_ins(cls, aln, closeness_to_end=200, min_split_size=500, check_end=None, use_sa=True):
        clipped_end = None
        partner_start = None
        if use_sa and aln.has_tag('SA'):
            clipped_end, partner_start = cls.is_split_aln_potential_ins_with_SA(aln,
                                                                                closeness_to_end=closeness_to_end,
                                                                                min_split_size=min_split_size,
                                                                                check_end=check_end)
        if clipped_end is None:
            clipped_end, partner_start = cls.is_split_aln_potential_ins_without_SA(aln,
                                                                                closeness_to_end=closeness_to_end,
                                                                                min_split_size=min_split_size,
                                                                                check_end=check_end)
        return clipped_end, partner_start

    @classmethod
    def is_split_aln_potential_ins_without_SA(cls, aln, closeness_to_end=200, min_split_size=500, check_end=None):
        def is_clipped(end):
            if end == 'start':
                i, j = 0, -1
            else:
                i, j = -1, 0

            if aln.cigartuples[i][0] >= 4 and aln.cigartuples[i][0] <=5 and\
                aln.cigartuples[i][1] >= min_split_size:
                if not check_end:
                    if not (aln.cigartuples[j][0] >= 4 and aln.cigartuples[j][0] <=5 and\
                            aln.cigartuples[j][1] >= min_split_size) and\
                        aln.cigartuples[j][1] <= closeness_to_end:
                            return True
                else:
                    return True
                return False

        clipped_end = None
        partner_start = None

        if (check_end is None or check_end == 'start') and is_clipped('start'):
            clipped_end = 'start'
            partner_start = aln.reference_start
        elif (check_end is None or check_end == 'end') and is_clipped('end'):
            clipped_end = 'end'
            partner_start = aln.reference_end

        return clipped_end, partner_start

    @classmethod
    def is_split_aln_potential_ins_with_SA(cls, aln, closeness_to_end=200, min_split_size=500, closeness_in_size=40000, check_end=None):
        """ aln must have been checked if it has SA """
        clipped = None
        partner_start = None
        clipped_size1 = None
        clipped_size2 = None

        clipped_start_regex = re.compile('^(\d+)[S|H]')
        clipped_end_regex = re.compile('(\d+)[S|H]$')

        closeness_to_each = aln.infer_read_length()
        sas = cls.get_secondary_alignments(aln)
        if aln.reference_name in sas:
            strand1 = '+' if not aln.is_reverse else '-'
            if (check_end is None or check_end == 'end') and\
                aln.cigartuples[-1][0] >= 4 and\
                aln.cigartuples[-1][1] >= min_split_size and\
                aln.query_alignment_start <= closeness_to_end:
                clipped_size1 = aln.cigartuples[-1][1]
                for sa in sas[aln.reference_name]:
                    if abs(sa[0] - aln.reference_start) < closeness_to_each and\
                       sa[1] == strand1:
                        match_start = clipped_start_regex.search(sa[2])
                        match_end = clipped_end_regex.search(sa[2])
                        if match_start and\
                           (match_end is None or int(match_end.group(1)) <= closeness_to_end):
                            clipped_size2 = int(match_start.group(1))
                            if abs(clipped_size1 - clipped_size2) <= closeness_in_size:
                                clipped = 'end'
                                partner_start = sa[0]
                                break

            if clipped is None and (check_end is None or check_end == 'start'):
                if aln.cigartuples[0][0] >= 4 and\
                   aln.cigartuples[0][1] >= min_split_size and\
                   aln.query_alignment_end >= aln.query_length - closeness_to_end:
                    clipped_size1 = aln.cigartuples[0][1]
                    for sa in sas[aln.reference_name]:
                        if abs(sa[0] - aln.reference_start) < closeness_to_each and\
                           sa[1] == strand1:
                            match_end = clipped_end_regex.search(sa[2])
                            match_start = clipped_start_regex.search(sa[2])
                            if match_end and\
                               (match_start is None or int(match_start.group(1)) <= closeness_to_end):
                                clipped_size2 = int(match_end.group(1))
                                if abs(clipped_size1 - clipped_size2) <= closeness_in_size:
                                    clipped = 'start'
                                    partner_start = sa[0]
                                    break

        return clipped, partner_start

    @classmethod
    def get_secondary_alignments(cls, aln):
        mappings = defaultdict(list)
        if aln.has_tag('SA'):
            for mapping in aln.get_tag('SA').split(';'):
                if mapping:
                    chrom, start, strand, cigar, mapq, nm = mapping.split(',')
                    mappings[chrom].append([int(start), strand, cigar, int(mapq), int(nm)])

        return mappings
        
    def extract_ins(self, aln, region, reads_fasta=None, out=None, w=100):
        ins_found = []

        ins_tuples = [ct for ct in aln.cigartuples if ct[0] == 1 and ct[1] >= self.min_size]
        chrom = aln.reference_name
        for tup in ins_tuples:
            size = tup[1]
            index_tup = aln.cigartuples.index(tup)

            offset_query = sum([c[1] for c in aln.cigartuples[:index_tup] if c[0] in (0, 1, 4)])
            offset_target = sum([c[1] for c in aln.cigartuples[:index_tup] if c[0] in (0, 2, 3)])
            rpos = offset_query
            gpos = aln.reference_start + offset_target - 1
            if gpos < int(region[1]) and gpos > int(region[2]):
                continue

            if rpos is not None:
                if not reads_fasta:
                    ins_seq = aln.query_sequence[rpos:rpos + size]
                else:
                    ins_seq = self.get_seq(reads_fasta, aln.query_name, aln.is_reverse, [rpos, rpos + size])
                ins = [chrom, gpos, gpos + 1, aln.query_name, rpos, ins_seq, 'ins', None]
                if not self.reads_fasta:
                    ins[7] = INS.extract_neighbour_seqs(aln.query_sequence, rpos, len(ins_seq), self.w)
                else:
                    ins[7] = INS.extract_neighbour_seqs(self.get_seq(reads_fasta, aln.query_name, aln.is_reverse), rpos, len(ins_seq), self.w)
                ins_found.append(ins)
                        
        return ins_found

    def extract_neighbour_seqs(self, read_seq, rpos, size, ins_len):
        return read_seq[rpos - size : rpos] + \
               read_seq[rpos + ins_len : rpos + ins_len + size]
 
    def create_tmp_bed(self, ins_list, include_read=False):
        ins_bed = ''
        for ins in ins_list:
            ins_bed += '{}\n'.format(INS.to_bed(ins, include_read=include_read))

        ins_bed_file = create_tmp_file(ins_bed)
        self.tmp_files.add(ins_bed_file)

        return ins_bed_file

    def merge(self, ins_list, d=50):
        ins_bed_file = self.create_tmp_bed(ins_list, include_read=True)
        if self.debug:
            print('ins all {}'.format(ins_bed_file))

        ins_all = BedTool(ins_bed_file)
        ins_merged = ins_all.sort().merge(d=d, c='4,5', o='distinct,count_distinct')
        if self.debug:
            ins_merged.saveas('ins_merged.bed')

        eids_merged = set()
        for c in ins_merged:
            if int(c[4]) >= self.min_support:
                eids_merged |= set(c[3].split(','))

        return eids_merged

    def screen(self, ins_list):
        ins_bed_file = self.create_tmp_bed(ins_list)

        ins_all = BedTool(ins_bed_file)
        exclusions = BedTool(self.exclude)

        ins_removed = ins_all.intersect(self.exclude)
        eids_removed = set([ins[3] for ins in ins_removed])
        ins_kept = [ins for ins in ins_list if INS.eid(ins) not in eids_removed]
        
        return ins_kept

    def split_region(self, region, max_size):
        subregions = []
        start = int(region[1])
        end = int(region[2])
        while start + max_size - 1 < end:
            subregions.append([region[0], start, start + max_size - 1])
            start = start + max_size
        subregions.append([region[0], start, end])
        return subregions

    def get_examine_regions(self, bam, chroms, chunk_size=1000000, regions_bed_file=None):
        regions = []
        if regions_bed_file:
            regions_bed = BedTool(regions_bed_file)
            for region in regions_bed.window_maker(regions_bed, w=chunk_size, s=chunk_size+1):
                regions.append((region[0], int(region[1]), int(region[2])))

        else:
            chrom_lens = dict((bam.references[i], bam.lengths[i]) for i in range(len(bam.references)))

            if self.exclude:
                chrom_spans = ''
                for chrom in sorted(chroms):
                    chrom_spans += '{}\n'.format('\t'.join(map(str, [chrom, 0, chrom_lens[chrom]])))
                chrom_bed_file = create_tmp_file(chrom_spans)
                self.tmp_files.add(chrom_bed_file)

                chrom_bed = BedTool(chrom_bed_file)
                exclude_bed = BedTool(self.exclude)
                regions_kept = map(tuple, chrom_bed.subtract(exclude_bed))

                for region in regions_kept:
                    regions.extend(self.split_region(region, chunk_size))

            else:
                for chrom in chroms:
                    i = 0
                    while i < chrom_lens[chrom]:
                        regions.append((chrom, i, min(i + chunk_size, chrom_lens[chrom])))
                        i += chunk_size + 1
 
        return regions
