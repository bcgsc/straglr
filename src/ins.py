import pysam
from expanded_repeat_finder import ExpandedRepeatFinder
from pybedtools import BedTool
import sys
import os
from utils import create_tmp_file, parallel_process, combine_batch_results
from intspan import intspan
from collections import defaultdict
import re

class INS:
    def __init__(self, chrom, gpos, read, rpos, seq, etype):
        self.chrom = chrom
        self.gpos = gpos
        self.gend = gpos + 1
        self.read = read
        self.rpos = rpos
        self.seq = seq
        self.etype = etype
        self.eid = '_'.join(map(str, [self.read, self.rpos, len(self.seq)]))

        self.neighbour_seq = None
        self.repeat_unit = None

    def extract_neighbour_seqs(self, read_seq, size):
        self.neighbour_seq = read_seq[self.rpos - size : self.rpos] + \
            read_seq[self.rpos + len(self.seq) : self.rpos + len(self.seq) + size]

    def as_bed(self):
        cols = (self.chrom, self.gpos, self.gpos + 1, self.eid)
        return '\t'.join(map(str, cols))

class INSFinder:
    def __init__(self, bam, genome_fasta, min_size, w=100, exclude=None, chroms=None, nprocs=None):
        self.bam = bam
        self.genome_fasta = genome_fasta

        self.min_size = min_size
        
        # for extracting neighbor sequence - for tre search
        self.w = w

        # bed file of exclude regions
        self.exclude = exclude
        
        self.erf = ExpandedRepeatFinder(self.bam, self.genome_fasta)
        self.erf.flank_len = self.w

        self.chroms = chroms
        self.nprocs = nprocs

        self.erf = ExpandedRepeatFinder(self.bam, self.genome_fasta)
        self.erf.flank_len = self.w

        self.check_split_alignments = False

    def find_ins(self):
        bam = pysam.Samfile(self.bam, 'rb')
        chroms = [chrom for chrom in bam.references if not '_' in chrom and not 'M' in chrom]
        if self.chroms:
            chroms = self.chroms

        regions = self.get_examine_regions(chroms)
        all_ins = []
        if regions:
            if self.nprocs > 1:
                batched_results = parallel_process(self.examine_region, regions, self.nprocs)
                all_ins = combine_batch_results(batched_results, type(batched_results[0]))
            else:
                all_ins = []
                for region in regions:
                    all_ins.extend(self.examine_region(region))
                    
        ins_cleared = all_ins
        if self.exclude:
            ins_cleared = self.screen(all_ins)        

        return ins_cleared

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
    def is_uniquely_mapped(cls, aln, closeness_to_end=100):
        if aln.query_alignment_start <= closeness_to_end and\
           aln.query_alignment_end >= aln.query_length - closeness_to_end and\
           not aln.has_tag('SA'):
            return True

    def examine_region(self, region):
        bam = pysam.Samfile(self.bam, 'rb')

        #skipped = set()
        ins_list = []
        clipped_pairs = defaultdict(dict)
        for aln in bam.fetch(region[0], int(region[1]), int(region[2])):
            #mapped_fractions = self.get_mapped_fractions(aln)
            #if aln.mapping_quality == 0 or self.has_secondary_alignment(aln):
                #skipped.add(aln.query_name)

            #if aln.has_tag('SA'):
                #clipped_end, partner_start = self.is_split_aln_potential_ins(aln)
                #if clipped_end is not None:
                    #clipped_pairs[aln.query_name][clipped_end] = (aln, partner_start)
            #else:
                #ins_list.extend(self.extract_ins(aln))

            if self.is_uniquely_mapped(aln):
                ins_list.extend(self.extract_ins(aln))
            elif self.check_split_alignments and aln.has_tag('SA'):
                clipped_end, partner_start = self.is_split_aln_potential_ins(aln, min_split_size=self.min_size)
                if clipped_end is not None:
                    clipped_pairs[aln.query_name][clipped_end] = (aln, partner_start)

        if clipped_pairs:
            ins_list.extend(self.extract_ins_from_clipped(clipped_pairs, bam))

        return ins_list

    @classmethod
    def find_missing_clipped_pairs(cls, clipped_pairs, bam):
        for read in clipped_pairs.keys():
            missing_end = None
            chrom = None
            coord = None
            partner_start = max
            if clipped_pairs[read].has_key('start') and not clipped_pairs[read].has_key('end'):
                missing_end = 'end'
                chrom = clipped_pairs[read]['start'][0].reference_name
                partner_start = clipped_pairs[read]['start'][0].reference_start
                coord = clipped_pairs[read]['start'][1]
            elif clipped_pairs[read].has_key('end') and not clipped_pairs[read].has_key('start'):
                missing_end = 'start'
                chrom = clipped_pairs[read]['end'][0].reference_name
                partner_start = clipped_pairs[read]['end'][0].reference_start
                coord = clipped_pairs[read]['end'][1]

            if missing_end:
                for aln in bam.fetch(chrom, coord-1, coord+1):
                    if aln.query_name == read and aln.reference_start == coord - 1:
                        clipped_pairs[read][missing_end] = (aln, partner_start)

    def extract_ins_from_clipped(self, clipped_pairs, bam):
        # find missing partners
        self.find_missing_clipped_pairs(clipped_pairs, bam)

        ins_from_clipped = []
        for read in clipped_pairs.keys():
            if clipped_pairs[read].has_key('start') and clipped_pairs[read].has_key('end'):
                aln1, aln2 = clipped_pairs[read]['end'][0], clipped_pairs[read]['start'][0]
                ins_seq = aln1.query_sequence[aln1.query_alignment_end:aln2.query_alignment_start]
                if len(ins_seq) >= self.min_size:
                    if aln1 != aln2:
                        print 'iii', aln1.query_name, aln1.reference_name, aln1.reference_start, aln2.reference_name, aln2.reference_start, ins_seq
                        ins_from_clipped.append(INS(aln1.reference_name,
                                                    aln1.reference_end,
                                                    aln1.query_name,
                                                    aln1.query_alignment_end,
                                                    ins_seq,
                                                    'ins'))
                else:
                    print 'iij', aln1.query_name, aln1.reference_name, aln1.reference_start, aln2.reference_name, aln2.reference_start

        return ins_from_clipped

    @classmethod
    def is_split_aln_potential_ins(cls, aln, closeness_to_end=100, min_split_size=500, closeness_to_each=5000, closeness_in_size=20000, max_outside_mappings=2):
        """ aln must have been checked if it has SA """
        clipped = None
        partner_start = None
        clipped_size1 = None
        clipped_size2 = None

        clipped_start_regex = re.compile('^(\d+)[S|H]')
        clipped_end_regex = re.compile('(\d+)[S|H]$')

        sas = cls.get_secondary_alignments(aln)
        if sas.has_key(aln.reference_name) and len(sas[aln.reference_name]) <= max_outside_mappings:
            strand1 = '+' if not aln.is_reverse else '-'
            if aln.cigartuples[-1][0] >= 4 and\
               aln.cigartuples[-1][1] >= min_split_size and\
               aln.query_alignment_start <= closeness_to_end:
                clipped_size1 = aln.cigartuples[-1][1]
                for sa in sas[aln.reference_name]:
                    if abs(sa[0] - aln.reference_start) < closeness_to_each and\
                       sa[1] == strand1:
                        #match = re.search('^(\d+)[S|H]', sa[2])
                        match_start = clipped_start_regex.search(sa[2])
                        match_end = clipped_end_regex.search(sa[2])
                        if match_start and\
                           (match_end is None or int(match_end.group(1)) <= closeness_to_end):
                            clipped_size2 = int(match_start.group(1))
                            if abs(clipped_size1 - clipped_size2) <= closeness_in_size:
                                clipped = 'end'
                                partner_start = sa[0]
                                #print 'ee', aln.query_name, aln.cigarstring, clipped_size2, sa[2], match_end
                                break

            elif aln.cigartuples[0][0] >= 4 and\
                 aln.cigartuples[0][1] >= min_split_size and\
                 aln.query_alignment_end >= aln.query_length - closeness_to_end:
                clipped_size1 = aln.cigartuples[0][1]
                for sa in sas[aln.reference_name]:
                    if abs(sa[0] - aln.reference_start) < closeness_to_each and\
                       sa[1] == strand1:
                        #match = re.search('(\d+)[S|H]$', sa[2])
                        match_end = clipped_end_regex.search(sa[2])
                        match_start = clipped_start_regex.search(sa[2])
                        if match_end and\
                           (match_start is None or int(match_start.group(1)) <= closeness_to_end):
                            clipped_size2 = int(match_end.group(1))
                            if abs(clipped_size1 - clipped_size2) <= closeness_in_size:
                                clipped = 'start'
                                partner_start = sa[0]
                                #print 'ss', aln.query_name, aln.cigarstring, clipped_size2, sa[2]
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
        
    def extract_ins(self, aln, out=None, w=100):
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

            if rpos is not None:
                ins_seq = aln.query_sequence[rpos:rpos + size]
                ins = INS(chrom, gpos, aln.query_name, rpos, ins_seq, 'ins')
                ins.extract_neighbour_seqs(aln.query_sequence, self.w)
                ins_found.append(ins)
                        
        return ins_found

    def extract_neighbour_seqs(self, read_seq, rpos, size, ins_len):
        return read_seq[rpos - size : rpos] + \
               read_seq[rpos + ins_len : rpos + ins_len + size]

    def screen(self, ins_list):
        ins_bed = ''
        for ins in ins_list:
            ins_bed += '%s\n' % ins.as_bed()

        ins_bed_file = create_tmp_file(ins_bed)

        ins_all = BedTool(ins_bed_file)
        exclusions = BedTool(self.exclude)
        ins_cleared = ins_all.subtract(exclusions, A=True)

        eids_cleared = set([c[3] for c in ins_cleared])
        return [ins for ins in ins_list if ins.eid in eids_cleared]

    def get_excluded_regions(self, chrom):
        skips = intspan()
        with open(self.exclude, 'r') as ff:
            for line in ff:
                cols = line.rstrip().split()
                if cols[0] == chrom:
                    skips |= '%s-%s' % (cols[1], cols[2])

        return skips

    def get_examine_regions(self, chroms):
        bam = pysam.Samfile(self.bam, 'rb')
        chrom_lens = dict((bam.references[i], bam.lengths[i]) for i in range(len(bam.references)))
        regions = []
        if self.exclude:
            chrom_spans = ''
            for chrom in sorted(chroms):
                chrom_spans += '%s\n' % '\t'.join(map(str, [chrom, 0, chrom_lens[chrom]]))
            chrom_bed_file = create_tmp_file(chrom_spans)

            chrom_bed = BedTool(chrom_bed_file)
            exclude_bed = BedTool(self.exclude)
            regions = map(tuple, chrom_bed.subtract(exclude_bed))

        else:
            size = 1000000
            for chrom in chroms:
                i = 0
                while i < chrom_lengths[chrom]:
                    regions.append((chrom, i, min(i + size, chrom_lens[chrom])))
                    print spans[-1]
                    i += size + 1

        return regions

