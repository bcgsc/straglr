import re
from collections import defaultdict, Counter
from operator import itemgetter

class VCF:
    file_format = '4.5'
    header = 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'
    info = (
            ('LOCUS', 1, 'String', 'Locus ID'),
            ('RUS_REF', 1, 'String', 'Repeat unit sequence in the reference orientation'),
            ('SVLEN', 'A', 'Integer', 'Length of structural variant'),
            ('RN', 'A', 'Integer', 'Total number of repeat sequences in this allele'),
            ('RUS', '.', 'String', 'Repeat unit sequence of the corresponding repeat sequence'),
            ('RUL', '.', 'Integer', 'Repeat unit length of the corresponding repeat sequence'),
            ('RUC', '.', 'Float', 'Repeat unit count of corresponding repeat sequence'),
            ('RUB', '.', 'Integer', 'Number of bases in each individual repeat unit'),
            ('RB', '.', 'Integer', 'Total number of bases in the corresponding repeat sequence'),
            ('CIRUC', '.', 'Float', 'Confidence interval around RUC'),
            ('CIRB', '.' 'Integer', 'Confidence interval around RB'),
            )
    format = (
              ('GT', 1, 'String', 'Genotype'),
              ('DP', 1, 'Integer', 'Read depth'),
             )
    filters = {
        ('UNMATCHED_MOTIF', 'Unmatched motif'): ('unmatched_motif'),
        ('PARTIAL_SPAN', 'Repeat not fully span locus'): ('partial_and_insufficient_span',
                                                         'insufficient_repeat_coverage',
                                                         'too_far_from_flank'),
        ('INVALID_MOTIF', 'Invalid motif'): ('motif_size_out_of_range'),
        ('CLUSTERING_FAILED', 'Clustering failed to formulate genotype'): ('clustering_failed'),
         }

    @classmethod
    def show_info_format(cls, name, data):
        lines = []
        keys = ('ID', 'Number', 'Type', 'Description')
        for row in data:
            pairs = []
            for key, col in zip(keys, row):
                if key == 'Description':
                    col = '"{}"'.format(col)
                pairs.append('{}={}'.format(key, col))
            lines.append('##{}=<{}>'.format(name, ','.join(pairs)))
        return '\n'.join(lines)

    @classmethod
    def show_meta(cls, sample, num_passes, contigs, ref=None, source=None, date=None, fails=None):
        lines = []
        # fileformat
        lines.append('##fileformat=VCFv{}'.format(cls.file_format))

        # fileDate
        if date is not None:
            lines.append('##fileDate={}'.format(date))

        # source
        if source is not None:
            lines.append('##source={}'.format(source))

        # reference
        if ref is not None:
            lines.append('##reference={}'.format(ref))

        # contigs
        if contigs:
            for contig, length in contigs:
                lines.append('##contig=<ID={},length={}>'.format(contig, length))
        # info
        lines.append(cls.show_info_format('INFO', cls.info))

        # filter
        filters = []
        if num_passes > 0:
            filters.append(('PASS', 'All filters passed'))
        if fails is not None:
            for mess in set(fails.values()):
                for filter in cls.filters:
                    if mess in cls.filters[filter]:
                        filters.append(filter)
        if filters:
            for id, desc in filters:
                lines.append('##FILTER=<ID={},Description="{}">'.format(id, desc))

        # format
        lines.append(cls.show_info_format('FORMAT', cls.format))

        # alt
        lines.append('##ALT=<ID=CNV:TR,Description="Tandem repeat determined base on DNA abundance">')

        # header
        cols = list(cls.header)
        if sample is not None:
            cols.append(sample)
        lines.append('#{}'.format('\t'.join(cols)))

        return '\n'.join(lines)
    
    @classmethod
    def extract_variant_genotypes(cls, variant, locus_id, genotype_in_size):
        def extract_dps():
            dps = {}
            for ad in variant[7].split(';'):
                a, d = ad.rstrip(')').split('(')
                dps[a] = d
            return dps
        
        # needs to be re-done for complex motif
        def extract_motifs(alleles, cn=None, size=None):
            ''' [(motif1, cn1), (motif2, cn2)]'''
            motifs = [a[2] for a in alleles]
            motif = Counter(motifs).most_common(1)[0][0]
            if cn is None and size is not None:
                cn = size / len(motif)
            return [(motif, cn)]

        ref_len = variant[2] - variant[1]
        info = defaultdict(list)
        genotype = defaultdict(list)
        if locus_id is not None:
            info['LOCUS'] = (locus_id,)
        info['RUS_REF'] = (variant[8],)
        
        dps = extract_dps()
        gts = defaultdict(list)
        for a in variant[3]:
            if a[9] is not None:
                gts[(a[9], a[11])].append(a)
        print('aa', locus_id, genotype_in_size, variant[4], variant[7], dps, len(gts), gts.keys())
        gts_sorted = sorted(gts.keys(), key=itemgetter(0))
        alts = []
        for gt, allele in gts_sorted:
            if gt > 0:
                alts.append('<CNV:TR>')
                cn, size = None, None
                if not genotype_in_size:
                    cn = allele
                else:
                    size = allele
                motifs = extract_motifs(gts[(gt, allele)], cn=cn, size=size)
                print('mm', motifs)
                
                info['RN'].append(len(motifs))
                info['RUS'].extend([m[0] for m in motifs])
                info['RB'].extend(['{:.0f}'.format(len(m[0]) * m[1]) for m in motifs])
                info['SVLEN'].append(ref_len)

            genotype['GT'].append(gt)
            genotype['DP'].append(dps[str(allele)])

        # INFO
        vals = []
        for i in cls.info:
            if i[0] in info:
                print('qq', i[0], info[i[0]])
                vals.append('{}:{}'.format(i[0], ','.join(map(str, info[i[0]]))))
        info_col = ';'.join(vals)
        print('ii', variant[:3], vals, info_col)
        
        # FORMAT
        vals = []
        for f in cls.format:
            if f[0] in genotype:
                vals.append((f[0], map(str, genotype[f[0]])))
        format_cols = ['\t'.join((':'.join([v[0] for v in vals]), ':'.join(['/'.join(v[1]) for v in vals])))]
        print('ff',  variant[:3], vals, format_cols)

        return ','.join(alts), info_col, format_cols

    @classmethod
    def extract_filter(cls, mess):
        for filter in cls.filters:
            if mess in cls.filters[filter]:
                return filter
        return '.'
