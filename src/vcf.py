class VCF:
    file_format = '4.2'
    header = 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'
    info = (
            ('LOCUS', 1, 'String', 'Locus ID'),
            ('END', 1, 'Integer', 'End position of repeat'), 
            ('MOTIF', 1, 'String', 'Reference motif'),
            ('COPIES', 1, 'Float', 'Reference copies')
            )
    format = (
              ('GT', 1, 'String', 'Genotype'),
              ('DP', 1, 'Integer', 'Read depth'),
              ('AL', '.', 'String', 'Allelic lengths'),
              ('ALR', '.', 'String', 'Allelic length ranges'),
              ('AC', '.', 'String', 'Allelic copies'),
              ('ACR', '.', 'String', 'Allelic copy ranges'),
              ('AD', '.', 'String', 'Allelic depths'),
              ('ALT_MOTIF', '.', 'String', 'Alternate motif(s)'),
             )
    filters = {
        ('UNMATCHED_MOTIF', 'Unmatched motif'): ('unmatched_motif'),
        ('PARTIAL_SPAN', 'Repeat not fully span locus'): ('partial_and_insufficient_span',
                                                         'insufficient_repeat_coverage',
                                                         'too_far_from_flank'),
        ('INVALID_MOTIF', 'Invalid motif'): ('motif_size_out_of_range'),
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

        # header
        cols = list(cls.header)
        if sample is not None:
            cols.append(sample)
        lines.append('#{}'.format('\t'.join(cols)))

        return '\n'.join(lines)
    
    @classmethod
    def extract_variant_info(cls, variant, locus_id):
        pairs = []
        if locus_id is not None:
            pairs.append('LOCUS={}'.format(locus_id))

        # 2 == chr, 8 == repeat
        for key, i in zip(cls.info[1:], (2,8,10)):
            pairs.append('{}={}'.format(key[0], variant[i]))

        return ';'.join(pairs)

    @classmethod
    def extract_variant_gt(cls, variant, gt_size, gt_cn, supports, size_ranges, cn_ranges, alt_motifs):
        vals = {}

        # 5 = coverage
        gts = sorted(list(set([a[9] for a in variant[3] if a[9] is not None and a[-2] == 'full'])))
        if gts:
            vals['GT'] = '/'.join(map(str, gts))

        for label, val in zip(('AL', 'AC', 'AD', 'ALR', 'ACR'), (gt_size, gt_cn, supports, size_ranges, cn_ranges)):
            if val:
                val_array = val
                if len(val) == 1:
                    val_array.append(val[0])
                vals[label] = '/'.join(map(str, val_array))

        if variant[5]:
            vals['DP'] = str(variant[5])

        if alt_motifs:
            vals['ALT_MOTIF'] = alt_motifs

        format = []
        for fmt in cls.format:
            if fmt[0] in vals:
                format.append(fmt[0])

        return [':'.join(format), ':'.join([vals[field[0]] for field in cls.format if field[0] in vals])]

    @classmethod
    def extract_filter(cls, mess):
        for filter in cls.filters:
            if mess in cls.filters[filter]:
                return filter
        return '.'
