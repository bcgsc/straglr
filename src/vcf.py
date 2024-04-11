class VCF:
    file_format = '4.2'
    header = 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'
    info = (
            ('LOCUS', 1, 'String', 'Locus ID'),
            ('END', 1, 'Integer', 'End position of repeat'), 
            ('MOTIF', 1, 'String', 'Reference motif')
            )
    format = (
              ('GT', 1, 'String', 'Genotype'),
              ('DP', 1, 'Integer', 'Read depth'),
              ('AL', '.', 'String', 'Allelic length(s)'),
              ('AD', '.', 'String', 'Allelic depth(s)'),
              ('AM', '.', 'String', 'Alternate motif(s)'),
             )

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
    def show_meta(cls, sample='.'):
        lines = []
        # fileformat
        lines.append('##fileformat=VCFv{}'.format(cls.file_format))

        # info
        lines.append(cls.show_info_format('INFO', cls.info))

        # filter

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
        for key, i in zip(cls.info[1:], (2,8)):
            pairs.append('{}={}'.format(key[0], variant[i]))

        return ';'.join(pairs)

    @classmethod
    def extract_variant_gt(cls, variant, gt, alt_motifs):
        vals = {}

        # 5 = coverage
        gts = sorted(list(set([a[9] for a in variant[3] if a[9] is not None and a[-2] == 'full'])))
        if gts:
            if len(gts) == 1:
                gts.append(gts[0])
            vals['GT'] = '/'.join(map(str, gts))
            vals['AL'] = '/'.join([str(a[0]) for a in gt])
            vals['AD'] = '/'.join([str(a[1]) for a in gt])

        if variant[5]:
            vals['DP'] = str(variant[5])

        if alt_motifs:
            vals['AM'] = alt_motifs

        format = []
        for fmt in cls.format:
            if fmt[0] in vals:
                format.append(fmt[0])

        return [':'.join(format), ':'.join([vals[field[0]] for field in cls.format if field[0] in vals])]
