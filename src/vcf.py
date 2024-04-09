class VCF:
    file_format = '4.1'
    header = 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'
    info = (
            ('END', 1, 'Integer', 'End position of repeat'), 
            ('MOTIF', 1, 'String', 'Motif')
            )
    format = (
              ('GT', 1, 'String', 'Genotype'),
              ('DP', 1, 'Integer', 'Read depth'),
              ('AL', '.', 'String', 'Allelic length(s)'),
              ('AD', '.', 'String', 'Allelic depth(s)'),
              ('AM', '.', 'String', 'Alternate motif(s)'),
              ('AMAD', '.', 'String', 'Alterate motif allelic depth(s)'),
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
    def create_variant_format(cls, variant):
        pairs = []

        # 2 == chr, 4 == repeat
        for key, i in zip(cls.info, (2,4)):
            pairs.append('{}={}'.format(key[0], variant[i]))

        format_str = ':'.join([f[0] for f in cls.format])
        return '{} {}'.format(';'.join(pairs), format_str)

    @classmethod
    def extract_variant_gt(cls, variant, gt, alt_motifs):
        cols = []
        vals = {}
        for field in cls.format:
            vals[field[0]] = '.'
        
        # 5 = coverage
        vals['DP'] = str(variant[5])

        vals['AL'] = '/'.join([str(a[0]) for a in gt])
        vals['AD'] = '/'.join([str(a[1]) for a in gt])

        vals['AM'] = alt_motifs[0]
        vals['AMAD'] = alt_motifs[1]

        return ':'.join([vals[field[0]] for field in cls.format])


        
