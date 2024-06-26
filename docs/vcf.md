## VCF output from Straglr
- [VCFv4.2 specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf) followed
- Sample name (`--sample`) must be provided (otherwise `.` will be used as sample name)
- One TR sequence from support reads will be picked - closest in size (bp) to reported (mean) value - to represent each non-reference allele in `ALT`

## FILTER
Loci requested for genotyping (specified in `--loci`) may fail and one of the filters (besides `PASS`) was assigned
| Name | Description |
| ------ | ------ |
| PASS | All filters passed |
| UNMATCHED_MOTIF | motif detected is not the same as the specified motif |
| PARTIAL_SPAN | repeat detected not spanning majority (90%) locus |
| INVALID_MOTIF | motif detected outside motif size range |

## INFO
| Field | Description |
| ------ | ------ |
| LOCUS | Locus ID: fourth column (after chrom, statrt, end, motif) in BED file from parameter `--loci` e.g. gene name |
| END | End position of repeat |
| MOTIF | Reference motif: motif detected by TRF on reference sequence |
| COPIES | Reference copies: extracted from TRF result on reference sequence |

## FORMAT
| Field | Description |
| ------ | ------ |
| GT | Genotype: reference (`0`) assigned when reference size (bp) is within size range (10% added to boundaries) of allele closest in size to reference |
| DP | Read depth: number of reads spanning locus (= `coverage` column in TSV output) |
| AL | Allelic lengths: length (bp) of each allele separated by "`/`" |
| ALR | Allelic length ranges: length range (min-max) of each allele separated by "`/`" |
| AC | Allelic copies: copies of each allele separated by "`/`" |
| ACR | Allelic copy ranges: copies range (min-max) of each allele separated by "`/`" | |
| AD | Allelic depths: number of support reads for each allele separated by "`/`" |
| ALT_MOTIF | Alternate motif(s): alternate motif(s) (followed by number of reads in brackets) observed for each allele separated by "`/`" |

