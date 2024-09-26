## VCF output from Straglr
- [VCFv4.5 specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf) followed
- Sample name (`--sample`) must be provided (otherwise `.` will be used as sample name)
- `<CNV:TR>` will be used for all ATL alleles, no `<CNV:TR>` will be reported for REF alleles
- `RUS` reported for STRs (motif length 2-6bp) and `RUL` for VNTRs, following v4.5 recommendation
- `RB` (allele length) will be reported if `--genotype_in_size` was specified, otherwise `RUC` (copy number)  will be reported
- `SVLEN` equals length (bp) of the reference allele, not the difference between start and end coordinate
- `RUS_REF` and `RUL_REF` added to `INFO` fields to indicate motif sequence (for STRs) and motif lengths (for VNTRs), as variants may exhibit different motifs
- `CIRUB` and `CIRB` were based on minimum and maximum values detected in all support reads for each allele (reported in TSV file)
- currently `RN` is always 1 as complex TRs are not supported

## FILTER
Loci requested for genotyping (specified in `--loci`) may fail and one of the filters (besides `PASS`) was assigned
| Name | Description |
| ------ | ------ |
| PASS | All filters passed |
| UNMATCHED_MOTIF | motif detected is not the same as the specified motif |
| PARTIAL_SPAN | repeat detected not spanning majority (90%) locus |
| INVALID_MOTIF | motif detected outside motif size range |
| CLUSTERING_FAILED | clustering failed to formulate genotype |

## INFO
| Field | Description |
| ------ | ------ |
| LOCUS | Locus ID: fourth column (after chrom, statrt, end, motif) in BED file from parameter `--loci` e.g. gene name |
| RUS_REF | Motif sequence detected by TRF on reference sequence if STR reported |
| RUL_REF | Motif length (bp) detected by TRF on reference sequence if VNTR reported |
| SVLEN | Allele length of reference sequence as determined by TRF |
| RN | Number of motifs in each allele, currently always 1 as complex motifs are not supported |
| RUS | Motif sequences of each allele, currently alway one sequence for each allele |
| RUL | Length of motif sequences reported in RUS |
| RUC | Copy number of each motif in RUS, reported when `--genotype_in_size` is not used |
| RB | Number of bases for each motif in RUS, reported when `--genotype_in_size` is used |
| CIRUC | Differences between lower and upper bounds in all support reads from RUC in each RUS |
| CIRB | Differences between lower and upper bounds in all support reads from RB in each RUS |

## FORMAT
| Field | Description |
| ------ | ------ |
| GT | Genotype: reference (`0`) assigned when reference size (determined by TRF, not SVLEN) is within size range (10% added to boundaries) of allele closest in size to reference |
| DP | Read depth: number of reads spanning locus (= `coverage` column in TSV output) |
| AD | Allelic depths: number of support reads for each allele |

