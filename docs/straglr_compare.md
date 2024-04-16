## Comparing Straglr results

`straglr_compare.py` is provided to compare Straglr results from test sample against parent(s) or normal control(s) for detection of _de novo_ expansions:
* corresponding alleles between test and control are found by intersecting their Straglr `.tsv` outputs using BEDTools. Loci that overlap at least 90% with each other are considered the sa
me
* actual sizes or copy numbers (`actual_repeat` column in `.tsv`) detected in each supporting read for each test allele are compared against their counterparts in the control using a 2-sam
ple t-test
* test alleles deemed to be larger than their counterparts in all control(s) by t-test are reported
* expanded alleles only reported in the test sample but not the control(s) are also reported in the output
* annotation file (GTF) and promoter/enhancer coordinates (BED) can be provided by user so that only expansions overlapping the provided features are reported

### Usage
```
straglr_compare.py <test_straglr.tsv> <control_straglr.tsv> [control_straglr.tsv ...] 
```
| optional parameter  | description  | Default | 
| ---------- | - | - |
| use_copy_number   | use copy number instead of size | False (i.e. size used)
| min_expansion | minimum expansion | 100 |
| min_support | minimum number of support reads in sample event | 4 |
| skip_chroms | skip results from chromosomes specified (e.g. chrY) | None
| pval_cutoff |  p-value cutoff for testing T-test hypothesis | 0.001 |
| gtf | GTF file of annotation | None |
| promoters | BED file of promoters | None |
| enhancers | BED file of enhancers | None |
| old_version | results from original version (v1.0) of Straglr used (`tsv` output different) | False |

### Output
| column | description |
| ------- | - |
|  chrom  | chromosome of test allele |
|  start  | start coordinate of test allele |
|  end  | end coordinate of test allele |
|  repeat  | repeat motif of test allele |
|  ref_size  | reference size (`start` - `end` + 1) |
|  test_allele  | allele size (or copy number) in test sample |
|  test_allele_support  | number of support reads for test allele in sample |
|  expansion  | `test_allele` - `ref_size` + 1 |
|  control_alleles  | control allele size(s)/copy numbers(s) e.g. c1a,c1b;c2a,c2b where c1a and c2a are alleles from control1; c2a,c2b are alleles from control2, etc. "-" if test allele no
t found in control(s) |
|  pvals  | t-test p-val of test allele against every control allele. Output format the same as `control_alleles`. "-" if test allele not found in control(s) |

