v1.0.0
- first public release

v1.1.0
- added `--simple` option to report genotype result for each repeat locus on a single line as an alternative output format
- fixed `setup.py` and `requirements.txt` to make sure `pip install` using Github URL works

v1.1.1
- bugfix in `tre.py`

v1.2.0
- bugfix: forgot to pass `reads_fasta` to `extract_missed_clipped()`
- bugfix: `min_support` checking should use "greater than or equal to"
- removed dependency on `intspan`
- reduced `min_len` from 20 to 15 for comparing motifs using `blastn`
- bugfix: keep track of genomic coordinates of all alignments to detect reads that truly span long reference repeats
- assigned `min_span` explicitly in conditional instead of relying on initialized value - `min_span` does not change in some cases when runing genome scan for some unknown reason
- added `check_same_pats()` in `is_same_repeat()` to check if one motif is subsequence of another in `blastn` matches
- test data added
- output BED file by default in addition to TSV (user provide output prefix instead of name). BED output simplied without details such as read names.
- added `--regions` to specify (in BED format) specific regions for scanning
- added `--min_cluster_size` (default=2) to separate from `--min_support` so that alleles with fewer read support can be segregated
- remove optional usage of dbscan for clustering
- improve boundary definition of long(kb range) repeat loci in genome scan
- added `--tmpdir` to allow user to specify TEMP location for generating tmp files
- allowed `--min_cluster_size` to be bigger than `--min_suport` for genotyping mode only

v1.3.0
- `strand` column added to `.tsv` output to indicated strand of support read, `read_start` is recalculated based on read instead of alignment coordinate for negative strand sequences
- `blastn` jobs deployed per locus for compareing target and detected motifs combined to one job per process
- minimum coverarge, percent identity, and e-value added to `blastn` jobs when launched
- when searching for repeat in insertion site, `target_flank` increased from 2000 to 3000
- `bedtools` merge parameter `d` for merging detected insertion loci increased from 50 to 100
- reads identified as `clipped` at repeat locus will not be examined again for same locus as if it is a full alignment
- complete removal of temporary bed files unless `--debug` specified
- `straglr_compare.py` added for comparing straglr results from test against control samples
- `Python`3.8 now needed becasuse `Scipy`1.8 is needed for t-test in `straglr_compare.py`
- new parameter `trf_args` to allow user specify his own preferred TRF settings
- changed seeding method of Gaussian Mixture Model to "k-means++" to reduce CPU overhead
- `--working_directory` added to allow user specify working directory, otherwise current directory will be used

v1.4.0
- `locus`, `coverage`, `actual_repeat`, `read_status` columns added, `repeat_unit` renamed as `target_repeat` in `.tsv` output 
- report failed reads in `.tsv` output with reason given in `read_status` column
- added reporting of unpaired clipped alignments in genotyping results
- allow unpaired clipped alignments to initiate TRE search in genome scan
- added `--include_alt_chroms` to include ATL chromosomes in genome scan; added code to ensure chromosome boundaries are not violated
- get rid of `--working_directory`, all temporary files will be found in `$TMPDIR` (or location specified by `--tmpdir`) if `--debug` is turned on

v1.4.1
- detection and reporting of partial repeats in genotyping now an option (`--include_partials`)
- bugfix: forgot to add read status to output vector from regex (not TRF) examination of repeats
- bugfix: assign partials to allele if size is between maximum and minimum read sizes, do not assign if size cannot be placed
- changed `closeness_to_end` from 10 to 1kb for the "unclipped" end of clipped alignments
- always go ahead to genotyping phase even if repeat in insertion cannot be matched against locus in reference genome (behaviour before `1.4.0`)
- skip read if it has multiple clipped alignments to same end
- fixed `test.bam` (duplicate records for some reads) and updated test results

v1.5.0
- output VCF file
- added `--sex` parameter: when `m` is specified, clustering of alleles is limited to one cluster for chrX and genotype is automatically homozygous in VCF for chrX loci
- support CRAM file

v1.5.1
- option to put `-` to NOT specify intended motif (4th column in BED file) in genotype mode, detected motif will not be screened for agreement
- `--symbolic` parameter to report ALT alleles in symbolic form - only enforced when motif of alternative alleles is the same as reference
- `CLUSTERING_FAILED` added as (default) reason for failure to detect/report genotype
- ALT assignment conditions changed to using interquartile range instead of all support read sizes and at least one copy number difference from REF
- bugfix: alternative motifs checked for equivalence for proper count reports in VCF output
- bugfix: size and copy range of alleles only include "full" supporting reads even if partials were included, but "partial" reads will be reported in support counts
- bugfix: allele motif equivalence checked to perform proper tallying for determining consensus motif in variant
- used read name as tiebreaker for choosing support read sequence as ALT
- add extensions (`.fa`, `.blastn`, `.bed`) to temporary files to help debugging
- `MOTIF` renamed to `RU`, `COPIES` to `REF` in VCF fields

v1.5.2
- post-processing of GMM clustering results modified
- VCFv4.5 on tandem repeats (`<CNV:TR>`) followed

v1.5.3
- fixed VCF bugs
	1. INFO field separator should be "=" instead of ":"
	2. AD separator should be "," instead of "/"
	3. CIRC definition typo
