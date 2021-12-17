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
- added enforcement in genotyping such that repeat sequences must extend within 200bp of locus boundary - SVA insertions, which have repeat sequence flanked by non-repeat sequences (including poly-A), will not be retained anymore
- added `--tmpdir` to allow user to specify TEMP location for generating tmp files
