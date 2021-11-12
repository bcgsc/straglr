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
