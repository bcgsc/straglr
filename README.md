# Straglr - *S*hort-*t*andem *r*epe*a*t *g*enotyping using *l*ong *r*eads

Straglr is a tool that can be used for genome-wide scans for tandem repeat(TR) expansions or targeted genotyping using long-read alignments.

## Installation
Straglr is implemented in Python 3.8 and has been tested in Linux environment.

Straglr depends on [Tandem Repeat Finder(TRF)](https://tandem.bu.edu/trf/trf.html) for identifying TRs and [blastn](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) for motif matching. (TRF and blastn executables must be in `$PATH`). Other Python dependencies are listed in `requirements.txt`.

The file `environment.yaml` can by used by conda to create an environment with all dependencies installed:
```
conda env create --name straglr --file=environment.yaml
```
Straglr can be added to the environment via `pip`,
```
pip install git+https://github.com/bcgsc/straglr.git@v1.3.0#egg=straglr
```
(for example to install v1.3.0), or run directly from the cloned repository:
```
conda activate straglr
./straglr.py
```

## Input
Long read alignments sorted by genomic coordinates in BAM format against the reference genome. Suggested aligner: [Minimap2](https://github.com/lh3/minimap2) **-- Please use the option `-Y` to enable soft-clipping so that read sequences can be assessed directly from the BAM file.** 

## Usage
```
python straglr.py <mm2.bam> <reference_fasta> <output_prefix> [--loci loci.bed] [--exclude skip_regions.bed] [--chroms chr] [--regions regions.bed] [--min_support N] [--min_ins_size N] [--min_str_len N] [--max_str_len N] [--nprocs N] [--genotype_in_size] [--max_num_clusters N] [--min_cluster_size N] [--working_dir] [--tmpdir] [--debug]
```

Some common parameters:

`--loci`: a BED file containing loci to be genotyped. 4 column BED format: chromosome start end repeat

`--exclude`: a BED file containing regions to be skipped in genome-scan (e.g. long segmental duplications or pericentromeric regions) 

`--chroms`: space-separated list of specific chromosomes for genome-scan

`--regions`: a BED file containing regions to be used only in genome-scan

`--include_alt_chroms`: include ALT chromosomes (chromosomes with "_" in names)  in genome scan (Default: NOT included)

`--use_unpaired_clips`: include examination of unpaired clipped alignments in genome scan to detect expansion beyond read size (Default:NOT used)

`--min_support`: minimum number of support reads for an expansion to be captured in genome-scan (Default:2)

`--min_ins_size`: minimum increase in size (relative to the reference genome) for an expansion to be captured in genome-scan (Default:100)

`--min_str_len`: minimum length of repeat-motif for an expansion to be captured in genome-scan (Default:2)

`--max_str_len`: maximum length of repeat-motif for an expansion to be captured in genome-scan (Default:50)

`--nprocs`: number of processes to use in Python's multiprocessing (Default:1)

`--genotype_in_size`: report genotype (column 5 of TSV output) in terms of allele sizes instead of copy numbers

`--max_num_clusters`: maximum number of clusters to be tried in Gaussian Mixture Model (GMM) clustering (Default:2)

`--min_cluster_size`: minimum number of reads required to constitute a cluster (allele) in GMM clustering (Default:2)

`--trf_args`: TRF arguments (Default:2 5 5 80 10 10 500)

`--tmpdir`: user-specified directory for holding temporary files

#### Example application: genome scan to detect TRs longer than the reference genome by 100bp:
The most common use of Straglr is for detecting TR expansions over the reference genome by a defined size threshold. This will save computations spent on genotyping the majority of TRs in the human genome with no substantial change in lengths. The identified expanded alleles can then be screened for pathogenicity by comparing against known TR polymorphisms. A sample Straglr run to detect expansions larger than the reference alleles by 100 bp on TR loci 2-100bp in motif length:
```
straglr.py sample.mm2.bam hg38.fa straglr_scan --min_str_len 2 --max_str_len 100 --min_ins_size 100 --genotype_in_size --exclude hg38.exclude.bed --min_support 2 --max_num_clusters 2 --nprocs 32
```
Highly repetitive genomic regions may be problematic for aligners and give rise to questionable genotyping results. They can be skipped over in Straglr's genome scan. To generate a bed file that contains all segmental duplications, centromeric and gap regions for exclusion from a Straglr run:
```
(cut -f1-3 hg38.segdups.bed;awk '$3-$2>=10000' hg38.simple_repeats.bed | cut -f1-3;cat hg38.centromeres.bed hg38_gaps.bed) | bedtools sort -i - | bedtools merge -i - -d 1000 > hg38.exclude.bed
```

#### Example application: genome-wide genotyping
1. Download UCSC Simple Repeats track in output format `all fields from selected table` using the online `Table Browser` tool
2. Convert downloaded table into bed format, skipping all homopolymers, with last field specifying motif sequence, e.g.:
	```
	grep -v '^#' simple_repeats.downloaded.tsv | awk -vOFS='\t' 'length($17)>1 {print $2,$3,$4,$17}' > simple_repeats.bed
	```
3. Split whole-genome `bed` file into batches with smaller numbers of loci (e.g. 10,000), e.g.:
	```
	split -l 10000 -d -a 4 --additional-suffix=.bed simple_repeats.bed batch
	```
4. Run Straglr on Minimap2 alignments for each batch of TR loci in parallel on, for example, computing cluster, e.g.:
	```
	straglr.py sample.mm2.bam hg38.fa batch0001 --genotype_in_size --min_support 2 --loci batch0001.bed --max_str_len 100 --max_num_clusters 2 --nprocs 32
	```

## Output
1. \<output_prefix>.tsv - detailed output one support read per line 
	* chrom - chromosome name
	* start - start coordinate of locus
	* end - end coordinate of locus
	* target_repeat - consensus(shortest) repeat motif from genome scan or target motif in genotyping
	* locus - locus in UCSC format (chrom:start-end)
	* coverage - coverage depth of locus
	* genotype - copy numbers (default) or sizes (`--genotype_in_size`) of each allele detected for given locus, separate by semi-colon(";") if multiple alleles detected, with number of support reads in bracket following each allele copy number/size. An example of a heterozygous allele in size: `990.8(10);30.9(10)` (Alleles preceded by `>` indicate minimum values, as full alleles are not captured in any support reads)
	* actual_repeat - actual repeat motif detected in mapped read
	* read_name - mapped read name
	* copy_number - number of copies of repeat in allele
	* size - size of allele
	* read_start - start position of repeat in support read
	* strand - strand of reference genome from which read originates
	* allele - allele to which support read is assigned
	* read_status - classification of mapped read
		* "full": read captures entire repeat (counted as support read)
		* "partial": read does not capture entire repeat (counted as support read)
		* "skipped": "not_spanning" - read does not span across locus (NOT counted as support read)
		* "failed" - read not used for genotyping (NOT counted as support read). Reasons are indicated with following descriptors:
		
		| Descriptor | Explanation |
		| --- | :--- |
		| cannot_extract_sequence | cannot extract repeat sequence, could be because the repeat is deleted for the read in question, or regions flanking motif are deleted |
		| motif_size_out_of_range | motif size detected outside specified size range |
		| insufficient_repeat_coverage | repeat detected does not cover enough (50%) of expansion/insertion sequence |
		| partial_and_insufficient_span | repeat not covering enough (90%) query minus flanking sequences |
		| unmatched_motif | no repeat found matching target motif |

2. \<output_prefix>.bed - summarized genotypes one locus per line
	* chrom - chromosome name
	* start - start coordinate of locus
	* end - end coordinate of locus
	* repeat_unit - repeat motifi
	* allele\<N>.size, where N={1,2,3...} depending on `--max_num_clusters` e.g. N={1,2} if `--max_num_clusters`==2 (default)
	* allele\<N>.copy_number
	* allele\<N>.support

3.  [\<output_prefix>.vcf](docs/vcf.md) 

## Utilities

- `straglr_compare.py`: [compare Straglr's results](docs/straglr_compare.md)
- `extract_repeats.py`: extract repeat sequences from alignment BAM given Straglr's TSV output

## Contact
[Readman Chiu](mailto:rchiu@bcgsc.ca)

## Citation
Chiu R, Rajan-Babu IS, Friedman JM, Birol I. Straglr: discovering and genotyping tandem repeat expansions using whole genome long-read sequences. *Genome Biol* 22, 224 (2021). https://doi.org/10.1186/s13059-021-02447-3

