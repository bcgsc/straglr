# Straglr - *S*hort-*t*andem *r*epe*a*t *g*enotyping using *l*ong *r*eads

Straglr is a tool that can be used for genome-wide scans for tandem repeat(TR) expansions or targeted genotyping using long-read alignments.

## Installation
Straglr is implemented in Python 3.7 and has been tested in Linux environment.
Straglr can be installed via pip:

```
pip install git+https://github.com/bcgsc/straglr.git#egg=straglr
```
Straglr depends on [Tandem Repeat Finder(TRF)](https://tandem.bu.edu/trf/trf.html) for identifying TRs and [blastn](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) for motif matching. (TRF and blastn executables must be in `$PATH`)

## Input
Long read alignments sorted by genomic coordindates in BAM format against the reference genome. Suggested aligner: [Minimap2](https://github.com/lh3/minimap2) **-- Please use the option `-Y` to enable soft-clipping so that read sequences can be assessed directly from the BAM file.** 

## Usage
```
python straglr.py <mm2.bam> <reference_fasta> <output_prefix> [--loci loci.bed] [--exclude skip_regions.bed] [--chroms chr] [--regions regions.bed] [--min_support N] [--min_ins_size N] [--min_str_len N] [--max_str_len N] [--nprocs N] [--genotype_in_size] [--max_num_clusters N] [--min_cluster_size N] [--debug]
```

Some common parameters:

`--loci`: a BED file containing loci to be genotyped. 4 column BED format: chromosome start end repeat

`--exclude`: a BED file containing regions to be skipped in genome-scan (e.g. long segmental duplications or pericentromeric regions) 

`--chroms`: space-separated list of specific chromosomes for genome-scan

`--regions`: a BED file containing regions to be used only in genome-scan

`--min_support`: minimum number of suppport reads for an expansion to be captured in genome-scan (Default:2)

`--min_ins_size`: minimum increase in size (relative to the reference genome) for an expansion to be captured in genome-scan (Default:100)

`--min_str_len`: minimum length of repeat-motif for an expansion to be captured in genome-scan (Default:2)

`--max_str_len`: maximum length of repeat-motif for an expansion to be captured in genome-scan (Default:50)

`--nprocs`: number of processes to use in Python's multiprocessing

`--genotype_in_size`: report genotype (column 5 of TSV output) in terms of allele sizes instead of copy numbers

`--max_num_clusters`: maximum number of clusters to be tried in Gausssian Mixture Model (GMM) clustering (Default:2)

`--min_cluster_size`: minimum number of reads required to constitute a cluster (allele) in GMM clustering (Default:2)

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
	* repeat_unit - repeat motif
	* genotype - copy numbers (default) or sizes (`--genotype_in_size`) of each allele detected for given locus, separate by semi-colon(";") if multiple alleles detected, with number of support reads in bracket following each allele copy number/size. An example of a heterozygyous allele in size: `990.8(10);30.9(10)`
	* read - name of support read
	* copy_number - number of copies of repeat in allele
	* size - size of allele
	* read_start - start position of repeat in support read
	* allele - allele that support read is assigned to

2. \<output_prefix>.bed - summarized genotypes one locus per line
	* chrom - chromosome name
	* start - start coordinate of locus
	* end - end coordinate of locus
	* repeat_unit - repeat motifi
	* allele\<N>.size, where N={1,2,3...} depending on `--max_num_clusters` e.g. N={1,2} if `--max_num_clusters`==2 (default)
	* allele\<N>.copy_number
	* allele\<N>.support

## Contact
[Readman Chiu](mailto:rchiu@bcgsc.ca)

## Citation
Chiu R, Rajan-Babu IS, Friedman JM, Birol I. Straglr: discovering and genotyping tandem repeat expansions using whole genome long-read sequences. *Genome Biol* 22, 224 (2021). https://doi.org/10.1186/s13059-021-02447-3

