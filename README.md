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
python straglr.py <mm2.bam> <reference_fasta> <output.tsv> [--loci loci.bed] [--exclude skip_regions.bed] [--chroms chr] [--min_ins_size N] [--nprocs N]
```

Some common parameters:

`--loci`: a BED file containing loci for which only genotyping is performed. A 4 column BED format: chromosome start end repeat

`--exclude`: a BED file containing regions such as segmental duplications or pericentromeric regions where alignment is less reliable and analysis is preferrably skipped. This can be compiled using UCSC's "Table Browser" tool 

`--chroms`: space-separated list of chromosomes for which results are only obtained

`--min_ins_size`: when used for searching repeat expansions, minimum insertion size for detection (<50 is not desirable as long reads are prone to small indels)

`--nprocs`: number of processes to use in Python's multiprocessing

#### Example application: genome scan to detect TRs longer than the reference genome by 100bp:
The most common use of Straglr is for detecting TR expansions over the reference genome by a defined size threshold. This will save computations spent on genotyping the majority of TRs in the human genome with no significant change in lengths. The identified expanded alleles can then be screened for pathogenicity by comparing against known TR polymorphisms. A sample Straglr run to detect expansions larger than the reference alleles by 100 bp on TR loci 2-100bp in motif length:
```
straglr.py sample.mm2.bam hg38.fa straglr_scan.tsv --min_str_len 2 --max_str_len 100 --min_ins_size 100 --genotype_in_size --exclude hg38.exclude.bed --min_support 2 --max_num_clusters 2 --nprocs 32
```
Highly repetitive genomic regions are likely to be problematic for aligners and give rise to unreliable genotyping results. They can be skipped over in Straglr's genome scan. To generate a bed file that contains all segmental duplications, centromeric and gap regions for exclusion from a Straglr run:
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
straglr.py sample.mm2.bam hg38.fa batch0001.tsv --genotype_in_size --min_support 2 --loci batch0001.bed --max_str_len 100 --max_num_clusters 2 --nprocs 36
```

## Output
Default TSV outputs information of 1 support read per line. Description of columns in output:
* chrom - chromosome name
* start - start coordinate of locus
* end - end coordinate of locus
* repeat_unit - repeat motif sequence
* genotype - copy numbers (default) or sizes (`--genotype_in_size`) of each allele detected for given locus, separate by semi-colon(";") if multiple alleles detected, with number of support reads in bracket following each allele copy number/size. An example of a heterozygyous allele in size: `990.8(10);30.9(10)`
* read - name of support read
* copy_number - number of copies of repeat in allele
* size - size of allele
* read_start - start position of repeat in support read
* allele - allele that support read is assigned to

The first 5 columns can be cut and collapsed (using Unix `uniq`) into essentially a bed file without the detailed read information for downstream analysis. Alternatively, Straglr can be run with `--simple` to generate an output that collapses output into one locus per line, with the first 5 columns same as the default output, and the additional columns:
* reads - comma-separated list of support read names separated by semi-colon for each allele
* copy_numbers: comma-separated list of copy numbers in support reads separated by semi-colon for each allele
* sizes: comma-separated list of allele sizes in support reads separated by semi-colon for each allele
* read_starts: comma-separated list of repeat start coordinates in support reads separated by semi-colon for each allele

## Contact
[Readman Chiu](mailto:rchiu@bcgsc.ca)

## Citation
Chiu R, Rajan-Babu IS, Friedman JM, Birol I. Straglr: discovering and genotyping tandem repeat expansions using whole genome long-read sequences. *Genome Biol* 22, 224 (2021). https://doi.org/10.1186/s13059-021-02447-3

