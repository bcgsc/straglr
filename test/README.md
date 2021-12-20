# Test data

## Genotype mode (locus given)
```
straglr.py test.bam /path/to/hg38.fa <output_prefix> --loci test.bed
```
test locus BED (`test.bed`):
```
chr22	45795354	45795424	ATTCT
```
expected output: `genotype.tsv`, `genotype.bed`
## Genome-scan mode
```
straglr.py test.bam /path/to/hg38.fa <output_prefix>
```
expected output: `genome_scan.tsv`, `genome_scan.bed`
