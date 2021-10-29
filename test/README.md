# Test data

## Genotype mode (locus given)
```
straglr.py test.bam /path/to/hg38.fa <output.tg.tsv> --loci test.bed
```
test locus BED (`test.bed`):
```
chr22	45795354	45795424	ATTCT
```
expected output: `test.gt.tsv`
## Genome-scan mode
```
straglr.py test.bam /path/to/hg38.fa <output.gs.tsv>
```
expected output: `test.gs.tsv`
