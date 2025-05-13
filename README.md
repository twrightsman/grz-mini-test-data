# Tiny synthetic test datasets for genomDE Model Project GRZs

## Environment

```shell
pixi shell
./build.sh
```

## Genomes

Each genome is three "chromosomes" with a single gene centered in each.


## Panel

- single, tumor-only, GRCh37
- 225x paired-end run of gene2 and gene3


## WGS

- trio, germline-only (rare disease), GRCh37
- 50X index paired-end in two 25X runs
- 40X mother single-end run
- 20X father paired-end
  - should fail because under threshold


## WGS long-read

- single, tumor+germline, GRCh38
- 110X Nanopore tumor in multiple FASTQs
- 40X PacBio HiFi germline in multiple BAMs
