# Tiny synthetic test datasets for genomDE Model Project GRZs

## Genomes

Each genome is three "chromosomes" with a single gene centered in each.

```
mkdir -p references/{GRCh37,GRCh38}
# Credit: https://stackoverflow.com/a/69982328
awk '$1 == "chr1"' < hg19_439_omim_genes.bed | shuf --random-source=<(yes 42) | head -n3 > references/GRCh37/test_regions.bed
awk '$1 == "chr1"' < hg38_440_omim_genes.bed | shuf --random-source=<(yes 42) | head -n3 > references/GRCh38/test_regions.bed
curl https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.1.fa.gz | gzip -cd | sed 's/^>1/>chr1/' | bgzip -@4 -c > references/GRCh37/chr1.fa.gz
samtools faidx references/GRCh37/chr1.fa.gz
bedtools getfasta -fi references/GRCh37/chr1.fa.gz -bed <(bedtools slop -b 500 -g references/GRCh37/chr1.fa.gz.fai -i references/GRCh37/test_regions.bed) | fold -w 60 | awk 'BEGIN{i=1} {if ($0 ~ /^>chr/) {print ">chr" i "_GRCh37"; i+=1} else {print}}' | bgzip -@4 -c > references/GRCh37/mini.fa.gz
samtools faidx references/GRCh37/mini.fa.gz
awk -v OFS=$'\t' '{print $1, $2 - $2 + 500, $3 - $2 + 500}' < references/GRCh37/test_regions.bed | awk 'BEGIN{i=1} {sub(/^chr[0-9]+/, "chr" i "_GRCh37"); i+=1; print}' > references/GRCh37/mini.genes.bed
curl https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz | gzip -cd | sed 's/^>1/>chr1/' | bgzip -@4 -c > references/GRCh38/chr1.fa.gz
samtools faidx references/GRCh38/chr1.fa.gz
bedtools getfasta -fi references/GRCh38/chr1.fa.gz -bed <(bedtools slop -b 500 -g references/GRCh38/chr1.fa.gz.fai -i references/GRCh38/test_regions.bed) | fold -w 60 | awk 'BEGIN{i=1} {if ($0 ~ /^>chr/) {print ">chr" i "_GRCh38"; i+=1} else {print}}' | bgzip -@4 -c > references/GRCh38/mini.fa.gz
samtools faidx references/GRCh38/mini.fa.gz 
awk -v OFS=$'\t' '{print $1, $2 - $2 + 500, $3 - $2 + 500}' < references/GRCh38/test_regions.bed | awk 'BEGIN{i=1} {sub(/^chr[0-9]+/, "chr" i "_GRCh38"); i+=1; print}' > references/GRCh38/mini.genes.bed
```

## Panel

- single, tumor-only, GRCh37
- 225x paired-end run of gene2 and gene3

```
mkdir -p samples/panel
simuG -seed 42 -refseq references/GRCh37/mini.fa.gz -snp_count 20 -indel_count 10 -prefix samples/panel/tumor
tail -n2 < references/GRCh37/mini.genes.bed > samples/panel/targets.bed
bedtools getfasta -fi samples/panel/tumor.simseq.genome.fa -bed samples/panel/targets.bed | fold -w 60 > samples/panel/targets.fa
art_illumina --rndSeed 42 --seqSys HS25 --in samples/panel/targets.fa --paired --len 150 --fcov 225 --mflen 200 --sdev 10 --out samples/panel/tumor_
gzip samples/panel/*.fq
mkdir submissions/panel/files
mv samples/panel/{targets.bed,*.fq.gz} submissions/panel/files/
gzip submissions/panel/files/targets.bed
touch submissions/panel/files/panel.vcf.gz
```

## WGS

- trio, germline-only (rare disease), GRCH37
- 50X index paired-end in two 25X runs
- 40X mother single-end run
- 20X father paired-end, GRCh38
  - should fail because under threshold

## WGS long-read

- single, tumor+germline, GRCh38
- 110X Nanopore tumor
- 110X PacBio HiFi tumor in two BAM files
- 40X PacBio HiFi germline in single BAM
