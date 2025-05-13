#!/usr/bin/env bash

set -euxo pipefail


# Clean
rm -rf submissions/*/{files,logs}


# Genomes
mkdir -p references/{GRCh37,GRCh38}
## GRCh37
# Credit: https://stackoverflow.com/a/69982328
awk '$1 == "chr1"' < hg19_439_omim_genes.bed | shuf --random-source=<(yes 42) | head -n3 > references/GRCh37/test_regions.bed
GRCh37_OUT="references/GRCh38/chr1.fa.gz"
if [ ! -f "$GRCh37_OUT" ]; then
  curl https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.1.fa.gz | gzip -cd | sed 's/^>1/>chr1/' | bgzip -@4 -c > "$GRCh37_OUT"
fi
samtools faidx references/GRCh37/chr1.fa.gz
bedtools getfasta -fi references/GRCh37/chr1.fa.gz -bed <(bedtools slop -b 500 -g references/GRCh37/chr1.fa.gz.fai -i references/GRCh37/test_regions.bed) | fold -w 60 | awk 'BEGIN{i=1} {if ($0 ~ /^>chr/) {print ">chr" i "_GRCh37"; i+=1} else {print}}' | bgzip -@4 -c > references/GRCh37/mini.fa.gz
samtools faidx references/GRCh37/mini.fa.gz
awk -v OFS=$'\t' '{print $1, $2 - $2 + 500, $3 - $2 + 500}' < references/GRCh37/test_regions.bed | awk 'BEGIN{i=1} {sub(/^chr[0-9]+/, "chr" i "_GRCh37"); i+=1; print}' > references/GRCh37/mini.genes.bed
## GRCh38
awk '$1 == "chr1"' < hg38_440_omim_genes.bed | shuf --random-source=<(yes 42) | head -n3 > references/GRCh38/test_regions.bed
GRCh38_OUT="references/GRCh38/chr1.fa.gz"
if [ ! -f "$GRCh38_OUT" ]; then
  curl https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz | gzip -cd | sed 's/^>1/>chr1/' | bgzip -@4 -c > "$GRCh38_OUT"
fi
samtools faidx references/GRCh38/chr1.fa.gz
bedtools getfasta -fi references/GRCh38/chr1.fa.gz -bed <(bedtools slop -b 500 -g references/GRCh38/chr1.fa.gz.fai -i references/GRCh38/test_regions.bed) | fold -w 60 | awk 'BEGIN{i=1} {if ($0 ~ /^>chr/) {print ">chr" i "_GRCh38"; i+=1} else {print}}' | bgzip -@4 -c > references/GRCh38/mini.fa.gz
samtools faidx references/GRCh38/mini.fa.gz
awk -v OFS=$'\t' '{print $1, $2 - $2 + 500, $3 - $2 + 500}' < references/GRCh38/test_regions.bed | awk 'BEGIN{i=1} {sub(/^chr[0-9]+/, "chr" i "_GRCh38"); i+=1; print}' > references/GRCh38/mini.genes.bed


# Panel
mkdir -p samples/panel
simuG -seed 42 -refseq references/GRCh37/mini.fa.gz -snp_count 20 -indel_count 10 -prefix samples/panel/tumor
tail -n2 < references/GRCh37/mini.genes.bed > samples/panel/targets.bed
bedtools getfasta -fi samples/panel/tumor.simseq.genome.fa -bed samples/panel/targets.bed | fold -w 60 > samples/panel/targets.fa
art_illumina --rndSeed 42 --seqSys HS25 --in samples/panel/targets.fa --paired --len 150 --fcov 225 --mflen 500 --sdev 50 --out samples/panel/tumor_
gzip -nf samples/panel/*.fq
mkdir submissions/panel/files
mv samples/panel/{targets.bed,*.fq.gz} submissions/panel/files/
gzip -n submissions/panel/files/targets.bed
touch submissions/panel/files/panel.vcf.gz


# WGS
mkdir -p samples/wgs
simuG -seed 24-index -refSeq references/GRCh37/mini.fa.gz -snp_count 10 -indel_count 5 -prefix samples/wgs/index
simuG -seed 24-mother -refSeq references/GRCh37/mini.fa.gz -snp_count 10 -indel_count 5 -prefix samples/wgs/mother
simuG -seed 24-father -refSeq references/GRCh37/mini.fa.gz -snp_count 10 -indel_count 5 -prefix samples/wgs/father
art_illumina --rndSeed 24-index --seqSys HS25 --in samples/wgs/index.simseq.genome.fa --paired --len 150 --fcov 25 --mflen 500 --sdev 50 --out samples/wgs/index_1_
art_illumina --rndSeed 24-index-second --seqSys HS25 --in samples/wgs/index.simseq.genome.fa --paired --len 150 --fcov 25 --mflen 500 --sdev 50 --out samples/wgs/index_2_
art_illumina --rndSeed 24-mother --seqSys HS25 --in samples/wgs/mother.simseq.genome.fa --len 150 --fcov 40 --out samples/wgs/mother
art_illumina --rndSeed 24-father --seqSys HS25 --in samples/wgs/father.simseq.genome.fa --paired --len 150 --fcov 20 --mflen 500 --sdev 50 --out samples/wgs/father_
gzip -nf samples/wgs/*.fq
mkdir submissions/wgs/files
mv samples/wgs/*.fq.gz submissions/wgs/files/
touch submissions/wgs/files/{index,mother,father}.vcf.gz


# WGS long-read
mkdir -p samples/wgs_lr
simuG -seed 24-germline-lr -refSeq references/GRCh38/mini.fa.gz -snp_count 10 -indel_count 5 -prefix samples/wgs_lr/germline
simuG -seed 24-tumor-lr -refSeq samples/wgs_lr/germline.simseq.genome.fa -snp_count 30 -indel_count 10 -prefix samples/wgs_lr/tumor
pbsim --strategy wgs --seed 24-tumor-lr-nano --prefix samples/wgs_lr/tumor --genome samples/wgs_lr/tumor.simseq.genome.fa --depth 110 --method qshmm --qshmm "${PIXI_PROJECT_ROOT}/.pixi/envs/default/data/QSHMM-ONT-HQ.model" --accuracy-mean 0.90 --difference-ratio '39:24:36'
pbsim --strategy wgs --seed 24-germline-lr-pr --prefix samples/wgs_lr/germline --genome samples/wgs_lr/germline.simseq.genome.fa --depth 40 --method errhmm --errhmm "${PIXI_PROJECT_ROOT}/.pixi/envs/default/data/ERRHMM-SEQUEL.model" --pass-num 10
ccs samples/wgs_lr/germline_0001.bam samples/wgs_lr/germline_0001.hifi.bam
ccs samples/wgs_lr/germline_0002.bam samples/wgs_lr/germline_0002.hifi.bam
ccs samples/wgs_lr/germline_0003.bam samples/wgs_lr/germline_0003.hifi.bam
mkdir submissions/wgs_lr/files
mv samples/wgs_lr/{*.hifi.bam,*.fq.gz} submissions/wgs_lr/files/
touch submissions/wgs_lr/files/{tumor,germline}.vcf.gz
