#!/bin/bash

# id=NSM345
# id=NSM346
# id=NSM351
# id=NSM332

genome1=$1
genome2=$2
genome3=$3


# wget http://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
# gunzip Mus_musculus.GRCm38.dna.toplevel.fa.gz
# mv Mus_musculus.GRCm38.dna.toplevel.fa mouse.fa
# ln data/Homo_sapiens.GRCh38.dna.toplevel.fa human.fa
# ln data/Macaca_mulatta.Mmul_10.dna.toplevel.fa macaque.fa

bbsplit.sh ref=${genome1}.fa,${genome2}.fa,${genome3}.fa build=1

# if [ -d ref/index/2 ]; then
# bbsplit.sh in=fastq/${id}.fastq.gz basename=fastq/bbsplit/${id}_%.fq outu=fastq/bbsplit/${id}_unmapped.fq ambig2=split build=2
# else
# bbsplit.sh in=fastq/${id}.fastq.gz ref=${genome1}.fa,${genome2}.fa basename=fastq/bbsplit/${id}_%.fq outu=fastq/bbsplit/${id}_unmapped.fq ambig2=split build=2
# fi