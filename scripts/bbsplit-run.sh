#!/bin/bash

# id=NSM345
# id=NSM346
# id=NSM351
# id=NSM332

id=$1
genome1=$2
genome2=$3
genome3=$4


# wget http://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
# gunzip Mus_musculus.GRCm38.dna.toplevel.fa.gz
# mv Mus_musculus.GRCm38.dna.toplevel.fa mouse.fa
# ln data/Homo_sapiens.GRCh38.dna.toplevel.fa human.fa
# ln data/Macaca_mulatta.Mmul_10.dna.toplevel.fa macaque.fa

mkdir -p fastq

if [ ! -f fastq/${id}.fastq.gz ]; then

run_id=$(zcat data/rna_run-data.txt.gz | grep ${id} | cut -f 1)

if [ $run_id = 'Snyder-Mackler_RNA3-019' ]; then
fq_dir=/data/CEM/smacklab/data/bbi/${run_id}-nova_data-rerun/fastqs
else
fq_dir=/data/CEM/smacklab/data/bbi/${run_id}-nova_data/fastqs
fi

cat ${fq_dir}/$(echo ${id} | sed 's/NSM0*/NSM/g')-L00[1-4].fastq.gz > fastq/${id}.fastq.gz
fi

if [ ! -f fastq/${id}-sampled.fastq.gz ]; then
zcat fastq/${id}.fastq.gz | head -n $((10000000 * 4)) | gzip -c > fastq/${id}-sampled.fastq.gz
fi

mkdir -p fastq/bbsplit

bbsplit.sh in=fastq/${id}-sampled.fastq.gz build=1 basename=fastq/bbsplit/${id}_%.fq outu=fastq/bbsplit/${id}_unmapped.fq ambig2=split

# bbsplit.sh in=fastq/${id}-sampled.fastq.gz ref=${genome1}.fa,${genome2}.fa,${genome3}.fa basename=fastq/bbsplit/${id}_%.fq outu=fastq/bbsplit/${id}_unmapped.fq ambig2=split build=$(echo $id | sed 's/NSM0*//g')

# if [ -d ref/index/2 ]; then
# bbsplit.sh in=fastq/${id}.fastq.gz basename=fastq/bbsplit/${id}_%.fq outu=fastq/bbsplit/${id}_unmapped.fq ambig2=split build=2
# else
# bbsplit.sh in=fastq/${id}.fastq.gz ref=${genome1}.fa,${genome2}.fa basename=fastq/bbsplit/${id}_%.fq outu=fastq/bbsplit/${id}_unmapped.fq ambig2=split build=2
# fi