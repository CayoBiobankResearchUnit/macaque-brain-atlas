#!/bin/bash

mkdir -p datasets

prefix=cortex-dev

wget -O datasets/${prefix}-expr.tsv.gz https://cells.ucsc.edu/cortex-dev/exprMatrix.tsv.gz
wget -O datasets/${prefix}-meta.tsv https://cells.ucsc.edu/cortex-dev/meta.tsv

prefix=human-cortex

wget -O datasets/${prefix}-expr.tsv.gz https://cells.ucsc.edu/allen-celltypes/human-cortex/various-cortical-areas/exprMatrix.tsv.gz
wget -O datasets/${prefix}-meta.tsv https://cells.ucsc.edu/allen-celltypes/human-cortex/various-cortical-areas/meta.tsv

prefix=dev-brain-regions

wget -O datasets/${prefix}-expr.tsv.gz https://cells.ucsc.edu/dev-brain-regions/wholebrain/exprMatrix.tsv.gz
wget -O datasets/${prefix}-meta.tsv https://cells.ucsc.edu/dev-brain-regions/wholebrain/meta.tsv

prefix=adult-brain-vasc-peri

wget -O datasets/${prefix}-expr.tsv.gz https://cells.ucsc.edu/adult-brain-vasc/perivascular/exprMatrix.tsv.gz
wget -O datasets/${prefix}-meta.tsv https://cells.ucsc.edu/adult-brain-vasc/perivascular/meta.tsv

prefix=adult-brain-vasc-endo

wget -O datasets/${prefix}-expr.tsv.gz https://cells.ucsc.edu/adult-brain-vasc/endothelial/exprMatrix.tsv.gz
wget -O datasets/${prefix}-meta.tsv https://cells.ucsc.edu/adult-brain-vasc/endothelial/meta.tsv

prefix=vascular-dev

wget -O datasets/${prefix}-expr.mtx.gz https://cells.ucsc.edu/vascular-dev/matrix.mtx.gz
wget -O datasets/${prefix}-features.tsv.gz https://cells.ucsc.edu/vascular-dev/features.tsv.gz
wget -O datasets/${prefix}-barcodes.tsv.gz https://cells.ucsc.edu/vascular-dev/barcodes.tsv.gz
wget -O datasets/${prefix}-meta.tsv https://cells.ucsc.edu/vascular-dev/meta.tsv

prefix=early-brain

wget -O datasets/${prefix}-expr.tsv.gz https://cells.ucsc.edu/early-brain/exprMatrix.tsv.gz
wget -O datasets/${prefix}-meta.tsv https://cells.ucsc.edu/early-brain/meta.tsv

prefix=brain-vasc-atlas

wget -O datasets/${prefix}-expr.mtx.gz https://cells.ucsc.edu/brain-vasc-atlas/matrix.mtx.gz
wget -O datasets/${prefix}-features.tsv.gz https://cells.ucsc.edu/brain-vasc-atlas/features.tsv.gz
wget -O datasets/${prefix}-barcodes.tsv.gz https://cells.ucsc.edu/brain-vasc-atlas/barcodes.tsv.gz
wget -O datasets/${prefix}-meta.tsv https://cells.ucsc.edu/brain-vasc-atlas/meta.tsv

prefix=myeloid-neuroinflam

wget -O datasets/${prefix}-expr.tsv.gz https://cells.ucsc.edu/myeloid-neuroinflam/all/exprMatrix.tsv.gz
wget -O datasets/${prefix}-meta.tsv https://cells.ucsc.edu/myeloid-neuroinflam/all/meta.tsv

prefix=mouse-drg-injury

wget -O datasets/${prefix}-expr.tsv.gz https://cells.ucsc.edu/mouse-drg-injury/peripheral-central-injury/exprMatrix.tsv.gz
wget -O datasets/${prefix}-meta.tsv https://cells.ucsc.edu/mouse-drg-injury/peripheral-central-injury/meta.tsv

prefix=dev-inhibitory-neurons-macaque

wget -O datasets/${prefix}-expr.tsv.gz https://cells.ucsc.edu/dev-inhibitory-neurons/macaque/exprMatrix.tsv.gz
wget -O datasets/${prefix}-meta.tsv https://cells.ucsc.edu/dev-inhibitory-neurons/macaque/meta.tsv



# Dryad

wget -O datasets/dryad/H18.06.006.MTG.4000.expand.rep1.features.csv https://datadryad.org/stash/downloads/file_stream/1641043
wget -O datasets/dryad/H18.06.006.MTG.4000.expand.rep1.genes.csv https://datadryad.org/stash/downloads/file_stream/1640993
wget -O datasets/dryad/H18.06.006.MTG.4000.expand.rep1.matrix.csv https://datadryad.org/stash/downloads/file_stream/1640994

wget -O datasets/dryad/H18.06.006.MTG.4000.expand.rep2.features.csv https://datadryad.org/stash/downloads/file_stream/1641044
wget -O datasets/dryad/H18.06.006.MTG.4000.expand.rep2.genes.csv https://datadryad.org/stash/downloads/file_stream/1640996
wget -O datasets/dryad/H18.06.006.MTG.4000.expand.rep2.matrix.csv https://datadryad.org/stash/downloads/file_stream/1640997

wget -O datasets/dryad/H18.06.006.MTG.4000.expand.rep3.features.csv https://datadryad.org/stash/downloads/file_stream/1641045
wget -O datasets/dryad/H18.06.006.MTG.4000.expand.rep3.genes.csv https://datadryad.org/stash/downloads/file_stream/1640999
wget -O datasets/dryad/H18.06.006.MTG.4000.expand.rep3.matrix.csv https://datadryad.org/stash/downloads/file_stream/1641000

wget -O datasets/dryad/H19.30.001.STG.4000.expand.rep1.features.csv https://datadryad.org/stash/downloads/file_stream/1641047
wget -O datasets/dryad/H19.30.001.STG.4000.expand.rep1.genes.csv https://datadryad.org/stash/downloads/file_stream/1641005
wget -O datasets/dryad/H19.30.001.STG.4000.expand.rep1.matrix.csv https://datadryad.org/stash/downloads/file_stream/1641006

wget -O datasets/dryad/H19.30.001.STG.4000.expand.rep2.features.csv https://datadryad.org/stash/downloads/file_stream/1641048
wget -O datasets/dryad/H19.30.001.STG.4000.expand.rep2.genes.csv https://datadryad.org/stash/downloads/file_stream/1641008
wget -O datasets/dryad/H19.30.001.STG.4000.expand.rep2.matrix.csv https://datadryad.org/stash/downloads/file_stream/1641009

wget -O datasets/dryad/H20.30.001.STG.4000.expand.rep1.features.csv https://datadryad.org/stash/downloads/file_stream/1641051
wget -O datasets/dryad/H20.30.001.STG.4000.expand.rep1.genes.csv https://datadryad.org/stash/downloads/file_stream/1641017
wget -O datasets/dryad/H20.30.001.STG.4000.expand.rep1.matrix.csv https://datadryad.org/stash/downloads/file_stream/1641018

wget -O datasets/dryad/H20.30.001.STG.4000.expand.rep2.features.csv https://datadryad.org/stash/downloads/file_stream/1641052
wget -O datasets/dryad/H20.30.001.STG.4000.expand.rep2.genes.csv https://datadryad.org/stash/downloads/file_stream/1641020
wget -O datasets/dryad/H20.30.001.STG.4000.expand.rep2.matrix.csv https://datadryad.org/stash/downloads/file_stream/1641021

wget -O datasets/dryad/H20.30.001.STG.4000.expand.rep3.features.csv https://datadryad.org/stash/downloads/file_stream/1641053
wget -O datasets/dryad/H20.30.001.STG.4000.expand.rep3.genes.csv https://datadryad.org/stash/downloads/file_stream/1641023
wget -O datasets/dryad/H20.30.001.STG.4000.expand.rep3.matrix.csv https://datadryad.org/stash/downloads/file_stream/1641024

wget -O datasets/dryad/H22.26.401.MTG.4000.expand.rep1.features.csv https://datadryad.org/stash/downloads/file_stream/1641055
wget -O datasets/dryad/H22.26.401.MTG.4000.expand.rep1.genes.csv https://datadryad.org/stash/downloads/file_stream/1641029
wget -O datasets/dryad/H22.26.401.MTG.4000.expand.rep1.matrix.csv https://datadryad.org/stash/downloads/file_stream/1641030

wget -O datasets/dryad/H22.26.401.MTG.4000.expand.rep2.features.csv https://datadryad.org/stash/downloads/file_stream/1641056
wget -O datasets/dryad/H22.26.401.MTG.4000.expand.rep2.genes.csv https://datadryad.org/stash/downloads/file_stream/1641032
wget -O datasets/dryad/H22.26.401.MTG.4000.expand.rep2.matrix.csv https://datadryad.org/stash/downloads/file_stream/1641033

for i in H18.06.006.MTG.4000.expand.rep1 H18.06.006.MTG.4000.expand.rep2 H18.06.006.MTG.4000.expand.rep3 H19.30.001.STG.4000.expand.rep1 H19.30.001.STG.4000.expand.rep2 H20.30.001.STG.4000.expand.rep1 H20.30.001.STG.4000.expand.rep2 H20.30.001.STG.4000.expand.rep3 H22.26.401.MTG.4000.expand.rep1 H22.26.401.MTG.4000.expand.rep2; do
echo $i
n_nonzero=$(tail -n+2 datasets/dryad/${i}.matrix.csv | wc -l)
n_genes=$(tail -n+2 datasets/dryad/${i}.genes.csv | wc -l)
n_cells=$(tail -n+2 datasets/dryad/${i}.features.csv | wc -l)
cat <(echo '%%MatrixMarket matrix coordinate integer general') <(echo ${n_genes} ${n_cells} ${n_nonzero}) <(tail -n+2 datasets/dryad/${i}.matrix.csv | sed 's/,/ /g') | gzip -c > datasets/dryad/${i}.matrix.mtx.gz
done