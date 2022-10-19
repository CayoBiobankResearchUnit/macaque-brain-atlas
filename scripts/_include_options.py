#!/usr/bin/env python
# 🐛

import sys
sys.path.insert(0,'/home/klchiou/.local/lib/python3.7/site-packages')

# wget -O data/Macaca_mulatta.Mmul_10.105.gtf.gz http://ftp.ensembl.org/pub/release-105/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.105.gtf.gz
gtf_file = 'data/Macaca_mulatta.Mmul_10.105.gtf'

# wget -O data/Macaca_mulatta.Mmul_10.105.gff.gz http://ftp.ensembl.org/pub/release-105/gff3/macaca_mulatta/Macaca_mulatta.Mmul_10.105.gff3.gz
gff_file = 'data/Macaca_mulatta.Mmul_10.105.gff'

fai_file = 'data/Macaca_mulatta.Mmul_10.dna.toplevel.fa.fai'

# minimum number of features for a cell
min_features = 1000

# minimum number of cells sharing a feature
min_cells = 5

# Number of principal components
n_pcs = 50

glue_n_neighbors = 30