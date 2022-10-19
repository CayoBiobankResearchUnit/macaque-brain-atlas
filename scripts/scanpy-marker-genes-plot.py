#!/usr/bin/env python
# 🐛

import sys
import pandas as pd
import scipy as sp, scipy.io as si
import anndata as ad
import scanpy as sc
import scanpy.external as sce
import numpy as np
import gzip as gz
import os
import matplotlib.pyplot as plt
import math
import pyreadr as pyr

exec(open('scripts/_include_options.py').read())
exec(open('scripts/_include_functions.py').read())
sc._settings.settings._vector_friendly=True

arguments = sys.argv
# arguments = ['scanpy-marker-genes.py','rna','rna']

prefix = arguments[1]
analysis = arguments[2]

# Read in main data
adata = ad.read_h5ad('hdf5/anndata_'+prefix+'_all_markergenes.h5ad')
push_status(prefix+' read_h5ad')

# # https://github.com/scverse/scanpy/issues/2239#issuecomment-1104178881
adata.uns['log1p']['base'] = None

# sc.pl.rank_genes_groups(adata,gene_symbols='gene_short_name',n_genes=20,fontsize=8,save='_'+prefix+'_marker_genes_summary.pdf')

# sc.pl.rank_genes_groups(adata,gene_symbols='gene_short_name',n_genes=20,fontsize=8,save='_'+prefix+'_marker_genes_celltype_summary.pdf')

# # https://github.com/scverse/scanpy/issues/2239#issuecomment-1104178881
# adata.uns['log1p']['base'] = None
# sc.tl.filter_rank_genes_groups(adata,min_in_group_fraction=0.1)
# 
# os.makedirs('figures/marker_genes',exist_ok=True)
# sc.pl.rank_genes_groups(adata,gene_symbols='gene_short_name',n_genes=20,fontsize=8,save='_'+prefix+'_marker_genes_celltype_filtered.pdf')
# 
# sc.pl.rank_genes_groups(adata,key='rank_genes_groups_filtered',fontsize=8,save='_'+prefix+'_marker_genes_celltype_filtered.pdf')

# sc.pl.rank_genes_groups(adata,key='rank_genes_groups_filtered',n_genes=20,fontsize=8,save='_'+prefix+'_marker_genes_celltype_filtered.pdf')

# Bring in ortholog symbols (genes with no symbol but which have 1:1 human orthologs with symbols)
ortholog_symbols = pd.read_csv('biomart/'+'rna'+'-orthologs-with-symbols.txt',delimiter='\t',header=None,sep='\t',index_col=None)
# ortholog_symbols = pd.read_csv('biomart/'+prefix+'-orthologs-with-symbols.txt',delimiter='\t',header=None,sep='\t',index_col=None)

for i in range(0,len(ortholog_symbols[0])):
	gene = ortholog_symbols[0][i]
	symbol = ortholog_symbols[1][i] + '*'
	print('Renaming '+gene+' to '+symbol)
	adata.var['gene_short_name'].cat.rename_categories({gene:symbol},inplace=True)

def plot_umap(gene,cell_type=None):
	if not bool(np.isin(gene,adata.var.gene_short_name.to_list())):
		print('Cannot plot '+gene)
		return
	print(gene)
	a = sc.pl.umap(adata, gene_symbols='gene_short_name', color=gene,return_fig=True)
	a.axes[0].set_aspect(aspect=1)
	if cell_type is None:
		cell_suffix=''
	else:
		cell_suffix='_'+cell_type
	a.savefig('figures/marker_genes/umap_'+prefix+'_'+gene+cell_suffix+'.pdf')
	a.clear()
	plt.close('all')

# scanpy should be 1.9.1 for colorbar=None to work
def plot_umaps(genes,cell_type=None):
    good_genes = []
    for gene in genes:
            if bool(np.isin(gene,adata.var.gene_short_name.to_list())): good_genes.append(gene)
    ncolumns = math.ceil(math.sqrt(len(good_genes)))
    if ncolumns < 3:
    	ncolumns = 3
    a = sc.pl.umap(adata, gene_symbols='gene_short_name', color=good_genes, colorbar_loc=None,frameon=None,ncols=ncolumns,wspace=0,return_fig=True)
    for i in range(0,len(good_genes)):
            a.axes[i].set_aspect(aspect=1)
            a.axes[i].set_axis_off()
            a.axes[i].title.set_fontsize('24')
    if cell_type is None:
            raise ValueError('Parameter cell_type cannot be None')
    a.savefig('figures/marker_genes/umap_'+prefix+'_'+cell_type+'_all.pdf')
    a.clear()
    plt.close('all')


# Cell type	Markers
# Neuroepithelial cells	Nestin, SOX2, Notch1, HES1, HES3, E-cadherin, occludin.
# Radial glia	Vimentin, nestin, PAX6, HES1, HES5, GFAP, GLAST, BLBP, TN-C, N-cadherin, SOX2.
# Intermediate progenitors	TBR2, MASH1/Ascl1.
# Immature neurons	Doublecortin, beta III tubulin, NeuroD1, TBR1, stathmin 1.
# Oligodendrocyte precursor cells	PDGF receptor alpha, NG2.
# Mature oligodendrocytes	Olig 1, olig 2, olig 3, MBP, OSP, MOG, SOX10.
# Schwann cells	MPZ, NCAM, GAP43, S100, P75NTR.
# Astrocytes	GFAP, EAAT1/GLAST, EAAT2/GLT-1, glutamine synthetase, S100 beta, ALDH1L1.
# Microglia	TMEM119, CD11b, CD45, Iba1, CX3CR1, F4/80, CD68, CD40.
# Mature neurons	NeuN, MAP2, 160 kDa neurofilament medium, 200kDa neurofilament heavy, synaptophysin, PSD95.
# Glutamatergic neurons	VGLUT1, VGLUT2, NMDAR1, NMDAR2B, glutaminase, glutamine synthetase.
# GABAergic neurons	GABA transporter 1, GABAB receptors 1 and 2, GAD65, GAD67.
# Dopaminergic neurons	Tyrosine hydroxylase, dopamine transporter, FOXA2, GIRK2, Nurr1, LMX1B
# Serotonergic neurons	Tryptophan hydroxylase, serotonin transporter, Pet1.
# Cholinergic neurons	Choline acetyltransferase, vesicular acetylcholine transporter, acetylcholinesterase

# Smooth muscle cells can be identified by several markers of their differentiation/maturation, including α-SMA, calponin, SM22, and SM-MHC isoforms (SM1 and SM2). Expression of SM-MHC isoforms is strictly limited to fully differentiated SMCs.18 Ultrastructurally, SMCs have a complete external lamina, many pinocytotic vesicles, perinuclear mitochondria, and varying amounts of myofilaments and focal densities depending on their phenotypic and functional (contractile or synthetic) state. Smooth muscle cells exhibit extensive phenotypic modulation and plasticity.19,20
# ACTA2, CNN1, TAGLN, MYH11

favorite_markers = {
'microglia_fav': ['P2RY12','C3','PTPRC','ENTPD1','PLXDC2','SRGAP2'], # 'ATP8B4','PALD1',
'astrocyte_fav': ['SLC1A2','SLC1A3','RGS20','ADGRV1','GFAP','ALDH1A1'],
'oligodendrocyte_fav': ['MBP','PLP1','MOG','CLDN11','TF','SLC44A1'],
'oligodendrocyteprecursor_fav': ['NTN1','SEMA5A','PTPRZ1','UST','PCDH15','STK32A'],
'excitatoryneuron_fav': ['RBFOX3','FSTL4','CAMK2A','NELL2','CLSTN2','EPHA4'],
'inhibitoryneuron_fav': ['RBFOX3','GAD2','BTBD11','ZNF385D','GRIK1','ADARB2'],
'dopaminergicneuron_fav': ['DBH','SLC6A2','TH','DLK1','ROBO1','GNAS*'],
'vascular_fav': ['SLC7A1','ATP10A','ADGRL4','SLC6A20','PTGDS','CEMIP'],
'endothelial_fav': ['SLC7A1','RNASE1','ATP10A','ADGRL4','RERGL','ABCG2'], # ['RNASE1','EPAS1','SLC7A1','FLT1'],
'pericyte_fav': ['SLC6A20','SLC13A4','C7','OGN','PTGDS','CEMIP'],
'mesenchymalstemcell_fav': ['COL1A1','COL1A2','GAS1','S100A4','FN1','VIM'],
'radialglialcell_fav': ['ASPM','CENPF','MKI67','CENPE','ZNF43','IRS4'],
'serotinergicneuron_fav': ['TPH2','TRH','SLC6A4','DDC','TAC1','MAP1B'],
'mediumspinyneuron_fav': ['DRD2','CA12','PTPN7','CA12','PHACTR1','DGKB'], # PPP1R1B
'thalamicinterneuron_fav': ['RBFOX3','TAFA4','COL4A4','IRAG1','NTNG1','CACNA1C'],
'cerebellargranulecell_fav': ['RBFOX3','GABRA6','RAB37','GRM4','KCND2','FSTL5'],
'purkinjecell_fav': ['PCP2','ATP2A3','ITPR1','GRID2','ALDOC','SLC1A6'], #,'ATP2A3'],
'basketcell_fav': ['SORCS3','KIT','DUSP10','GRID2','FRMD4A','UNC5C'], #,'KIT'],
'ependymal_fav': ['GFAP','SPAG17','ADGB','LHB','ZBBX','DNAH9'],
'AHSGneuron_fav': ['RBFOX3','AHSG','APOB','HBE1','APOA2','APOH'], # unknown neurons (AHSG expression, RBFOX3 positive)
'F5neuron_fav': ['RBFOX3','F5','HS6ST3','CSMD1','KCNIP4','SCHIP1'], # unknown neurons (F5 expression, RBFOX3 positive)
'KIR3DL12neuron_fav': ['KIR3DL12','RBFOX3','KCNIP4','CTNNA2','MACF1','WFDC13'], # unknown neurons (KIR3DL1/KIR3DL2 expression, RBFOX3 positive)
'KIR3DL12microglia_fav': ['KIR3DL12','ARHGAP15','ENTPD1','PLXDC2','APBB1IP','ST6GALNAC3'], # unknown glia (KIR3DL1/KIR3DL2 expression, RBFOX3 negative)
}

for cell in favorite_markers:
    print(cell)
    plot_umaps(favorite_markers[cell],cell)


# https://www.abcam.com/neuroscience/neural-markers-guide
abcam_markers = {
'neuroepithelial_abcam': ['NES','SOX2','NOTCH1','HES1','HES3','CDH1','OCLN'], # neuroepithelial cells
'radialglia_abcam': ['VIM','NES','PAX6','HES1','HES5','GFAP','SLC1A3','FABP7','TNC','CDH2','SOX2'], # radial glial
'intermediateprogenitor_abcam': ['EOMES','ASCL1'], # intermediate progenitors
'immatureneuron_abcam': ['DCX','TUBB3','NEUROD1','TBR1','STMN1'], # immature neurons
'oligodendrocyteprecursor_abcam': ['PDGFRA','CSPG4'], # oligodendrocyte precursor cells
'oligodendrocyte_abcam': ['OLIG1','OLIG2','OLIG3','MBP','CLDN11','MOG','SOX10'], # oligodendrocyte
'schwanncell_abcam': ['MPZ','NCAM1','GAP43','S100A1','NGFR'], # Schwann cells
'astrocyte_abcam': ['GFAP','SLC1A3','SLC1A2','GLUL','S100B','ALDH1L1'], # astrocyte
'microglia_abcam': ['TMEM119','ITGAM','PTPRC','AIF1','CX3CR1','ADGRE1','CD68','CD40'], # microglia
'neuron_abcam': ['RBFOX3','MAP2','NEFM','NEFH','SYP','DLG4'], # mature neurons
'glutaminergicneuron_abcam': ['SLC17A7','SLC17A6','GRIN1','GRIN2B','GLS','GLUL'], # glutamatergic neurons
'GABAergicneuron_abcam': ['SLC6A1','GABBR1','GABBR2','GAD2','GAD1'], # GABAergic neurons
'dopaminergicneuron_abcam': ['TH','SLC6A3','FOXA2','KCNJ6','NR4A2','LMX1B'], # dopaminergic neurons
'serotonergicneuron_abcam': ['TPH1','SLC6A4','FEV'], # serotonergic neurons
'cholinergicneuron_abcam': ['CHAT','SLC18A3','ACHE'] # cholinergic neurons
}

# https://www.proteinatlas.org/humanproteome/single+cell+type
# Parse using parse-markers.R

proteinatlas_markers = {
'microglia_proteinatlas': ['FGF20','R3HDML','ATP8B4','PALD1','TBX18','RASGEF1C','ST8SIA6','KCNK13','P2RY12','NLRP10','SFMBT2','ATP10A','SRGAP2'],
'oligodendrocyte_proteinatlas': ['AGPAT4','VRK2','FUT8','PLEKHH1','SLC44A1','TESK2','PDE8A','FRYL','PLD1','PSEN1','SLCO1A2','NKAIN1','FOLH1','PRUNE2','ABCA2','MAN2A1','FRMD4B','STRN','PLCL1','DNAJC6','CDK18','MOB3B','PLP1','ATP8A1','ATG4C','TRPV5','HAPLN2','SORT1','PRR5L','MAP7','DOCK10','ERMN','ZNF189','ENPP2','TJAP1','TTLL7','SLAIN1','DOCK5','ST18','FAM13C','FAM124A','CNDP1','PIP4K2A','C9H10orf90','PDE1C','PIEZO2','SLC24A2','SLC5A11','SLC22A15','TMEM144','KIF6','JAM3','USP54','PXK','MOBP','SH3TC2','CARNS1','CNP','LRRC63','UGT8','PRR18','SLCO3A1','NT5DC1','TMTC2','TMEM151A','ZDHHC20','AATK','UNC5C','CTNNA3','CNTN2','CHRM5','GLDN','RNF220','COL4A5','DPYD','TMEM63A','MVB12B','OPALIN','MBP','TMEM235','MOG','TRIM59'],
'oligodendrocyteprecursor_proteinatlas': ['TLL1','NTN1','CACNG5','CACNG4','XYLT1','PTPRZ1','UST','COL9A1','SEMA5A','GPR17','LRRC4C','PCDH15','KLHL1','CA10','MMP16','USP24','STK32A','DSCAM','AMZ1','LRRN1','LUZP2','LHFPL3'],
'astrocyte_proteinatlas': ['PHKA1','COL5A3','MT3','MLC1','ACSBG1','EYA1','GLI3','GLIS3','SLC1A2','NCAN','SLC7A10','GFAP','NHSL1','BMPR1B','SLCO1C1','STON2','PRDM16','LRIG1','PSD2','RGS20','PAMR1','CACHD1','GRIN2C','SLC9C2','GABRG1','ETNPPL','RANBP3L','ADGRV1','NTSR2','HPSE2','ACBD7','GPC5','RNF182','ZNRF3','RYR3','TRIL'],
'inhibitoryneuron_proteinatlas': ['TAC1','TRPC7','HRH3','GLRA2','TACR1','GRPR','TRPC4','GAD2','HUNK','NXPH2','DLX1','GRIP2','SLC10A4','VAX1','BTBD11','ZNF385D','SLC35F4','MEPE','GRIP1','SYNPR','C8H8orf34','ENTPD3','GPR149','NHS','GRIN3A','SP9','TENM3','PLSCR5','PTCHD4'],
'excitatoryneuron_proteinatlas': ['SLC6A7','FSTL4','CAMK2A','HTR2A','EPHA4','ANO3','CBLN2','CDH22','MCHR2','GPR26','CLSTN2','NEUROD6','OR1Q1','COL24A1','GAP43','NWD2','ADAMTSL1','NELL2','RALYL','OR14I1','CPNE4']
}

deg_markers = {
'excitatoryneuron_deg': ['HS3ST4','POU6F2','MLIP','CSMD1','KCNIP4','KALRN','PHACTR1','DLGAP2','NRG3','SCHIP1'], # excitatory neurons
'mediumspinyneuron_deg': ['DRD2','CA12','CCDC88C','SMOC2','PHACTR1','DGKB','GPC5','PDE10A','CACNA2D3'], # medium spiny neurons
'inhibitoryneuron_deg': ['GRIK1','NXPH1','BTBD11','ADARB2','SOX6','CNTNAP2','ERBB4','GRIK2'], # inhibitory neurons
'dopaminergicneuron_deg': ['DBH','DLK1','TH','ROBO1','GNAS*','UCHL1'], # dopaminergic neurons
'serotinergicneuron_deg': ['TPH2','TRH','SLC6A4','TAC1','MAP1B','GNAS*'], # dopaminergic neurons
'thalamicneuron_deg': ['TAFA4','COL4A4','IRAG1','NTNG1','CACNA1C','RBFOX1'], # thalamic neurons
'granulecell_deg': ['RAB37','GRM4','ENSMMUG00000061812','KCND2','RBFOX1','FSTL5'], # cerebellar granule cells
'purkinjecell_deg': ['PCP2','ATP2A3','ITPR1','GRID2','ALDOC','PRKG1','SLC1A6'], # Purkinje cells
'basketneuron_deg': ['KIT','DUSP10','GRID2','FRMD4A','UNC5C','DPP6'], # basket neurons
'astrocyte_deg': ['ITGB4','ETNPPL','ALDH1A1','SLC1A2','GPC5','GPM6A','ADGRV1'], # astrocytes
'oligodendrocyte_deg': ['MBP','SLC5A11','CNDP1','TF','HAPLN2','QKI','RNF220','CTNNA3'], # oligodendrocytes
'oligodendrocyteprecursor_deg': ['STK32A','SLC22A3','C1QL1','VCAN','SUSD5','PCDH15','LHFPL3','TNR','XYLT1','DSCAM','LRP1B','PTPRZ1'], # oligodendrocyte precursor cells
'endothelial_deg': ['ADGRL4','RERGL','ABCG2','EMCN','PTPRG','DLC1','ATP10A','RBMS3'], # endothelial cells
'pericyte_deg': ['SLC6A20','SLC13A4','C7','OGN','PTGDS','ATP1A2','CEMIP','LAMA2'], # pericytes
'microglia_deg': ['STAB1','MRC1','CMAH','ENTPD1','PLXDC2','ST6GALNAC3','SRGAP2','APBB1IP','INPP5D','ENTPD1','ASTN2','SLC9A9'], # microglia
'radialglia_deg': ['ASPM','CENPF','CPA6','MKI67','CENPE','ZNF43','SLIT2','C4BPA','IRS4','TPR'], # radial glial cells
'mesenchymalstemcell_deg': ['GAS1','COL1A1','LOX','EMP1*','S100A4','COL1A2','FN1','VIM','PABPC1','DDX5'], # mesenchymal stem cells
'ependymal_deg': ['SPAG17','ADGB','ZBBX','CFAP100','LHB','DNAH9','CFAP299','HYDIN','TMEM232','CFAP54','CFAP43'], # ependymal cells
'AHSGneuron_deg': ['AHSG','APOB','ALB','HBE1','APOA2','APOH','KALRN','BEND5','KCNMA1','SCHIP1','HS6ST3','DLGAP2','KHDRBS2'], # unknown neurons (AHSG expression, RBFOX3 positive)
'F5neuron_deg': ['F5','HS6ST3','CSMD1','KCNIP4','SCHIP1','DLGAP2'], # unknown neurons (F5 expression, RBFOX3 positive)
'KIR3DL12neuron_deg': ['KIR3DL12','ENSMMUG00000028447','WFDC13','MDN1','KCNIP4','RBFOX3','CTNNA2','MACF1'], # unknown neurons (KIR3DL1/KIR3DL2 expression, RBFOX3 positive)
'KIR3DL12microglia_deg': ['KIR3DL12','ARHGAP15','ENTPD1','DOCK8','PLXDC2','APBB1IP','ASTN2','ST6GALNAC3'], # unknown glia (KIR3DL1/KIR3DL2 expression, RBFOX3 negative)
}

# deg_markers = {
# 'excitatoryneuron_deg': ['TBC1D10C','SPN','RAC2','IDO1','IL2RB','FCGR3','FCN1','CD3D','SELP','PRF1','GZMB'], # excitatory neurons
# 'inhibitoryneuron_deg': ['TAC3','SST','NPR3','ADAMTSL5'], # inhibitory neurons
# 'mediumspinyneuron_deg': ['DRD2','PTPN7','CA12','ADGRG7'], # medium spiny neurons
# 'dopaminergicneuron_deg': ['DBH','SLC6A2','TH','DLK1'], # dopaminergic neurons
# 'thalamicneuron_deg': ['RGS16','SHOX2','TAFA4','COL4A4'], # thalamic neurons
# 'purkinjecell_deg': ['PCP2','ATP2A3','ZNF385C','HOMER3'], # Purkinje cells
# 'basketneuron_deg': ['LIPG','TRIM67','TRAC','KIT','DUSP10'], # basket neurons
# 'granulecell_deg': ['GLI1','LCAT','PAX3','SDS'], # cerebellar granule cells
# 'neuralprogenitor_deg': ['ASPM','CENPF','CPA6','MKI67','C4BPA'], # neural progenitors
# 'radialglia_deg': ['COL1A1','GAS1','S100A4','ENSMMUG00000062238','COL1A2'], # radial glial cells
# 'astrocyte_deg': ['SERPINA3','IL33','ITPRID1','Metazoa_SRP'], # astrocytes
# 'oligodendrocyte_deg': ['OPALIN','SLC5A11','SLC45A3','CNDP1','CD22'], # oligodendrocytes
# 'oligodendrocyteprecursor_deg': ['NKAIN4','STK32A','SLC22A3','C1QL1','VCAN'], # oligodendrocyte precursor cells
# 'ependymal_deg': ['SPAG17','ADGB','LHB','ZBBX'], # ependymal cells
# 'endothelial_deg': ['MYOC','ITIH3','SFRP2','SLC6A20','FMO3'], # endothelial cells
# 'microglia_deg': ['STAB1','HPGDS','LYVE1','MRC1','MARCO','CD84','DAB2'] # microglia
# }

mckenzie_markers = {
'astrocyte_mckenzie': ['AQP4','BMPR1B','ETNPPL','GJB6','GJA1','FGFR3','SLC25A18','SLC1A2','SDC4','GFAP','EDNRB','ALDH1L1','CHI3L1','CLDN10','AGT','SLCO1C1'],
'endothelial_mckenzie': ['APOLD1','CD34','TGM2','IFI27','TM4SF1','ITIH5','SELE','TM4SF18','MECOM','VWF','ANXA3','ITGA1','CYSLTR2','ATP10A','EDN3','ERG'],
'microglia_mckenzie': ['CCL3','ITGAX','CD74','C1QB','FOLR2','TLR1','SLA','DHRS9','P2RX4','ARHGAP25','KBTBD8','TNFSF18','IL1A','HAVCR2','MSR1','FPR1'],
'neuron_mckenzie': ['SYNPR','RELN','CNR1','GAD2','OPRK1','GABRB2','RAB3C','SYT1','KCNC2','ZMAT4','RIMBP2','CHGB','GABRA1','MYT1L','GAD1','PTHLH'],
'oligodendrocyte_mckenzie': ['UGT8','PLP1','ERMN','CNDP1','CLDN11','CDH19','TF','FOLH1','KLK6','CNTN2','MOBP','SH3TC2','ST18','ERBB3','MYRF','MOG']
}

darmanis_markers = {
'oligodendrocyteprecusor1_darmanis': ['B3GNT7','BLM','CRISPLD2','GAB3','GALR1','HAS2','IL1RAP','LAMB4','LHFPL3','LRRK2','LUZP2','MEGF11','PCDH15','PDGFRA','RAB31','SEMA5A','SH2D4A','TNR','TREM1','XYLT1'],
'oligodendrocyteprecusor2_darmanis': ['C12H2orf88','DNER','FERMT1','GPNMB','GPR17','HR','HS6ST2','KY','LIMS2','NRIP3','PLD5','PLP1','SCG2','SENP8','SLC2A13','ZNF268'],
'oligodendrocyte_darmanis': ['ABCA8','CAPN3','CARNS1','CLDN11','CNDP1','DBNDD2','ENPP2','ERMN','FOLH1','GJB1','HHIP','KLK6','MOBP','OPALIN','RNASE1','TF','TMEM144','UGT8'],
'astrocyte1_darmanis': ['ACSBG1','AGT','AQP4','ATP13A4','BMPR1B','FGFR3','GFAP','GJA1','GJB6','GPR37L1','MGST1','RANBP3L','SDC4','SFXN5','SLC25A18','SLC39A12','SLC4A4','SLCO1C1'],
'microglia_darmanis': ['ALOX5AP','C3','C3AR1','CCL3','CD53','CD74','DHRS9','GPR183','IL1B','ITGAX','LAPTM5','LCP1','LILRB4','MSR1','OLR1','PLEK','PTPRC','RGS10'],
'astrocyte2_darmanis': ['AQP4','CALB2','CALN1','CLDN18','COL21A1','EDNRB','GJA1','HTR2C','IGFBP5','KCNT2','KDM5D','MTUS2','NLGN4Y','PROX1','PTHLH','SLC10A4','TGFB2','USP9Y','VIP'],
'neuron_darmanis': ['CAMK2A','CCK','CELF4','CHGB','CIT','DNM1','GABRA1','GABRB2','GABRG2','MAP7D2','NRXN3','SCG2','SCN2A','SNAP25','SYNPR','THY1','TMEM130','UNC80','VSNL1'],
'endothelial_darmanis': ['ANXA1','APOLD1','CEBPD','EDN3','ESAM','FOXC1','GNG11','ICAM1','IFI27','IL6','ITGA1','ITIH5','MSX1','PEAR1','TBX3','TGM2','TM4SF1']
}

zhu_markers = {
'Astro1_zhu': ['SLC1A2','ADGRV1','SLC1A3','NKAIN3','GPC5','GPM6A','APOE','NRXN1','ALDH1A1','MSI2'],
'Astro2_zhu': ['SLC1A2','SLC1A3','ADGRV1','NKAIN3','ENSMMUG00000000428','GPC5','NCKAP5','CDH20','MAG','PLP1'],
'Endo_zhu': ['RNASE1','EPAS1','EBF1','FLT1','SLC7A1','ITM2A','SLC2A1','GSN','ENSMMUG00000015144','MECOM'],
'ExN1_zhu': ['NAA15','SATB2','LDB2','NRG1','ARPP21','RTN1','HECW1','CTIF','PRKCB','KCNQ5'],
'ExN10_zhu': ['ENSMMUG00000044284','ENSMMUG00000040417','RPS11','KCNIP4','ARPP21','KCNQ5','NRG1','STXBP5L','ENSMMUG00000003851','OPCML'],
'ExN2_zhu': ['FAM19A1','DPYD','WFS1','PTPRK','ZBTB18','EPHA6','SERPINE2','GRIA4','CDH12','CACNA2D1'],
'ExN3_zhu': ['RORB','IL1RAPL2','TSHZ2','POU6F2','FSTL5','PDZRN4','CLSTN2','DCC','CPNE4','HECW1'],
'ExN4_zhu': ['KCNIP4','NRG1','PHACTR1','STXBP5L','KCNQ5','NKAIN2','RGS6','LRRTM4','LDB2','OPCML'],
'ExN5_zhu': ['SYT1','SNAP25','ENSMMUG00000028701','NRGN','CALM2','CALM1','CAMK2A','NEFL','CHN1','ENSMMUG00000020894'],
'ExN6_zhu': ['ENSMMUG00000009066','RPL26','ENSMMUG00000012537','KCNIP4','NRG1','NKAIN2','PHACTR1','KCNQ5','ARPP21','LDB2'],
'ExN7_zhu': ['TLE4','RXFP1','ASIC2','HTR2C','SEMA3E','EPHA5','SEMA3A','KIAA1217','ENSMMUG00000039690','FOXP2'],
'ExN8_zhu': ['ENSMMUG00000005253','C14ORF93','ENSMMUG00000009272','KCNIP4','KCNQ5','NKAIN2','ARPP21','CAMK2A','SATB2','PHACTR1'],
'ExN9_zhu': ['SLC1A2','SLC1A3','ADGRV1','NKAIN3','GPC5','APOE','ALDH1A1','SLC4A4','PDZRN4','MSI2'],
'InN1_zhu': ['ERBB4','NXPH1','ZNF385D','MYO16','ZNF804A','BTBD11','KCNC2','GAD1','DPP10','ANK1'],
'InN2_zhu': ['SST','GRIK1','NXPH1','SPOCK3','GRIK2','KIAA1217','SLC24A3','GRIN3A','NRXN3','KIF26B'],
'InN3_zhu': ['ADARB2','ERBB4','VIP','ZNF804A','TAC3','THSD7A','GALNTL6','CRH','CNTNAP2','SLC24A3'],
'InN4_zhu': ['ADARB2','CXCL14','RELN','ERBB4','GRIK2','GRIK1','CNTN5','GALNTL6','FSTL5','GAD2'],
'InN5_zhu': ['FBXL7','ENSMMUG00000012484','PTCHD4','SLC6A1','GRIK1','EYA4','NXPH1','CDH13','SGCZ','ADARB2'],
'Oligo_zhu': ['PLP1','MBP','MAG','RNF220','MOG','FA2H','TF','OPALIN','ANLN','ST18'],
'OPC_zhu': ['LHFPL3','PTPRZ1','APOD','PCDH15','COL9A1','SEMA5A','LUZP2','XYLT1','SOX6','VCAN'],
'Peri_zhu': ['PTGDS','SLC13A4','MFAP4','IGFBP5','BNC2','DSP','SLC6A20','VIM','CPAMD8','SLC22A8']
}

kang_celltypes = {
'astrocytes_kangcelltypes': ['ALDOC','GFAP','S100B'],
'corticalgabainterneurons_kangcelltypes': ['CALB2','CALB1','NOS1','PVALB','CCK','VIP','DLX1','DLX2','NKX2-1','ASCL1','GAD1','GAD2'],
'corticalglutamatergicneurons_kangcelltypes': ['RELN','CUX1','UNC5D','RORB','BCL11B','ETV1','FEZF2','OTX1','FOXP2','NTSR1','SOX5','SSTR2','TBR1','TLE4','ZFPM2','CTGF','UNC5C'],
'microglia_kangcelltypes': ['CFH','FCER1G','TNIP2'],
'oligodendrocytes_kangcelltypes': ['CNP','CSPG4','OLIG1','OLIG2','PDGFRA']
}

kang_neurodevelopmental = {
'cellproliferationmarker_kangneurodevelopmental': ['HES1','MKI67','NES','RRM1'],
'gabashift_kangneurodevelopmental': ['SLC12A1','SLC12A2','CA2','CA7','SLC12A4','SLC12A5','SLC12A6','SLC12A7'],
'myelination_kangneurodevelopmental': ['C11orf9','MAG','MBP','MOG','PLP1'],
'neurondevelopment_kangneurodevelopmental': ['CNTN2','CAMK2A','MAP1A','MAPT','MAP1B','MAP2','TUBB','DCX','SYN1','SYP','SYPL1','SYPL2']
}

kang_neurotransmitters = {
'adenosinergicsystem_kangneurotransmitters': ['ADORA1','ADORA2A','ADORA2B','ADORA3'],
'cholinergicsystem_kangneurotransmitters': ['CHAT','ACHE','CHRNA1','CHRNA2','CHRNA3','CHRNA4','CHRNA5','CHRNA6','CHRNA7','CHRNA9','CHRNA10','CHRNB1','CHRNB2','CHRNB3','CHRNB4','CHRND','CHRNE','CHRNG','CHRM1','CHRM2','CHRM3','CHRM4','CHRM5'],
'dopaminergicsystem_kangneurotransmitters': ['TH','DDC','COMT','DRD1','DRD5','DRD2','DRD3','DRD4'],
'gabaergicsystem_kangneurotransmitters': ['GAD1','GAD2','GABRA1','GABRA2','GABRA3','GABRA4','GABRA5','GABRA6','GABRB1','GABRB2','GABRB3','GABRG1','GABRG2','GABRG3','GABRD','GABRE','GABRP','GABRQ','GABRR1','GABRR2','GABBR1','GABBR2'],
'glutamergicsystem_kangneurotransmitters': ['GLS','SLC1A3','SLC1A2','SLC1A1','SLC1A6','SLC1A7','SLC17A7','SLC17A6','GRIA1','GRIA2','GRIA3','GRIA4','GRIK1','GRIK2','GRIK3','GRIK4','GRIK5','GRIN1','GRIN2A','GRIN2B','GRIN2C','GRIN2D','GRIN3A','GRM1','GRM5','GRM2','GRM3','GRM4','GRM6','GRM7','GRM8'],
'glycinergicsystem_kangneurotransmitters': ['SHMT1','SHMT2','GLRA1','GLRA2','GLRA3','GLRB'],
'histaminergicsystem_kangneurotransmitters': ['HDC','HNMT','ABP1','HRH1','HRH2','HRH3','HRH4'],
'noradrenalinergicandadrenergicsystem_kangneurotransmitters': ['DBH','PNMT','SLC6A2','ADRA1A','ADRA1B','ADRA1D','ADRA2A','ADRA2B','ADRA2C','ADRB1','ADRB2','ADRB3'],
'serotoninergicsystem_kangneurotransmitters': ['TPH1','TPH2','DDC','MAOA','MAOB','HTR1A','HTR1B','HTR1D','HTR1E','HTR1F','HTR2A','HTR2B','HTR2C','HTR3A','HTR3B','HTR3C','HTR3D','HTR3E','HTR4','HTR5A','HTR6','HTR7']
}

# litsearch

lit_markers = {
'mediumspinyneuron_lit': ['PPP1R1B'], # https://brainxell.com/medium-spiny-neurons
'microglia_lit': ['TREM2'], # https://brainxell.com/microglia
'glutamatergicneuron_lit': ['BCL11B'], # https://brainxell.com/cortical-glutamatergic-neurons (layer V [internal pyramidal layer] - specific) https://brainxell.com/layer-v-glutamatergic-neurons
'granulecell_lit': ['RBFOX3','GABRA6','CBLN3'], # https://elifesciences.org/articles/37551; GABRA6 and CBLN3 from Hannah's mousebrain file, #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7918299/ (CNTN3 and TUBB3 unhelpful)
'basketneuron_lit': ['SORCS3'], # https://elifesciences.org/articles/37551
'purkinjecell_lit': ['ITPR1','SLC1A6','CA8'] # https://elifesciences.org/articles/37551, ADCY1,ITPKA,CSDC2, minimally helpful  (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7918299/)
}


for cell in abcam_markers:
    print(cell)
    plot_umaps(abcam_markers[cell],cell)

for cell in proteinatlas_markers:
    print(cell)
    plot_umaps(proteinatlas_markers[cell],cell)

for cell in deg_markers:
    print(cell)
    plot_umaps(deg_markers[cell],cell)

for cell in mckenzie_markers:
    print(cell)
    plot_umaps(mckenzie_markers[cell],cell)

for cell in zhu_markers:
    print(cell)
    plot_umaps(zhu_markers[cell],cell)

for cell in kang_neurotransmitters:
    print(cell)
    plot_umaps(kang_neurotransmitters[cell],cell)

for cell in kang_celltypes:
    print(cell)
    plot_umaps(kang_celltypes[cell],cell)

for cell in kang_neurodevelopmental:
    print(cell)
    plot_umaps(kang_neurodevelopmental[cell],cell)

for cell in darmanis_markers:
    print(cell)
    plot_umaps(darmanis_markers[cell],cell)

for cell in lit_markers:
    print(cell)
    plot_umaps(lit_markers[cell],cell)



# Thalamus is comprised of 3 basic cell types: relay cells, interneurons, and cells of the thalamic reticular nucleus (Fig. 1) (for details, see Sherman and Guillery, 1996; Sherman and Guillery, 2013). Each of these may be further subdivided, but the complete classification of these cell types has yet to be done. Relay cells are glutamatergic, but interneurons and reticular cells are GABAergic, providing a major inhibitory input to relay cells. Interneurons are located within the main dorsal thalamic relay nuclei, intermixed with relay cells, and the ratio of interneurons to relay cells is roughly 1–32.

# for cell in abcam_markers:
# 	print(cell)
# 	for gene in abcam_markers[cell]: plot_umap(gene,cell)
