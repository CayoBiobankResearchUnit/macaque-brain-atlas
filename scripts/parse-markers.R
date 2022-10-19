#!/usr/bin/env

cell.types = c('Microglial','Oligodendrocytes_Cell','Oligodendrocyte','Astrocytes_Cell','Inhibitory','Excitatory')
cell.names = c('microglia','oligodendrocyte','oligodendrocyteprecursor','astrocyte','inhibitoryneuron','excitatoryneuron')

cell.marker.files = file.path('data/proteinatlas',paste0('cell_type_category_rna_',cell.types,'.tsv'))

cell.markers = lapply(cell.marker.files,read.delim)

names(cell.markers) = cell.types

library(biomaRt)

hsap = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset=paste0('hsapiens','_gene_ensembl'))

cell.markers.out = lapply(cell.markers,function(x) {
	getBM(attributes=c('ensembl_gene_id','external_gene_name',
					paste0('mmulatta','_homolog_ensembl_gene'),
					paste0('mmulatta','_homolog_associated_gene_name'),
					paste0('mmulatta','_homolog_orthology_type')),
		filters = 'ensembl_gene_id',
		values = x$Ensembl,
		mart = hsap)
})

names(cell.markers.out) = cell.names

cell.markers.out = lapply(cell.markers.out,function(x) {
	subset(x,mmulatta_homolog_orthology_type=='ortholog_one2one')$mmulatta_homolog_associated_gene_name
})

cat(paste0('proteinatlas_markers = {\n',do.call(paste,c(
lapply(names(cell.markers.out),function(i) {
	paste0('\'',i,'_proteinatlas\': [',paste(paste0('\'',cell.markers.out[[i]],'\''),collapse=','),']')
}),
list(sep=',\n'))),'\n}'))

# Supplemental File 1, top_human_specificity sheet
# McKenzie et al. 2018, https://doi.org/10.1038/s41598-018-27293-5
mckenzie = read.delim('data/markergenes/mckenzie-human-specificity.txt')
mckenzie$Celltype = factor(mckenzie$Celltype)
mckenzie$Celltype = factor(mckenzie$Celltype,levels=c('ast','end','mic','neu','oli'),labels=c('astrocyte','endothelial','microglia','neuron','oligodendrocyte'))

mckenzie = do.call(rbind,lapply(split(mckenzie,mckenzie$Celltype),function(x) head(x,32)))

mckenzie.out = with(mckenzie,
	merge(getBM(attributes=c('ensembl_gene_id','external_gene_name',
					paste0('mmulatta','_homolog_ensembl_gene'),
					paste0('mmulatta','_homolog_associated_gene_name'),
					paste0('mmulatta','_homolog_orthology_type')),
		filters = 'external_gene_name',
		values = gene,
		mart = hsap),mckenzie,by.x='external_gene_name',by.y='gene',all.x=TRUE,all.y=FALSE)
)
mckenzie.out = subset(mckenzie.out,mmulatta_homolog_orthology_type == 'ortholog_one2one')
mckenzie.out = mckenzie.out[order(mckenzie.out$grand_mean,decreasing=TRUE),]
mckenzie.out = lapply(split(mckenzie.out,mckenzie.out$Celltype),function(x) head(x,16)$mmulatta_homolog_associated_gene_name)

cat(paste0('mckenzie_markers = {\n',do.call(paste,c(
lapply(names(mckenzie.out),function(i) {
	paste0('\'',i,'_mckenzie\': [',paste(paste0('\'',mckenzie.out[[i]],'\''),collapse=','),']')
}),
list(sep=',\n'))),'\n}'))

# Table S3
# Darmanis et al. 2015, https://doi.org/10.1073/pnas.1507125112
darmanis = read.delim('data/markergenes/darmanis-unbiased-clusters.txt')
darmanis = reshape2::melt(darmanis,id=0)
darmanis$cell = factor(darmanis$variable,
	levels=c('Cluster1','Cluster2','Cluster3','Cluster4','Cluster5','Cluster6','Cluster7','Cluster8','Cluster9','Cluster10'),
	labels=c('oligodendrocyteprecusor1','oligodendrocyteprecusor2','oligodendrocyte','astrocyte1','microglia','astrocyte2','neuron','endothelial','fetalcells1','fetalcells2')
)
darmanis$gene = darmanis$value

darmanis.out = with(darmanis,
	merge(getBM(attributes=c('ensembl_gene_id','external_gene_name',
					paste0('mmulatta','_homolog_ensembl_gene'),
					paste0('mmulatta','_homolog_associated_gene_name'),
					paste0('mmulatta','_homolog_orthology_type')),
		filters = 'external_gene_name',
		values = gene,
		mart = hsap),darmanis,by.x='external_gene_name',by.y='gene',all.x=TRUE,all.y=FALSE)
)

darmanis.out = droplevels(subset(darmanis.out,(!cell %in% c('fetalcells1','fetalcells2')) & mmulatta_homolog_orthology_type == 'ortholog_one2one',select=c('cell','mmulatta_homolog_associated_gene_name')))

darmanis.out = lapply(split(darmanis.out,darmanis.out$cell),function(x) x$mmulatta_homolog_associated_gene_name)

cat(paste0('darmanis_markers = {\n',do.call(paste,c(
lapply(names(darmanis.out),function(i) {
	paste0('\'',i,'_darmanis\': [',paste(paste0('\'',darmanis.out[[i]],'\''),collapse=','),']')
}),
list(sep=',\n'))),'\n}'))


zhu = data.frame(gene = c('SLC1A2', 'ADGRV1', 'SLC1A3', 'NKAIN3', 'GPC5', 'GPM6A', 'APOE', 
'NRXN1', 'ALDH1A1', 'MSI2', 'SLC1A2', 'SLC1A3', 'ADGRV1', 'NKAIN3', 
'ENSMMUG00000000428', 'GPC5', 'NCKAP5', 'CDH20', 'MAG', 'PLP1', 
'RNASE1', 'EPAS1', 'EBF1', 'FLT1', 'SLC7A1', 'ITM2A', 'SLC2A1', 
'GSN', 'ENSMMUG00000015144', 'MECOM', 'NAA15', 'SATB2', 'LDB2', 
'NRG1', 'ARPP21', 'RTN1', 'HECW1', 'CTIF', 'PRKCB', 'KCNQ5', 
'ENSMMUG00000044284', 'ENSMMUG00000040417', 'RPS11', 'KCNIP4', 
'ARPP21', 'KCNQ5', 'NRG1', 'STXBP5L', 'ENSMMUG00000003851', 'OPCML', 
'FAM19A1', 'DPYD', 'WFS1', 'PTPRK', 'ZBTB18', 'EPHA6', 'SERPINE2', 
'GRIA4', 'CDH12', 'CACNA2D1', 'RORB', 'IL1RAPL2', 'TSHZ2', 'POU6F2', 
'FSTL5', 'PDZRN4', 'CLSTN2', 'DCC', 'CPNE4', 'HECW1', 'KCNIP4', 
'NRG1', 'PHACTR1', 'STXBP5L', 'KCNQ5', 'NKAIN2', 'RGS6', 'LRRTM4', 
'LDB2', 'OPCML', 'SYT1', 'SNAP25', 'ENSMMUG00000028701', 'NRGN', 
'CALM2', 'CALM1', 'CAMK2A', 'NEFL', 'CHN1', 'ENSMMUG00000020894', 
'ENSMMUG00000009066', 'RPL26', 'ENSMMUG00000012537', 'KCNIP4', 
'NRG1', 'NKAIN2', 'PHACTR1', 'KCNQ5', 'ARPP21', 'LDB2', 'TLE4', 
'RXFP1', 'ASIC2', 'HTR2C', 'SEMA3E', 'EPHA5', 'SEMA3A', 'KIAA1217', 
'ENSMMUG00000039690', 'FOXP2', 'ENSMMUG00000005253', 'C14ORF93', 
'ENSMMUG00000009272', 'KCNIP4', 'KCNQ5', 'NKAIN2', 'ARPP21', 
'CAMK2A', 'SATB2', 'PHACTR1', 'SLC1A2', 'SLC1A3', 'ADGRV1', 'NKAIN3', 
'GPC5', 'APOE', 'ALDH1A1', 'SLC4A4', 'PDZRN4', 'MSI2', 'ERBB4', 
'NXPH1', 'ZNF385D', 'MYO16', 'ZNF804A', 'BTBD11', 'KCNC2', 'GAD1', 
'DPP10', 'ANK1', 'SST', 'GRIK1', 'NXPH1', 'SPOCK3', 'GRIK2', 
'KIAA1217', 'SLC24A3', 'GRIN3A', 'NRXN3', 'KIF26B', 'ADARB2', 
'ERBB4', 'VIP', 'ZNF804A', 'TAC3', 'THSD7A', 'GALNTL6', 'CRH', 
'CNTNAP2', 'SLC24A3', 'ADARB2', 'CXCL14', 'RELN', 'ERBB4', 'GRIK2', 
'GRIK1', 'CNTN5', 'GALNTL6', 'FSTL5', 'GAD2', 'FBXL7', 'ENSMMUG00000012484', 
'PTCHD4', 'SLC6A1', 'GRIK1', 'EYA4', 'NXPH1', 'CDH13', 'SGCZ', 
'ADARB2', 'PLP1', 'MBP', 'MAG', 'RNF220', 'MOG', 'FA2H', 'TF', 
'OPALIN', 'ANLN', 'ST18', 'LHFPL3', 'PTPRZ1', 'APOD', 'PCDH15', 
'COL9A1', 'SEMA5A', 'LUZP2', 'XYLT1', 'SOX6', 'VCAN', 'PTGDS', 
'SLC13A4', 'MFAP4', 'IGFBP5', 'BNC2', 'DSP', 'SLC6A20', 'VIM', 
'CPAMD8', 'SLC22A8'), type = c('Astro1', 'Astro1', 'Astro1', 'Astro1', 'Astro1', 'Astro1', 
'Astro1', 'Astro1', 'Astro1', 'Astro1', 'Astro2', 'Astro2', 'Astro2', 
'Astro2', 'Astro2', 'Astro2', 'Astro2', 'Astro2', 'Astro2', 'Astro2', 
'Endo', 'Endo', 'Endo', 'Endo', 'Endo', 'Endo', 'Endo', 'Endo', 
'Endo', 'Endo', 'ExN1', 'ExN1', 'ExN1', 'ExN1', 'ExN1', 'ExN1', 
'ExN1', 'ExN1', 'ExN1', 'ExN1', 'ExN10', 'ExN10', 'ExN10', 'ExN10', 
'ExN10', 'ExN10', 'ExN10', 'ExN10', 'ExN10', 'ExN10', 'ExN2', 
'ExN2', 'ExN2', 'ExN2', 'ExN2', 'ExN2', 'ExN2', 'ExN2', 'ExN2', 
'ExN2', 'ExN3', 'ExN3', 'ExN3', 'ExN3', 'ExN3', 'ExN3', 'ExN3', 
'ExN3', 'ExN3', 'ExN3', 'ExN4', 'ExN4', 'ExN4', 'ExN4', 'ExN4', 
'ExN4', 'ExN4', 'ExN4', 'ExN4', 'ExN4', 'ExN5', 'ExN5', 'ExN5', 
'ExN5', 'ExN5', 'ExN5', 'ExN5', 'ExN5', 'ExN5', 'ExN5', 'ExN6', 
'ExN6', 'ExN6', 'ExN6', 'ExN6', 'ExN6', 'ExN6', 'ExN6', 'ExN6', 
'ExN6', 'ExN7', 'ExN7', 'ExN7', 'ExN7', 'ExN7', 'ExN7', 'ExN7', 
'ExN7', 'ExN7', 'ExN7', 'ExN8', 'ExN8', 'ExN8', 'ExN8', 'ExN8', 
'ExN8', 'ExN8', 'ExN8', 'ExN8', 'ExN8', 'ExN9', 'ExN9', 'ExN9', 
'ExN9', 'ExN9', 'ExN9', 'ExN9', 'ExN9', 'ExN9', 'ExN9', 'InN1', 
'InN1', 'InN1', 'InN1', 'InN1', 'InN1', 'InN1', 'InN1', 'InN1', 
'InN1', 'InN2', 'InN2', 'InN2', 'InN2', 'InN2', 'InN2', 'InN2', 
'InN2', 'InN2', 'InN2', 'InN3', 'InN3', 'InN3', 'InN3', 'InN3', 
'InN3', 'InN3', 'InN3', 'InN3', 'InN3', 'InN4', 'InN4', 'InN4', 
'InN4', 'InN4', 'InN4', 'InN4', 'InN4', 'InN4', 'InN4', 'InN5', 
'InN5', 'InN5', 'InN5', 'InN5', 'InN5', 'InN5', 'InN5', 'InN5', 
'InN5', 'Oligo', 'Oligo', 'Oligo', 'Oligo', 'Oligo', 'Oligo', 
'Oligo', 'Oligo', 'Oligo', 'Oligo', 'OPC', 'OPC', 'OPC', 'OPC', 
'OPC', 'OPC', 'OPC', 'OPC', 'OPC', 'OPC', 'Peri', 'Peri', 'Peri', 
'Peri', 'Peri', 'Peri', 'Peri', 'Peri', 'Peri', 'Peri'),stringsAsFactors=FALSE)

# zhu.out = with(zhu,
# 	merge(getBM(attributes=c('ensembl_gene_id','external_gene_name',
# 					paste0('mmulatta','_homolog_ensembl_gene'),
# 					paste0('mmulatta','_homolog_associated_gene_name'),
# 					paste0('mmulatta','_homolog_orthology_type')),
# 		filters = 'external_gene_name',
# 		values = gene,
# 		mart = hsap),zhu,by.x='external_gene_name',by.y='gene',all.x=TRUE,all.y=FALSE)
# )

zhu.out = lapply(split(zhu,zhu$type),function(x) x$gene)

cat(paste0('zhu_markers = {\n',do.call(paste,c(
lapply(names(zhu.out),function(i) {
	paste0('\'',i,'_zhu\': [',paste(paste0('\'',zhu.out[[i]],'\''),collapse=','),']')
}),
list(sep=',\n'))),'\n}'))





# Kang et al. 2012, https://doi.org/10.1038/nature10523
# Table S12
class.names = c('neurodevelopmental','neurotransmitters','celltypes')

kang = do.call(rbind,lapply(class.names,function(i) {
	out = read.delim(file.path('data/markergenes',paste0('kang-',i,'.txt')))[c('Functional.group','Gene.symbol')]
	out$gene = out$Gene.symbol
	out$group = gsub(' ','',tolower(out$Functional.group))
	out$class = i
	out[c('class','group','gene')]
}))

kang.out = with(kang,
	merge(getBM(attributes=c('ensembl_gene_id','external_gene_name',
					paste0('mmulatta','_homolog_ensembl_gene'),
					paste0('mmulatta','_homolog_associated_gene_name'),
					paste0('mmulatta','_homolog_orthology_type')),
		filters = 'external_gene_name',
		values = gene,
		mart = hsap),kang,by.x='external_gene_name',by.y='gene',all.x=TRUE,all.y=FALSE)
)

kang.out = droplevels(subset(kang.out,mmulatta_homolog_orthology_type == 'ortholog_one2one'))

kang.out = lapply(split(kang,kang$class),function(x) {
	lapply(split(x,x$group),function(x) x$gene)
})


cat(do.call(paste,c(lapply(names(kang.out),function(i) {
	paste0(
		'kang_',i,' = {\n',
		do.call(paste,c(lapply(
			names(kang.out[[i]]),
			function(j) {
				paste0('\'',j,'_kang',i,'\': [',paste(paste0('\'',kang.out[[i]][[j]],'\''),collapse=','),']')
			}
		),list(sep=',\n'))),
		'\n}'
	)
}),list(sep='\n\n'))))



lake = structure(list(gene = c("TMEM90A", "SLC17A6", "GLRA2", "EBF3", 
"DLX5", "TBR1", "ISLR2", "LHX5", "TMEM130", "MAB21L1", "DLX6", 
"ST8SIA2", "KCNK9", "LHFPL5", "COL25A1", "ECEL1", "NKX2-1", "SRRM4", 
"NELL1", "MEG3", "SCUBE3", "NOS1", "GPR83", "KCNH7", "TACR1", 
"SLC32A1", "HCN4", "KLHL1", "VSTM2L", "RELN", "CACNA2D2", "ANKRD35", 
"DLX2", "CELF4", "NHLH2", "KCNJ5", "TRIM66", "ROBO2", "MAGEL2", 
"HS3ST5", "GREM2", "CELF6", "SST", "PCSK1", "NPPC", "DLX1", "DPYSL5", 
"NTRK1", "LHX8", "P2RX5", "SNHG11", "RIPK4", "NPAS4", "L1CAM", 
"NPY", "GDF5", "LHX1", "SP9", "BCL11A", "SLC10A4", "PNOC", "MARCH4", 
"DYNC1I1", "CHODL", "CRABP1", "BHLHE22", "IGFBPL1", "NXPH4", 
"GRM2", "LHX6", "SLC18A3", "SP8", "MRAP2", "CLSTN2", "FIBIN", 
"ECSCR", "RAB20", "CCL3", "CLDN11", "IRF5", "IL10RA", "PPP1R14A", 
"C1QC", "TGFB1", "RASAL3", "TNFRSF1B", "PLP1", "TMEM119", "TYROBP", 
"NCF4", "OPALIN", "TNF", "APOD", "UNC93B1", "CSF1R", "LAPTM5", 
"ERMN", "CD83", "LTC4S", "NCF1", "CX3CR1", "GNA15", "CTSS", "P2RY13", 
"NKX6-2", "CCRL2", "TLR7", "RGS1", "TLR2", "TREM2", "C5AR1", 
"P2RY6", "CAPN3", "GPR183", "GJB1", "C1QB", "CD33", "TMEM63A", 
"CD37", "INPP5D", "SLCO2B1", "SELPLG", "CNP", "CSF3R", "SLC25A45", 
"HMHA1", "VAV1", "FA2H", "LRCH4", "GDF15", "ENPP6", "CRYBB1", 
"GRAP", "CCL4", "ZC3H12A", "SLC11A1", "CCR1", "C1QA", "TMEM125", 
"EVI2A", "OLFML3", "TMEM88B", "ITGAM", "TEK", "MBP", "NFAM1", 
"MAG", "NOSTRIN", "IL1A", "NLRP3", "FCGR2B", "MOBP", "LRRC25", 
"GSN", "OSM", "CCR5", "TLR9", "SLC2A5", "MOG", "MAL", "GNG11", 
"GNGT2", "ASPA", "RCSD1", "PLD4", "SLC15A3"), type = structure(c(2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L), .Label = c("glia", "neuron"), class = "factor")), .Names = c("gene", 
"type"), row.names = c(NA, -166L), class = "data.frame")

spaethling = structure(list(gene = c("TMEM125", "NR0B1", "MOG", "HS3ST5", 
"ZNF488", "GAL3ST1", "CARNS1", "LCNL1", "TSPAN15", "PRKCQ", "LGI3", 
"GJB1", "MIR219A2", "ADAP1", "GJC2", "KCNS3", "RBP7", "SEC14L5", 
"CMTM5", "KLK6", "NINJ2", "PPP1R16B", "TMEM151A", "SPOCK3", "NKAIN2", 
"KLHL32", "SLC45A3", "RP11.334C17.5", "RP4.635A23.6", "RP4.555D20.2", 
"FA2H", "MAG", "CDKN1C", "C10orf90", "SOX10", "PTHLH", "GPD1", 
"RP11.358L22.3", "ST6GALNAC3", "RAPGEF3", "CNTNAP4", "DEPDC7", 
"CHADL", "MYRF", "RNASE1", "GALNT13", "UGT8", "FAM13A.AS1", "ENPP4", 
"AC018647.3", "FOLH1", "KCNH8", "RNF125", "ABCG1", "ADAMTS4", 
"ITM2A", "CNDP1", "LDB3", "MAL", "MTG1", "AATK", "ERMN", "SLC7A14", 
"MAP7", "PEX5L", "STXBP6", "CDH20", "CYB5R2", "SLC22A15", "APOD", 
"LINC.PINT", "PRR18", "CLDN11", "MOBP", "ATP10B", "CXorf57", 
"CDK18", "FAM124A", "SEMA4D", "SHROOM1", "PCSK6", "PLEKHG3", 
"PTPDC1", "MOB3B", "SOCS2", "TF", "ZNF365", "NKX6.2", "ENPP2", 
"MEGF10", "HDAC11", "SPNS2", "GPR62", "SLC31A2", "X4.Sep", "ARHGEF37", 
"TMEM63A", "BCAS1", "FRMD4B", "ST18", "HSPA2", "FAHD2B", "INPP5F", 
"PTGDS", "RAB33A", "FANCF", "RGS3", "CDKL5", "COBL", "ENOSF1", 
"FN1", "CSAG1", "TINAGL1", "HSPB2", "TBX3", "FAM26E", "FOXF2", 
"ACTA2.AS1", "CCDC3", "COL6A2", "NOX4", "RP11.383H13.1", "ISLR", 
"MARVELD1", "TMEM119", "PRKG1", "RP11.479G22.8", "MYL9", "VGLL3", 
"COX7A1", "LTBP2", "CCDC74A", "TAGLN", "TUSC3", "TGFBI", "COL7A1", 
"PAPSS2", "VDR", "ENPP1", "MRVI1", "PRELP", "COL12A1", "FHOD1", 
"SYNPO2", "LUM", "BGN", "PCDH18", "GPR176", "THBS1", "FHL2", 
"MAP3K7CL", "ITPRIP", "PLAC9", "CYR61", "NTN4", "FANCE", "IGFBP4", 
"LMO7", "MGP", "COL1A1", "COL3A1", "ATP8B1", "RGS5", "TPM2", 
"FOXC1", "NR2F2", "DCN", "PDLIM1", "SLC9A1", "ACTA2", "C11orf96", 
"PHACTR2", "COL1A2", "CYP1B1", "UACA", "CAV2", "HPGDS", "MPEG1", 
"LINC01272", "RGS1", "HCST", "RHBDF2", "MS4A7", "IL18", "CYTH4", 
"CD53", "PIK3AP1", "C15orf48", "PLEK", "TNFRSF1B", "NCKAP1L", 
"CD84", "CD83", "SLC1A5", "EGFR", "RP11.403A3.2", "HOXA4", "IL13RA2", 
"CHRDL1", "HOXD.AS2", "ATP13A4", "LAMP5", "PCOLCE2", "FZD5", 
"TRIL", "C21orf119", "SLC2A10", "ZFHX4.AS1", "ADORA1", "RP11.403A3.3", 
"SLC15A2", "C1QTNF5", "AFAP1L1", "HERC5", "KCNN2", "ROM1", "SERPINA3", 
"SOX21.AS1", "APLN", "LCTL", "GSG1L", "KATNAL2", "PPP1R1C", "SHC3", 
"CGNL1", "HOXA10", "C2orf88", "ADAMTS5", "HRH1", "NR2E1", "THSD1", 
"IL33", "ACKR3", "DIRAS3", "GRB14", "MAOB", "CPNE4", "DBNDD1", 
"SHQ1", "VAV3", "BBOX1", "COL22A1", "IFIT2", "NOV", "GLIDR", 
"RP11.421L21.3", "ST20", "HLF", "OLFM2", "SLITRK5", "GLI3", "MBOAT7", 
"EYA2", "C5orf56", "KCNN3", "GRIA1", "MARVELD3", "CTD.2044J15.2", 
"PSEN2", "EYA4", "NETO2", "CRISPLD1", "MASP1", "HMGN5", "TRIM5", 
"ASTN2", "FGFR3", "COX10", "LUZP2", "NCKAP5", "RP11.305K5.1", 
"PRKD1", "PPARGC1A", "TSPAN11", "C1S", "C2orf72", "LMAN2L", "SLC16A9", 
"EPHB1", "MIPEP", "TTC23", "TTLL4", "FGFBP3", "IGDCC4", "NKAIN3", 
"TP53I3", "ABHD17C", "AQP4", "FIBIN", "HEY1", "MORC4", "RILPL1", 
"UST", "FAM86DP", "GJA1", "LRRN1", "PXYLP1", "BEND6", "POSTN", 
"TRAF3IP2.AS1", "PIPOX", "KCNIP1", "GRIK3", "NEDD9", "FAM101B", 
"UG0898H09", "KIAA1549L", "RFX4", "SLC4A4", "CHST2", "JAM2", 
"GRIA3", "POMGNT2", "APCDD1", "FOXRED2", "ELN", "TMEM42", "PLCD3", 
"HES1", "TSPAN6", "TMEM8B", "RGMA", "MLC1", "PID1", "FAM198B", 
"SERPINE2", "CWC27", "SCRG1", "NCAN", "GFAP", "SOX11", "MET", 
"CELF3", "NKAIN1", "TTC9B", "VAT1L", "EYA1", "NEFL", "SHD", "PLK2", 
"NNAT", "GABRB3", "SLC16A14", "CTB.31O20.2", "CHGB", "TAGLN3", 
"JPH4", "RPH3A", "INSM1", "NXPH1", "KIF21B", "DLL3", "GAD1", 
"RNF165", "BEND5", "GDAP1L1", "ELAVL4", "SBK1", "DCX", "MPPED2", 
"SNAP25", "ATCAY", "CELF5", "MMP24", "NRXN1", "TSPAN13"), type = structure(c(5L, 
5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 
5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 
5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 
5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 
5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 
5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 
5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 
3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
3L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 
4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L), .Label = c("astrocyte", 
"endothelial", "microglial", "neuron", "oligodendrocyte"), class = "factor")), .Names = c("gene", 
"type"), row.names = c(NA, -366L), class = "data.frame")

marker_genes = structure(list(type = c("Oligodendrocyte", "Oligodendrocyte", 
"Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte", 
"Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte", 
"Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte", 
"Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte", 
"Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte", 
"Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte", "Oligodendrocyte", 
"Oligodendrocyte", "Oligodendrocyte", "Astrocyte", "Astrocyte", 
"Astrocyte", "Astrocyte", "Astrocyte", "Astrocyte", "Astrocyte", 
"Astrocyte", "Astrocyte", "Astrocyte", "Astrocyte", "Astrocyte", 
"Astrocyte", "Astrocyte", "Astrocyte", "Astrocyte", "Astrocyte", 
"Astrocyte", "Astrocyte", "Astrocyte", "Astrocyte", "Astrocyte", 
"Astrocyte", "Astrocyte", "Microglia", "Microglia", "Microglia", 
"Microglia", "Microglia", "Microglia", "Microglia", "Microglia", 
"Microglia", "Endothelial", "Endothelial", "Endothelial", "Endothelial", 
"OPC", "OPC", "OPC", "OPC", "Neuron", "Neuron", "Neuron", "Neuron", 
"Neuron", "Neuron", "Neuron", "Neuron", "Neuron", "Neuron", "Neuron", 
"Neuron", "Neuron", "Neuron", "Neuron", "Neuron", "Neuron", "Neuron", 
"Neuron", "Neuron", "Neuron", "Neuron", "Neuron", "Neuron", "Neuron", 
"Neuron", "Neuron"), gene = c("DAAM2", "ASPA", "MAL", "SEC14L5", 
"MAP6D1", "DPYD", "PPP1R14A", "GJB1", "FA2H", "MAG", "CDK18", 
"LGI3", "SHC4", "UGT8", "KLK6", "KCNH8", "MOBP", "LPAR1", "ADAMTS4", 
"ERMN", "OPALIN", "CLDN11", "PLEKHB1", "GSN", "GRM3", "CNP", 
"MBP", "PLP1", "SLC14A1", "GLIS3", "GLI3", "PPP1R3C", "CHRDL1", 
"CYBRD1", "CTH", "SORCS2", "ITGB4", "RNF43", "NWD1", "PAQR6", 
"C16orf89", "ALDH1L1", "TRIM66", "HGF", "CBS", "ITGA7", "SLC30A10", 
"SLC4A4", "FGFR3", "BMPR1B", "ATP13A4", "AQP4", "GPR183", "CCL4", 
"CD83", "LAPTM5", "CSF1R", "HLA-DRA", "BCL2A1", "CD14", "CCL2", 
"APOLD1", "TM4SF1", "FLT1", "A2M", "PDGFRA", "LHFPL3", "MEGF11", 
"PCDH15", "KCNK1", "KIAA1324", "LNX1", "NELL1", "COBL", "SLITRK1", 
"DPYSL5", "C14orf37", "DLX1", "DLX6", "GLRA2", "DLX2", "DLX5", 
"SLC10A4", "EGFR", "SST", "PNOC", "NXPH1", "BCL11A", "DCN", "TMEM130", 
"CNTN4", "CDO1", "NFASC", "LRRTM3", "GRIA3", "RELN")), class = "data.frame", row.names = c(NA, 
-96L))

zhu = data.frame(gene = c('SLC1A2', 'ADGRV1', 'SLC1A3', 'NKAIN3', 'GPC5', 'GPM6A', 'APOE', 
'NRXN1', 'ALDH1A1', 'MSI2', 'SLC1A2', 'SLC1A3', 'ADGRV1', 'NKAIN3', 
'ENSMMUG00000000428', 'GPC5', 'NCKAP5', 'CDH20', 'MAG', 'PLP1', 
'RNASE1', 'EPAS1', 'EBF1', 'FLT1', 'SLC7A1', 'ITM2A', 'SLC2A1', 
'GSN', 'ENSMMUG00000015144', 'MECOM', 'NAA15', 'SATB2', 'LDB2', 
'NRG1', 'ARPP21', 'RTN1', 'HECW1', 'CTIF', 'PRKCB', 'KCNQ5', 
'ENSMMUG00000044284', 'ENSMMUG00000040417', 'RPS11', 'KCNIP4', 
'ARPP21', 'KCNQ5', 'NRG1', 'STXBP5L', 'ENSMMUG00000003851', 'OPCML', 
'FAM19A1', 'DPYD', 'WFS1', 'PTPRK', 'ZBTB18', 'EPHA6', 'SERPINE2', 
'GRIA4', 'CDH12', 'CACNA2D1', 'RORB', 'IL1RAPL2', 'TSHZ2', 'POU6F2', 
'FSTL5', 'PDZRN4', 'CLSTN2', 'DCC', 'CPNE4', 'HECW1', 'KCNIP4', 
'NRG1', 'PHACTR1', 'STXBP5L', 'KCNQ5', 'NKAIN2', 'RGS6', 'LRRTM4', 
'LDB2', 'OPCML', 'SYT1', 'SNAP25', 'ENSMMUG00000028701', 'NRGN', 
'CALM2', 'CALM1', 'CAMK2A', 'NEFL', 'CHN1', 'ENSMMUG00000020894', 
'ENSMMUG00000009066', 'RPL26', 'ENSMMUG00000012537', 'KCNIP4', 
'NRG1', 'NKAIN2', 'PHACTR1', 'KCNQ5', 'ARPP21', 'LDB2', 'TLE4', 
'RXFP1', 'ASIC2', 'HTR2C', 'SEMA3E', 'EPHA5', 'SEMA3A', 'KIAA1217', 
'ENSMMUG00000039690', 'FOXP2', 'ENSMMUG00000005253', 'C14ORF93', 
'ENSMMUG00000009272', 'KCNIP4', 'KCNQ5', 'NKAIN2', 'ARPP21', 
'CAMK2A', 'SATB2', 'PHACTR1', 'SLC1A2', 'SLC1A3', 'ADGRV1', 'NKAIN3', 
'GPC5', 'APOE', 'ALDH1A1', 'SLC4A4', 'PDZRN4', 'MSI2', 'ERBB4', 
'NXPH1', 'ZNF385D', 'MYO16', 'ZNF804A', 'BTBD11', 'KCNC2', 'GAD1', 
'DPP10', 'ANK1', 'SST', 'GRIK1', 'NXPH1', 'SPOCK3', 'GRIK2', 
'KIAA1217', 'SLC24A3', 'GRIN3A', 'NRXN3', 'KIF26B', 'ADARB2', 
'ERBB4', 'VIP', 'ZNF804A', 'TAC3', 'THSD7A', 'GALNTL6', 'CRH', 
'CNTNAP2', 'SLC24A3', 'ADARB2', 'CXCL14', 'RELN', 'ERBB4', 'GRIK2', 
'GRIK1', 'CNTN5', 'GALNTL6', 'FSTL5', 'GAD2', 'FBXL7', 'ENSMMUG00000012484', 
'PTCHD4', 'SLC6A1', 'GRIK1', 'EYA4', 'NXPH1', 'CDH13', 'SGCZ', 
'ADARB2', 'PLP1', 'MBP', 'MAG', 'RNF220', 'MOG', 'FA2H', 'TF', 
'OPALIN', 'ANLN', 'ST18', 'LHFPL3', 'PTPRZ1', 'APOD', 'PCDH15', 
'COL9A1', 'SEMA5A', 'LUZP2', 'XYLT1', 'SOX6', 'VCAN', 'PTGDS', 
'SLC13A4', 'MFAP4', 'IGFBP5', 'BNC2', 'DSP', 'SLC6A20', 'VIM', 
'CPAMD8', 'SLC22A8'), type = c('Astro1', 'Astro1', 'Astro1', 'Astro1', 'Astro1', 'Astro1', 
'Astro1', 'Astro1', 'Astro1', 'Astro1', 'Astro2', 'Astro2', 'Astro2', 
'Astro2', 'Astro2', 'Astro2', 'Astro2', 'Astro2', 'Astro2', 'Astro2', 
'Endo', 'Endo', 'Endo', 'Endo', 'Endo', 'Endo', 'Endo', 'Endo', 
'Endo', 'Endo', 'ExN1', 'ExN1', 'ExN1', 'ExN1', 'ExN1', 'ExN1', 
'ExN1', 'ExN1', 'ExN1', 'ExN1', 'ExN10', 'ExN10', 'ExN10', 'ExN10', 
'ExN10', 'ExN10', 'ExN10', 'ExN10', 'ExN10', 'ExN10', 'ExN2', 
'ExN2', 'ExN2', 'ExN2', 'ExN2', 'ExN2', 'ExN2', 'ExN2', 'ExN2', 
'ExN2', 'ExN3', 'ExN3', 'ExN3', 'ExN3', 'ExN3', 'ExN3', 'ExN3', 
'ExN3', 'ExN3', 'ExN3', 'ExN4', 'ExN4', 'ExN4', 'ExN4', 'ExN4', 
'ExN4', 'ExN4', 'ExN4', 'ExN4', 'ExN4', 'ExN5', 'ExN5', 'ExN5', 
'ExN5', 'ExN5', 'ExN5', 'ExN5', 'ExN5', 'ExN5', 'ExN5', 'ExN6', 
'ExN6', 'ExN6', 'ExN6', 'ExN6', 'ExN6', 'ExN6', 'ExN6', 'ExN6', 
'ExN6', 'ExN7', 'ExN7', 'ExN7', 'ExN7', 'ExN7', 'ExN7', 'ExN7', 
'ExN7', 'ExN7', 'ExN7', 'ExN8', 'ExN8', 'ExN8', 'ExN8', 'ExN8', 
'ExN8', 'ExN8', 'ExN8', 'ExN8', 'ExN8', 'ExN9', 'ExN9', 'ExN9', 
'ExN9', 'ExN9', 'ExN9', 'ExN9', 'ExN9', 'ExN9', 'ExN9', 'InN1', 
'InN1', 'InN1', 'InN1', 'InN1', 'InN1', 'InN1', 'InN1', 'InN1', 
'InN1', 'InN2', 'InN2', 'InN2', 'InN2', 'InN2', 'InN2', 'InN2', 
'InN2', 'InN2', 'InN2', 'InN3', 'InN3', 'InN3', 'InN3', 'InN3', 
'InN3', 'InN3', 'InN3', 'InN3', 'InN3', 'InN4', 'InN4', 'InN4', 
'InN4', 'InN4', 'InN4', 'InN4', 'InN4', 'InN4', 'InN4', 'InN5', 
'InN5', 'InN5', 'InN5', 'InN5', 'InN5', 'InN5', 'InN5', 'InN5', 
'InN5', 'Oligo', 'Oligo', 'Oligo', 'Oligo', 'Oligo', 'Oligo', 
'Oligo', 'Oligo', 'Oligo', 'Oligo', 'OPC', 'OPC', 'OPC', 'OPC', 
'OPC', 'OPC', 'OPC', 'OPC', 'OPC', 'OPC', 'Peri', 'Peri', 'Peri', 
'Peri', 'Peri', 'Peri', 'Peri', 'Peri', 'Peri', 'Peri'),stringsAsFactors=FALSE)