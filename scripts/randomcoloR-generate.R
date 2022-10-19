#!/usr/bin/env Rscript
# 🐛

arguments = commandArgs(trailing=TRUE)

k = as.integer(arguments[1])

# Grab the current color space
currentColorSpace = randomcoloR:::ourColorSpace@coords

# Perform UMAP
umapColorSpace = umap::umap(currentColorSpace)$layout

colors = unname(colorspace::hex(colorspace::LAB(matrix(currentColorSpace[cluster::pam(umapColorSpace,k)$id.med,],nrow=k))))

dir.create('palettes',showWarnings=FALSE)

write(paste(colors,collapse='\t'),file=file.path('palettes',paste0('randomcoloR-umap-',formatC(k,width=4,flag=0),'.txt')))
