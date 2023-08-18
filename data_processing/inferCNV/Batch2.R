library(infercnv)
library(Matrix)
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(stringr)
library(reticulate)

# Get count data
SeuratDisk::Convert('/wynton/home/ye/thomas_mazumder/data/aneuploidy/counts_Batch2.h5ad', dest = "h5seurat", overwrite = TRUE)
dat <- LoadH5Seurat('/wynton/home/ye/thomas_mazumder/data/aneuploidy/counts_Batch2.h5seurat')
counts = GetAssayData(object = dat, slot = "counts")

# Get gene data
genes <- read.table('/wynton/home/ye/thomas_mazumder/data/aneuploidy/genes_Batch2.csv', sep = ',', header = F, row.names = 1)
rownames(counts) = rownames(genes) # give gene names to counts

# Get annotations
annotations = read.table('/wynton/home/ye/thomas_mazumder/data/aneuploidy/annotations_Batch2.csv', sep = ',', header = F, row.names = 1)
colnames(counts) <- rownames(annotations)# give the column names of the count data the cell barcodes

# Get gene and chromosome info 
geneOrder = read.table('/wynton/home/ye/thomas_mazumder/data/aneuploidy/inferCNVgeneName.txt', sep = '\t', row.names = 1)

# Create inferCNV object
infercnv_obj <- CreateInfercnvObject(
raw_counts_matrix=counts,
annotations_file=annotations,
gene_order_file=geneOrder,
ref_group_names = c('CONTROL')
)

# Run inferCNV
infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir='/wynton/home/ye/thomas_mazumder/data/aneuploidy/results/Batch2/', 
    cluster_by_groups=TRUE, 
    denoise=TRUE,
    HMM=TRUE,
    num_threads=10,
)

cnvRes <- infercnv_obj@expr.data
# write results
write.table(
cnvRes,
file = stringr::str_interp('/wynton/home/ye/thomas_mazumder/data/aneuploidy/results/Batch2/results.txt'),
sep = '\t',
row.names = T,
col.names = T,
quote = F)
