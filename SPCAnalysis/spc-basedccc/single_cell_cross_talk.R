
parser = argparse::ArgumentParser(description = "Script to QC and Cluster scRNA data")
parser$add_argument('-s', '--sample', help = 'sample')
opts = parser$parse_args()

#load library
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(dplyr)
library(tidyr)
library(reshape2)
library(NICHES)
library(stringr)

source("./NICHES_3D.R")

custom_lrdb <- read.table('./CellChatDB.human.interaction.planarian.txt', sep='\t', header=T)
custom_lrdb <- custom_lrdb[c('gene_id.ligand', 'gene_id.receptor')]
custom_lrdb <- na.omit(custom_lrdb)
custom_lrdb <- custom_lrdb[!duplicated(custom_lrdb), ]
custom_lrdb

rds.result <- paste0(opts$sample, '.rds')
if (!file.exists(rds.result)){
    obj <- readRDS('/dellfsqd2/ST_OCEAN/USER/hankai/Project/10.smed/37.Nb2_niche/NsC.rds')
    obj <- subset(obj, timepoint == opts$sample)
    results <- RunNICHES(obj,
                     assay="SCT",
	                 LR.database="custom",
	                 min.cells.per.ident = 3,
	                 min.cells.per.gene = 3,
	                 meta.data.to.map = colnames(obj@meta.data),
	                 position.x = 'x',
	                 position.y = 'y',
	                 position.z = 'z',
	                 cell_types = 'hood_inte.seurat_clusters',
	                 custom_LR_database = custom_lrdb,
	                 k = NULL,
	                 rad.set = 60, # median_mindis * 5, five layer cells
	                 blend = 'mean',
	                 CellToCell = F,
	                 CellToSystem = F,
	                 SystemToCell = F,
	                 CellToCellSpatial = T,
	                 CellToNeighborhood = F,
	                 NeighborhoodToCell = F,
	                 output_format = "seurat"
	                 )

    results
    obj.result <- results$CellToCellSpatial
    saveRDS(obj.result, rds.result)
}else {
    obj.result <- readRDS(rds.result)
}
obj.loom <- as.loom(obj.result, filename = paste0(opts$sample, '.loom'), overwrite = T)
obj.loom$close_all()

