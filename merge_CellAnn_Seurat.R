#####
##### 
#####

merge_CellAnn_Seurat <- function(seurat_obj,cluster,cell_anno){
	library(Seurat)
	library(readr)
	####
	cell_anno_table = read_delim(cell_anno,delim="\t")
	####
	Meta = seurat_obj@meta.data
	k2 = which(colnames(Meta) == cluster)
	data_cluster = Meta[,k2]
	####
	m = match(data_cluster,cell_anno_table$cluster)
	cell_anno_table_merge = cell_anno_table[m,]
	####
	seurat_obj@meta.data = cbind(seurat_obj@meta.data,cell_anno_table_merge)
	####
	return(seurat_obj)
}
