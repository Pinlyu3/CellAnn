#########
#########
## pre-calculated data before the database ##
#########
#########

library(Seurat)
library(Matrix)

######### test /zp1/data/plyu3/SCPData #######
######### Rscript --vanilla sillyScript.R iris.txt out.txt ######

######### folder = '/zp1/data/plyu3/SCPData'
######### file = 'pmid32994417_seurat_obj'

Args <- commandArgs()

folder <- Args[6]
file <- Args[7]

setwd(folder)

seurat_obj = readRDS(file)

tag_pre = strsplit(file,split='_seurat_obj')
tag = tag_pre[[1]]

print(tag)

Output_dims <- function(seurat_obj,tag){
	######
	meta = seurat_obj@meta.data
	######
	###### key words: celltype,dim1,dim2,dim3 ########
	k = which(colnames(meta) %in% c('celltype','dim1','dim2','dim3') == T)
	meta_cl = meta[,k]
	###### change celltype to cluster ######
	k = which(colnames(meta_cl) == 'celltype')
	colnames(meta_cl)[k] = 'cluster'
	###### saveRDS #####
	FN = paste(tag,'_Dimplot',sep='')
	saveRDS(meta_cl,FN)
}

Output_dims(seurat_obj,tag)

Prepare_the_Seurat_objects_Step2 <- function(seurat_obj,tag){
	Idents(seurat_obj) = 'celltype'
	######
	###### calculate DEGs between cell types ##########
	celltype = seurat_obj$celltype
	celltype = as.character(celltype[!duplicated(celltype)])
	######
	DEGs_list = list()
	######
	for(i in 1:(length(celltype)-1)){
		for(j in (i+1):length(celltype)){
			print(paste(celltype[i],celltype[j],sep='@VS@'))
			library(future)
			plan("multicore", workers = 30)
			tmp_markers = FindMarkers(seurat_obj,ident.1=celltype[i],ident.2=celltype[j],test.use='MAST',logfc.threshold = 0.5)
			tmp_markers = data.frame(tmp_markers)
			tmp_markers$gene = rownames(tmp_markers)
			k = which(tmp_markers$p_val_adj < 0.01)
			tmp_markers = tmp_markers[k,]
			tmp_markers = list(tmp_markers)
			names(tmp_markers) = paste(celltype[i],celltype[j],sep='@VS@')
			DEGs_list = c(DEGs_list,tmp_markers)
		}
	}
	FN = paste(tag,'DEGs_long',sep='_')
	saveRDS(DEGs_list,file=FN)
}





