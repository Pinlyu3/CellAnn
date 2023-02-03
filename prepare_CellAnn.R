#####
##### prepare_CellAnn(seurat_obj,folder=getwd(),"test1","RNA","umap","seurat_clusters")
#####

prepare_CellAnn <- function(seurat_obj,folder=getwd(),sample_name='Human_pancreas',matrix_name='RNA',dims='UMAP',cluster='seurat_clusters'){
	###########
	library(Seurat)
	###########
	print('Start!')
	###########
	setwd(folder)
	###########
	Output_Step1_name = paste0(sample_name,'_CellAnn_Step1_input.txt')
	Output_Step4_name = paste0(sample_name,'_CellAnn_Step4_input.txt')
	###########
	AssayNames = names(seurat_obj@assays)
	k = which(AssayNames == matrix_name)
	if(length(k) == 1){
		CountMatrix = seurat_obj@assays[[k]]@counts
	}else{
		print(paste0('Error: Not find the assay: ',matrix_name,' in your seurat object'))
		return("check your seurat object")
	}
	###########
	Meta = seurat_obj@meta.data
	k2 = which(colnames(Meta) == cluster)
	data_cluster = Meta[,k2]
	###########
	AverageMatrix = prepare_CellAnn_Avg_Mat(CountMatrix,data_cluster)
	########### add column of Gene symbol ##############
	Gene = rownames(AverageMatrix)
	###########
	Out_matrix = cbind(data.frame(GENE=Gene),data.frame(AverageMatrix))
	###########
	colnames(Out_matrix) <- c('GENE',colnames(AverageMatrix))
	###########
	########### write table for step1 ################## 
	########### columns: GENE, cluster names ###########
	write.table(Out_matrix,file=Output_Step1_name,sep='\t',quote=F,row.names=F)
	########### write table for step4 ##################
	########### columns: cell, cluster, dim1, dim2 #####
    DimNames = names(seurat_obj@reductions)
    k3 = which(DimNames == dims)
    if(length(k3) == 1){
		DimMatrix = seurat_obj@reductions[[k3]]@cell.embeddings
	}else{
		print(paste0('Error: Not find the dims: ',dims,' in your seurat object'))
		return("check your seurat object")
	}
	###########
	cell = rownames(DimMatrix)
	dim_table = data.frame(DimMatrix)
	colnames(dim_table) <- paste0("dim",1:dim(dim_table)[2])
	###########
	dim_table_merge = cbind(data.frame(cell=cell),dim_table)
	########### then we add clusters ##########
	m = match(dim_table_merge$cell,rownames(Meta))
	dim_table_merge$cluster = Meta[m,k2]
	############
	write.table(dim_table_merge,file=Output_Step4_name,sep='\t',quote=F,row.names=F)
	############
	print("Done!")
	print(paste(Output_Step1_name,"are in",folder))
	print(paste(Output_Step4_name,"are in",folder))
	###########
}


####
####

prepare_CellAnn_Avg_Mat <- function(data_mat,data_cluster,scale_factor=10000){
	######
	tag_cluster = unname(data_cluster)
	tag_cluster_level = levels(as.factor(tag_cluster))
	###### normalized back datasets ######
	data_mat_exp = data_mat
	######
	###### print(paste('Sums:',head(colSums(data_mat_exp[,c(1:2)]))))
	###### data_mat_exp is 1e5 normalize #######
	merge_mat = c()
	for(i in 1:length(tag_cluster_level)){
		index = which(data_cluster %in% tag_cluster_level[i] == T)
		index_mat = data_mat_exp[,index]
		######
		index_sum = rowSums(index_mat)
		######
		merge_mat = c(merge_mat,index_sum)
	}
	###
	merge_mat = matrix(merge_mat,nrow=dim(data_mat)[1])
	###
	rownames(merge_mat) = rownames(data_mat)
	colnames(merge_mat) = tag_cluster_level
	### colSums(merge_mat)
	scale = colSums(merge_mat)/scale_factor
	merge_mat = sweep(merge_mat,2,scale,FUN='/')
	### default norm ####
	merge_mat = round(log(merge_mat+1),3)
	return(merge_mat)
	### return is a log transformed matrix ####
}