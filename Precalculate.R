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
######### nohup Rscript /zp1/data/plyu3/SCPData/precode/Precalculate.R /zp1/data/plyu3/SCPData pmid32994417_seurat_obj > pmid32994417_info &

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


#### OK next: #########
#### using the data matrix in seurat: ############

check_the_matrix <- function(seurat_obj){
	seurat_obj_mat = seurat_obj[['RNA']]@data
	seurat_obj_mat = as.vector(seurat_obj_mat)
	max = max(seurat_obj_mat)
	######
	max_biggest = log2(1e6+1)
	if(max < max_biggest){
		print('data matrix looks good!')
	}else{
		print('data matrix has something wrong')
	}
}

check_the_matrix(seurat_obj)


#### Next: get bulk data for each cell type ####

get_NonOverlap <- function(knnIdx,Sample_Numbers,non_overlapping_cutoff = 0.99){
	counts = 0
	index_list_ori = sample(1:dim(knnIdx)[1],1)
	i = 1
	while(i < (Sample_Numbers)){
		######
		tmp_index = sample(1:dim(knnIdx)[1],1)
		######
		counts = counts+1
		#####
		tmp_index_new_dat = as.vector(knnIdx[tmp_index,])
		tmp_index_old_dat = as.vector(knnIdx[index_list_ori,])
		######
		tmp_overlap = length(which(tmp_index_new_dat %in% tmp_index_old_dat == T)) / length(tmp_index_new_dat)
		######
		if(tmp_overlap < (1-non_overlapping_cutoff)){
			index_list_ori = c(index_list_ori,tmp_index)
			i = i+1
		}
		if(counts > 10000){
			break
		}
	}
	print(paste('clusters = ',length(index_list_ori)))
	return(index_list_ori)
}


get_smooth_gene_exp <- function(x,mat){
	#####
	mat_cl = mat[,x]
	#####
	vector = rowMeans(mat_cl)
	#####
	return(vector)
}


Prepare_the_Seurat_objects_Step3 <- function(seurat_obj,matrix_tag = 'data',KNN = 5, non_overlapping_cutoff = 0.99, tag='pmid30348985'){
	seurat_obj$pesudo_index = 'Non_Selected'
	################
	if(matrix_tag == 'data'){
		mat = seurat_obj[['RNA']]@data
	}
	if(matrix_tag == 'counts'){
		mat = seurat_obj[['RNA']]@counts
	}
	Idents(seurat_obj) = 'celltype'
	###### added the pseudo-bulk to the single-cell RNA-seq datasets #######
	###### calculate DEGs between cell types ##########
	celltype = seurat_obj$celltype
	celltype = as.character(celltype[!duplicated(celltype)])
	dims_tab = seurat_obj@meta.data
	k = which(colnames(dims_tab) %in% c('dim1','dim2','dim3'))
	dims_tab = dims_tab[,k]
	######
	###### using UMAPs as dim #################
	All_mat = list()
	######
	for(i in 1:(length(celltype))){
		#print(i)
		#############
		tmp_k = which(seurat_obj$celltype == celltype[i])
		#print(length(tmp_k))
		#####
		tmp_mat = mat[,tmp_k]
		#####
		tmp_dims_tab = dims_tab[tmp_k,]
		tmp_dims_tab = as.matrix(tmp_dims_tab)
		############# decide KNN ##########
		print(paste('KNN =',KNN))
		#############
		#############
		knnObj <- FNN::get.knn(data = tmp_dims_tab,k = KNN)
		knnIdx <- knnObj$nn.index
		rownames(knnIdx) = 1:dim(knnIdx)[1]
		############# sample cells #########################
		Sample_Numbers = round(length(tmp_k)/(KNN*5))
		if(Sample_Numbers > 300){Sample_Numbers=300}
		if(Sample_Numbers < 5){Sample_Numbers=5}
		#print(paste('Sample_Numbers',Sample_Numbers))
		#############
		Sample_index = get_NonOverlap(knnIdx,Sample_Numbers,non_overlapping_cutoff)
		#############
		sample_index_mat = knnIdx[Sample_index,]
		rownames(sample_index_mat) = paste(celltype[i],1:dim(sample_index_mat)[1],sep='@psedobulk')
		############# problems with index ###############
		#############
		tmp_dims_tab = data.frame(tmp_dims_tab)
		tmp_dims_tab$pesudo_index = 'Non_Selected'
		for(j in 1:dim(sample_index_mat)[1]){
			tmp_k = sample_index_mat[j,]
			tmp_dims_tab$pesudo_index[tmp_k] = rownames(sample_index_mat)[j]
			tmp_k_total = match(rownames(tmp_dims_tab),colnames(seurat_obj))
			# print(head(tmp_dims_tab$pesudo_index))
			seurat_obj$pesudo_index[tmp_k_total] = tmp_dims_tab$pesudo_index
		}
		#############
		sample_index_mat_res = apply(sample_index_mat,1,get_smooth_gene_exp,mat=tmp_mat)
		sample_index_mat_res = as.matrix(sample_index_mat_res,nrow=dim(tmp_mat)[1])
		#############
		rownames(sample_index_mat_res) = rownames(tmp_mat)
		colnames(sample_index_mat_res) = paste(celltype[i],1:dim(sample_index_mat_res)[2],sep='@pseudobulk')
		#############
		sample_index_mat_res = list(sample_index_mat_res)
		#############
		All_mat = c(All_mat,sample_index_mat_res)
	}
	##### 
	#####
	All_mat = do.call('cbind',All_mat)
	All_mat = round(All_mat,3)
	##### save the pseudo-matrix #########
	FN_mat = paste(tag,'pseudo_bulk',sep='_')
	save(All_mat,file=FN_mat)
	##### save the seuratobj #############
	FN = paste(tag,'seurat_obj',sep='_')
	saveRDS(seurat_obj,file=FN)
	#####
	print('Done!')
}

####
Prepare_the_Seurat_objects_Step3(seurat_obj,matrix_tag = 'data',KNN = 5, non_overlapping_cutoff = 0.99, tag=tag)


#### Next calculate the total average expression by cell types ####



Prepare_the_Seurat_objects_Step1 <- function(seurat_obj,matrix_tag = 'data',tag='pmid32386599'){
	####### average expression ########
	####### generate average expression for each celltype #########
	####### key: celltype ######
	if(matrix_tag == 'data'){
		mat = seurat_obj[['RNA']]@data
	}
	if(matrix_tag == 'counts'){
		mat = seurat_obj[['RNA']]@counts
	}
	####### 
	celltype = seurat_obj$celltype
	celltype = as.character(celltype[!duplicated(celltype)])
	#######
	avg_total = c()
	for(i in 1:length(celltype)){
		#######
		print(celltype[i])
		k = which(seurat_obj$celltype == celltype[i])
		#######
		mat_cl = mat[,k]
		#######
		avg = Matrix::rowMeans(mat_cl)
		avg_total = c(avg_total,avg)
	}
	avg_total = matrix(avg_total,ncol=length(celltype))
	colnames(avg_total) = celltype
	rownames(avg_total) = rownames(seurat_obj)
	k = which(rowSums(avg_total) == 0)
	if(length(k) > 0){avg_total = avg_total[-k,]}
	print(dim(avg_total))
	#####
	FN = paste(tag,'_avg_expmat',sep='')
	saveRDS(avg_total,FN)
}



Prepare_the_Seurat_objects_Step1(seurat_obj,tag = tag)


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


Prepare_the_Seurat_objects_Step2(seurat_obj,tag = tag)



