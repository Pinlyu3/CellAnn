######
######
######
######
######
Example datasets prepare 2 
######
######
######
######
######


Res_mat_highest_celltype <- function(res_mat){
	res_list = list()
	for(i in 1:dim(res_mat)[1]){
		res_mat_tmp = res_mat[i,]
		k = which(res_mat_tmp == max(res_mat_tmp))
		max_cor = res_mat_tmp[k]
		##### #######
		k2 = which(res_mat_tmp <= max_cor)
		##### #######
		res_mat_tmp_k2 = res_mat_tmp[k2]
		res_mat_tmp_k2 = sort(res_mat_tmp_k2,decreasing=T)
		#####
		if(length(res_mat_tmp_k2) >3){
			res_mat_tmp_k2 = res_mat_tmp_k2[1:3]
		}
		#####
		res_list = c(res_list,list(names(res_mat_tmp_k2)))
	}
	names(res_list) = rownames(res_mat)
	return(res_list)
}


cor_res = cor_res

collapse_cor = function(cor_res){
	######
	cols = colnames(cor_res)
	######
	cols_index = strsplit(cols,split='@')
	cols_index = sapply(cols_index,function(x) x[[1]])
	cols_index_levels = cols_index[!duplicated(cols_index)]
	######
	mat_list = c()
	for(i in cols_index_levels){
		k = which(cols_index == i)
		if(length(k) == 1){
			tmp_mat = cor_res[,k]
		}
		if(length(k) > 1){
			tmp_mat = cor_res[,which(cols_index == i)]
			tmp_mat = apply(tmp_mat,1,median)
		}
		mat_list = c(mat_list,tmp_mat)
	}
	mat_list = matrix(mat_list,ncol=length(cols_index_levels))
	colnames(mat_list) = cols_index_levels
	rownames(mat_list) = rownames(cor_res)
	return(mat_list)
}

######
###### load datasets using DEGs to plot correlations ######
######
i = 25

Main_compare_process_New2 <- function(compare_df,method,folder=folder,mode='easy'){
	library(Seurat)
	##########
	out_table_list = list()
	##########
	for(i in 1:dim(compare_df)[1]){
		print(paste("num",i))
		print(paste('query:',compare_df[i,1],'  ','ref:',compare_df[i,2],sep=''))
		TAG = paste('query:',compare_df[i,1],'-->','ref:',compare_df[i,2],sep='')
		#######
		query = compare_df[i,1]
		ref = compare_df[i,2]
		query_index = paste0(query,'_test_input')
		ref_index = paste0(ref,'_Ref')
		####### load ref and query seurat object !!!! #######
		query_seurat = loadRData(query_index)
		ref_seurat = loadRData(ref_index)
		#######	OK!
		####### we next find how is the overlap between ref and query datasets ##########
		query_seurat_ct = levels(as.factor(query_seurat$celltype))
		ref_seurat_ct = levels(as.factor(ref_seurat$celltype))
		####### then we see how many overlap cell types between query and ref ###########
		overlap_ct = query_seurat_ct[which(query_seurat_ct %in% ref_seurat_ct == T)]
		#######
		if(length(overlap_ct) == 0){
			print('Warning: no overlaps between ref and query datasets !')
		}
		####### then we load query input matrix !!!! #####
		query_avg = paste(query,'_test_input_ClusterAvg.txt',sep='')
		query_mat = read.table(query_avg,header=T)
		query_mat = df_to_mat(query_mat)
		####### ref input mat :: ref_mat #######
		####### ref input cell types:: ref_label
		ref_mat = ref_seurat[['RNA']]@data
		ref_label = unname(ref_seurat$celltype)
		####### then we will filter query datasets if the mode == 'easy' !!!! ####
		if(mode == 'easy'){
			#### we will filter the query matrix !!! #######
			#### we need to load the GroundTruth to get the cell types of the query matrix ########
			query_cluster_truth_index = paste(query,'_test_input_GroundTruthLabel.txt',sep='')
			query_cluster_truth = read.table(query_cluster_truth_index,sep='\t',header=T)
			#### clusters in query must be the overlapped ct #####
			query_seurat_cluster = query_cluster_truth$cluster[which(query_cluster_truth$ground.truth %in% overlap_ct == T)]
			#### Then we update the query matrix !!!! ####
			query_mat = query_mat[,which(colnames(query_mat) %in% query_seurat_cluster == T)]
		}
		########
		if(dim(query_mat)[2] < 2){
			print('Warning:: the query matrix must have at least 2 clusters !!!! Skip !!!!')
    		next
		}
		########## load DEGs #######
		ref = compare_df[i,2]
		ref_DEGs_names = paste(ref,'_Ref_CTDEG_ORDER150',sep='')
		ref_DEGs_input = readRDS(ref_DEGs_names)
		##########
		ref_DEGs_names = paste(ref,'_Ref_CTsubDEG_ORDER150',sep='')
		ref_subDEGs_input = readRDS(ref_DEGs_names)
		########## load ref ########	
		ref_mat_input_name = paste(ref,'_Ref_SubCTAvg',sep='')
		ref_mat_input = readRDS(ref_mat_input_name)
		query_mat_input = query_mat
		##########
		cor_res = Out_Cor2(ref_mat_input,query_mat_input,ref_DEGs_input)
		#hist(cor_res,breaks = 50)
		cor_res2 = collapse_cor(cor_res)
		##### see the positive and negative samples 
		res_mat_v = as.vector(cor_res2)
		res_mat_v = sort(res_mat_v,decreasing=T)
		model <- mclust::densityMclust(as.vector(res_mat_v),G=1:5)
		#####
		cutoffs = Get_model_cutoff1_simple(model)
		candidate_align = Res_mat_highest_celltype(cor_res)
		#####
		res_max = apply(cor_res,1,max)
		print(paste('res_max',min(res_max)))
		print(paste('cutoff',cutoffs))
		Unassigned_index = which(res_max < cutoffs)
		#####
		res_tab = as.character(sapply(candidate_align,function(x) x[[1]]))
			res_tab = gsub('@(.+)','',res_tab)
			for(j in 1:length(candidate_align)){
				tmp = candidate_align[[j]]
				tmp_index = sapply(strsplit(tmp,split='@'),function(x) x[[1]])
				if(length(levels(as.factor(tmp_index))) == 1){
					next
				}else{
					print('compare!')
					ct_list = list()
					for(ii in 1:length(tmp)){
						sub_ct = tmp[ii]
						subDEGs = ref_subDEGs_input
						sub_ct_DEGs_index = which(colnames(subDEGs) == sub_ct)
						sub_ct_DEGs = subDEGs[,c(sub_ct_DEGs_index,sub_ct_DEGs_index+1)]
						sub_ct_DEGs = sub_ct_DEGs[1:20,1]
						sub_query_mat_input = query_mat_input[,names(candidate_align)[j]]
						query_sub_ct = sub_query_mat_input[which(names(sub_query_mat_input) %in% sub_ct_DEGs == T)]
						ct_list = c(ct_list,list(query_sub_ct))
					}
					#### then get DEGs #####
					indexJ = Compare3_corr(ct_list,tmp_index)
					res_tab[j] = indexJ
				}
			}
			#####
			#####
			if(length(Unassigned_index) > 0){
				res_tab[Unassigned_index] = 'unassigned'
			}
			if(length(Unassigned_index) > 0){
				res_tab[Unassigned_index] = 'unassigned'
			}
			res_table = data.frame(cluster=colnames(query_mat_input),result=res_tab)
			res_table = list(res_table)
			names(res_table) = TAG
			out_table_list = c(out_table_list,res_table)
		}
		return(out_table_list)
}



folder = "C:/Users/plyu3/Desktop/CellAnn_methods_test/Tabula_Muris_mouse_data_prepare"
Mouse_self_res = Main_compare_process_New2(Self_compare_list,method='CellAnnV5',folder=folder,mode='easy')
Mouse_self_res_v = Visualize_res(Mouse_self_res)
Mouse_self_res_v_p = Visualize_res_plot(Mouse_self_res_v)
Mouse_self_res_v_p = Visualize_res_plot_process(Mouse_self_res_v_p)

ggplot(Mouse_self_res_v_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')

folder = "C:/Users/plyu3/Desktop/CellAnn_methods_test/Tabula_Muris_mouse_data_prepare"
Mouse_self_res = Main_compare_process_New2(mouse_df_batch,method='CellAnnV5',folder=folder,mode='easy')
Mouse_self_res_v = Visualize_res(Mouse_self_res)
Mouse_self_res_v_p = Visualize_res_plot(Mouse_self_res_v)
Mouse_self_res_v_p = Visualize_res_plot_process(Mouse_self_res_v_p)

ggplot(Mouse_self_res_v_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')



Selected_DEGs <- function(ref_DEGs_input){
	dims = dim(ref_DEGs_input)[2]/2
	####
	all_DEGs = c()
	for(i in 1:dims){
		tmp_dat = ref_DEGs_input[,c(i*2-1,i*2)]
		tmp_dat = tmp_dat[order(tmp_dat[,2],decreasing=T),]
		DEGs_1 = tmp_dat[1:20,]
		DEGs_2 = tmp_dat[dim(tmp_dat)[1]:(dim(tmp_dat)[1]-19),]
		all_DEGs = c(all_DEGs,DEGs_1[,1],DEGs_2[,1])
	}
	all_DEGs = all_DEGs[!duplicated(all_DEGs)]
	return(all_DEGs)
}


#apply(merge_mat,2,function(x) return(length(which(x == 0))))


#x = corr_tab_cl[9,]
#vector1 = query_mat_input_cl[,5]
#vector2 = ref_mat_input_cl[,10]
#cor(vector1,vector2)
#tab = data.frame(v1=vector1,v2=vector2)
#k = which(vector1 < 0.1 & vector2 < 0.1)


#calculate_tab <- function(x,query_mat_input_cl,ref_mat_input_cl){
#	vector1 = query_mat_input_cl[,which(colnames(query_mat_input_cl) == x[1])]
#	vector2 = ref_mat_input_cl[,which(colnames(ref_mat_input_cl) == x[2])]
#	k = which(vector1 == 0 & vector2 == 0)
#	if(length(k) > 0){
#		vector1 = vector1[-k]
#		vector2 = vector2[-k]
#	}
#	res = cor(vector1,vector2)
#	return(res)
#}


calculate_Cor2 <- function(query_mat_input,ref_mat_input,DEGs_overlap){
	######
	query_mat_input_cl = query_mat_input[which(rownames(query_mat_input) %in% DEGs_overlap == T),]
	ref_mat_input_cl = ref_mat_input[which(rownames(ref_mat_input) %in% DEGs_overlap == T),]
	######
	m1 = match(DEGs_overlap,rownames(query_mat_input_cl))
	m2 = match(DEGs_overlap,rownames(ref_mat_input_cl))
	query_mat_input_cl= query_mat_input_cl[m1,]
	ref_mat_input_cl= ref_mat_input_cl[m2,]
	######
	merge_mat = cbind(query_mat_input_cl,ref_mat_input_cl)
	corr_res0 = cor(merge_mat)
	query_dim = dim(query_mat_input_cl)[2]
	corr_res0 = corr_res0[,-c(1:query_dim)]
	corr_res0 = corr_res0[c(1:query_dim),]
	######
	#query_row = data.frame(query=colnames(query_mat_input_cl))
	#ref_row = data.frame(ref=colnames(ref_mat_input_cl))
	#corr_tab = merge(query_row,ref_row)
	######
	#corr_res = apply(corr_tab,1,calculate_tab,query_mat_input_cl,ref_mat_input_cl)
	#corr_tab$corr = corr_res
	### corr_tab[which(corr_tab$query=='C13_1'),]
	### corr_tab_cl = corr_tab[which(corr_tab$query=='C13_1'),]
	#Cor_res <- cor(merge_mat, method = 'spearman')
	###### split the Cor_res #######
	#query_dim = dim(query_mat_input_cl)[2]
	#ref_dim = dim(ref_mat_input_cl)[2]
	######
	#Cor_res = Cor_res[,-c(1:query_dim)]
	#Cor_res = Cor_res[c(1:query_dim),]
	###### Then we output the most largest clusters ########
	return(corr_res0)
}



Out_Cor2 <- function(ref_mat_input,query_mat_input,ref_DEGs_input){
	#######
	genes_overlap = rownames(ref_mat_input)[which(rownames(ref_mat_input) %in% rownames(query_mat_input) == T)]
	#######
	m = match(genes_overlap,rownames(ref_mat_input))
	ref_mat_input2 = ref_mat_input[m,]
	#######
	m = match(genes_overlap,rownames(query_mat_input))
	query_mat_input2 = query_mat_input[m,]
	#######
	ref_mat_input3 = exp(ref_mat_input2)-1
	print(colSums(ref_mat_input3))
	query_mat_input3 = exp(query_mat_input2)-1
	print(colSums(query_mat_input3))
	#######
	combined_mat = cbind(query_mat_input3,ref_mat_input3)
	combined_mat_norm = limma::normalizeBetweenArrays(combined_mat,method="quantile")
	combined_mat_norm = log(combined_mat_norm+1)
  	#####
  	combined_mat_norm_rank = matrix_to_ranks(combined_mat_norm)
	####### selected DEGs #######
	DEGs_overlap = Selected_DEGs(ref_DEGs_input)
	#######
	query_mat_input4 = combined_mat_norm_rank[,1:dim(query_mat_input3)[2]]
	ref_mat_input4 = combined_mat_norm_rank[,-c(1:dim(query_mat_input3)[2])]
	#######
	cor_res = calculate_Cor2(query_mat_input4,ref_mat_input4,DEGs_overlap)
	####### return format is a vector of celltype names ######
	return(cor_res)
}

#######
#######

#print(apply(cor_res,1,max))
#test_v = apply(cor_res,1,max)

#no_Zero_genes = apply(query_mat_input4,2,function(x) length(which(x>0)))

#fold = no_Zero_genes / max(no_Zero_genes)

#norm(test_v,fold)

############
############
#query_mat_input_test = apply(query_mat_input4,2,function(x) return(length(which(x > 0))))

#query_mat_input_test = query_mat_input_test / max(query_mat_input_test)
############
############

cor_res = Out_Cor2(ref_mat_input,query_mat_input,ref_DEGs_input)
hist(cor_res,breaks = 50)

##### see the positive and negative samples 
res_mat_v = as.vector(cor_res)
res_mat_v = sort(res_mat_v,decreasing=T)
model <- mclust::densityMclust(as.vector(res_mat_v),G=1:5)

Get_model_cutoff1_simple(model)
print(apply(cor_res,1,max))

### cor_res[which(rownames(cor_res) == 'C13_1'),]

Get_model_cutoff1_simple <- function(model){
	########
	index1 = 1
	model1 = model$classification[index1]
	index2 = length(which(model$classification == model1)) + 1
	model2 = model$classification[index2]
	########
	print(paste("length1 = ",length(which(model$classification == model1))))
	cutoff1 = model$data[length(which(model$classification == model1))]
	cutoff2 = model$data[length(which(model$classification == model1)) +1]
	## return((cutoff1+cutoff2)/2)
	cutoff = cutoff1
	######## if model1 and model2 are too overlap ######
	model_mean_total = model$parameters$mean
	model_sd_total = model$parameters$variance$scale
	if(length(model_sd_total) == 1){
		model_No1_sd = model_sd_total
		model_No2_sd = model_sd_total
	}
	if(length(model_sd_total) > 1){
		model_No1_sd = model_sd_total[model1]
		model_No2_sd = model_sd_total[model2]
	}
	########
	model_No1_mean = model_mean_total[model1]
	model_No2_mean = model_mean_total[model2]
	########
	model_1_overlap = pnorm(cutoff,mean=model_No1_mean,sd=sqrt(model_No1_sd))
	model_2_overlap = 1-pnorm(cutoff,mean=model_No2_mean,sd=sqrt(model_No2_sd))
	#######
	model_1_overlap_cut = qnorm(0.1,mean=model_No1_mean,sd=sqrt(model_No1_sd))
	model_2_overlap_cut = qnorm(0.9,mean=model_No2_mean,sd=sqrt(model_No2_sd))
	if(model_1_overlap+model_2_overlap > 0.1){
		cutoff = qnorm(0.5,mean=model_No2_mean,sd=sqrt(model_No2_sd))
	}
	####
	####
	return(cutoff)
}


x <- seq(-1, 1, by = .01)
y <- dnorm(x, mean = model_No1_mean, sd = sqrt(model_No1_sd))
plot(x,y)

#####
#####
#####


library(reshape2)
cor_res_m = melt(cor_res)

cor_res_m$index2 = sapply(strsplit(as.character(cor_res_m$Var2),split='@'),function(x) x[[1]])
##### load the ground truth !!!!! #####
query_cluster_truth_index = paste(query,'_test_input_GroundTruthLabel.txt',sep='')
query_cluster_truth = read.table(query_cluster_truth_index,sep='\t',header=T)
m = match(cor_res_m$Var1,query_cluster_truth$cluster)
cor_res_m$index1 = query_cluster_truth$ground.truth[m]

k = which(cor_res_m$index1 == cor_res_m$index2)

cor_res_m$class = 'NO'
cor_res_m$class[k] = 'Yes'

library(ggplot2)
ggplot(cor_res_m,aes(x=value,fill=class)) + geom_histogram() + facet_grid(~class)

#####
##### 没有必要去改 cutoff ！！！ #######
#####

##### red ######
##### red ######
##### red ######
##### red ######
##### OK!!! first we write functions for other tools !!!!! ######
##### we use the mouse datasets first !!!!!! ####################
#####

##### we first send the folder and the mode #####
#####

#####
setwd("C:/Users/plyu3/Desktop/CellAnn_methods_test/Tabula_Muris_mouse_data_prepare")

files = list.files()

load("Self_compare_list")

folder = "C:/Users/plyu3/Desktop/CellAnn_methods_test/Tabula_Muris_mouse_data_prepare"

compare_df = Self_compare_list
compare_df = 


Main_compare_process_OtherTools <- function(compare_df,method='scmap-cluster',folder=folder,mode='easy'){
	library(Seurat)
	##########
	out_table_list = list()
	ref_label_list = list()
	##########
	for(i in 1:dim(compare_df)[1]){
		print(paste("num",i))
		print(paste('query:',compare_df[i,1],'  ','ref:',compare_df[i,2],sep=''))
		TAG = paste('query:',compare_df[i,1],'-->','ref:',compare_df[i,2],sep='')
		#######
		query = compare_df[i,1]
		ref = compare_df[i,2]
		query_index = paste0(query,'_test_input')
		ref_index = paste0(ref,'_Ref')
		####### load ref and query seurat object !!!! #######
		query_seurat = loadRData(query_index)
		ref_seurat = loadRData(ref_index)
		#######	OK!
		####### we next find how is the overlap between ref and query datasets ##########
		query_seurat_ct = levels(as.factor(query_seurat$celltype))
		ref_seurat_ct = levels(as.factor(ref_seurat$celltype))
		####### then we see how many overlap cell types between query and ref ###########
		overlap_ct = query_seurat_ct[which(query_seurat_ct %in% ref_seurat_ct == T)]
		#######
		if(length(overlap_ct) == 0){
			print('Warning: no overlaps between ref and query datasets !')
		}
		####### then we load query input matrix !!!! #####
		query_avg = paste(query,'_test_input_ClusterAvg.txt',sep='')
		query_mat = read.table(query_avg,header=T)
		query_mat = df_to_mat(query_mat)
		####### ref input mat :: ref_mat #######
		####### ref input cell types:: ref_label
		ref_mat = ref_seurat[['RNA']]@data
		ref_label = unname(ref_seurat$celltype)
		####### then we will filter query datasets if the mode == 'easy' !!!! ####
		if(mode == 'easy'){
			#### we will filter the query matrix !!! #######
			#### we need to load the GroundTruth to get the cell types of the query matrix ########
			query_cluster_truth_index = paste(query,'_test_input_GroundTruthLabel.txt',sep='')
			query_cluster_truth = read.table(query_cluster_truth_index,sep='\t',header=T)
			#### clusters in query must be the overlapped ct #####
			query_seurat_cluster = query_cluster_truth$cluster[which(query_cluster_truth$ground.truth %in% overlap_ct == T)]
			#### Then we update the query matrix !!!! ####
			query_mat = query_mat[,which(colnames(query_mat) %in% query_seurat_cluster == T)]
		}
		########
		if(dim(query_mat)[2] < 2){
			print('Warning:: the query matrix must have at least 2 clusters !!!! Skip !!!!')
    		next
		}
		###### then Next we filter the ref datasets if the mode == 'hard' #######
		if(mode == 'hard'){
			### trim ref ####
			k = which(ref_seurat$celltype %in% overlap_ct == T)
			ref_cl = ref_seurat[,k]
			##### then we update ref_mat and ref label !!!! #######
			ref_mat = ref_cl[['RNA']]@data
			ref_label = unname(ref_cl$celltype)
		}
		######################################################
		####### prepare the datasets in the folder ###########
		####### scmap-cluster ################################
		if(method=='scmap-cluster'){
			library(scmap)
			library(SingleCellExperiment)
			#### train and test should be normalized counts #######
			query_mat_input = exp(query_mat)-1
			ref_mat_input = exp(ref_mat)-1
			ref_label_input = ref_label
			res = CellAnn_scmapcluster(ref_mat_input,query_mat_input,ref_label_input,threshold = 0.5,time = F)
			#####
			res_table = data.frame(cluster=colnames(query_mat_input),result=res)
			res_table = list(res_table)
			names(res_table) = TAG
			out_table_list = c(out_table_list,res_table)
			ref_label_list = c(ref_label_list,list(ref_label))
		}
		if(method=='chetah'){
			library(scmap)
			library(SingleCellExperiment)
			#### train and test should be normalized counts #######
			query_mat_input = exp(query_mat)-1
			ref_mat_input = exp(ref_mat)-1
			ref_label_input = ref_label
			res = CellAnn_chetah(ref_mat_input,query_mat_input,ref_label_input,time = F)
			#####
			res_table = data.frame(cluster=colnames(query_mat_input),result=res)
			res_table = list(res_table)
			names(res_table) = TAG
			out_table_list = c(out_table_list,res_table)
			ref_label_list = c(ref_label_list,list(ref_label))
		}
		if(method=='seurat'){
			library(Seurat)
			#### train and test should be normalized counts #######
			query_mat_input = exp(query_mat)-1
			ref_mat_input = exp(ref_mat)-1
			ref_label_input = ref_label
			res = CellAnn_seurat(ref_mat_input,query_mat_input,ref_label_input,time = F)
			#####
			res_table = data.frame(cluster=colnames(query_mat_input),result=res)
			res_table = list(res_table)
			names(res_table) = TAG
			out_table_list = c(out_table_list,res_table)
			ref_label_list = c(ref_label_list,list(ref_label))
		}
		if(method=='scpred'){
			library(Seurat)
			#### train and test should be normalized counts #######
			query_mat_input = exp(query_mat)-1
			ref_mat_input = exp(ref_mat)-1
			ref_label_input = ref_label
			res = CellAnn_scpred(ref_mat_input,query_mat_input,ref_label_input,time = F)
			#####
			res_table = data.frame(cluster=colnames(query_mat_input),result=res)
			res_table = list(res_table)
			names(res_table) = TAG
			out_table_list = c(out_table_list,res_table)
			ref_label_list = c(ref_label_list,list(ref_label))
		}
		if(method == 'scClassify'){
			library(scClassify)
			#### log-transformed (size-factor normalized) matrices as query datasets #####
			query_mat_input = query_mat
			ref_mat_input = ref_mat
			ref_label_input = ref_label
			####
			res = CellAnn_scClassify(query_mat_input,ref_mat_input,ref_label_input,time = F,prob_threshold=0.5)
			####
			res_table = data.frame(cluster=colnames(query_mat_input),result=res)
			res_table = list(res_table)
			names(res_table) = TAG
			out_table_list = c(out_table_list,res_table)
			ref_label_list = c(ref_label_list,list(ref_label))
		}

	}
	####### Then Next tools !!!! ##########################
	#######
	out_table_list_v = Visualize_res(out_table_list,ref_label_list)
	#######
	return(out_table_list_v)
}

############### Then we first run cellAnn for scmap-cluster tools !!! ##############
library(ggplot2)
folder = "C:/Users/plyu3/Desktop/CellAnn_methods_test/Tabula_Muris_mouse_data_prepare"
scmap_cluster_res0.5_Mouse_Self = Main_compare_process_OtherTools(Self_compare_list,method='scmap-cluster',folder=folder,mode='easy')
scmap_cluster_res0.5_Mouse_Self_p = Visualize_res_plot(scmap_cluster_res0.5_Mouse_Self,'mouse')
ggplot(scmap_cluster_res0.5_Mouse_Self_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("scmap_cluster_res0.5_Mouse_Self.png",width=8,height=4)

############### Then we first run cellAnn for tools !!! ##############
library(ggplot2)
folder = "C:/Users/plyu3/Desktop/CellAnn_methods_test/Tabula_Muris_mouse_data_prepare"
chetah_Mouse_Self = Main_compare_process_OtherTools(Self_compare_list,method='chetah',folder=folder,mode='easy')
chetah_Mouse_Self_p = Visualize_res_plot(chetah_Mouse_Self,'mouse')
ggplot(chetah_Mouse_Self_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("chetah_Mouse_Self.png",width=8,height=4)

################ OK then we use next tool!!! ############################################
library(ggplot2)
folder = "C:/Users/plyu3/Desktop/CellAnn_methods_test/Tabula_Muris_mouse_data_prepare"
seurat_Mouse_Self = Main_compare_process_OtherTools(Self_compare_list,method='seurat',folder=folder,mode='easy')
seurat_Mouse_Self_p = Visualize_res_plot(seurat_Mouse_Self,'mouse')
ggplot(seurat_Mouse_Self_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("seurat_Mouse_Self.png",width=8,height=4)


################ OK then we use next tool!!! ############################################
library(ggplot2)
folder = "C:/Users/plyu3/Desktop/CellAnn_methods_test/Tabula_Muris_mouse_data_prepare"
scpred_Mouse_Self = Main_compare_process_OtherTools(Self_compare_list,method='scpred',folder=folder,mode='easy')
scpred_Mouse_Self_p = Visualize_res_plot(scpred_Mouse_Self,'mouse')
ggplot(scpred_Mouse_Self_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')

ggsave("scpred_Mouse_Self.png",width=8,height=4)

################ Next scClassify ##########################################################
#####
##### n = 14 ### we need test !!! #####

library(ggplot2)
folder = "C:/Users/plyu3/Desktop/CellAnn_methods_test/Tabula_Muris_mouse_data_prepare"
scClassify_Mouse_Self = Main_compare_process_OtherTools(Self_compare_list,method='scClassify',folder=folder,mode='easy')
scClassify_Mouse_Self_p = Visualize_res_plot(scClassify_Mouse_Self,'mouse')
ggplot(scClassify_Mouse_Self_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')

ggsave("scClassify_Mouse_Self.png",width=8,height=4)


################## OK then next !!!! ########################################################
################## Next we recalculate the marker genes for each cell type ##################
##################

###########
"yellow"
"yellow"
"yellow"
"yellow"
"yellow"
"yellow"
"yellow"
###########
########### Next we test the batch effect compares for Mouse datasets !!!! #############################
###########
files = list.files()
m1 = grep('_DGE',files)
m2 = grep('_auth',files)
m3 = grep('_Clu',files)
m4 = grep('_clu',files)
m5 = grep('_Ground',files)
m6 = grep('_Sub',files)

files = files[-c(m1,m2,m3,m4,m5,m6)]

m7 = grep('_Ref',files)
m8 = grep('_test',files)

files = files[-c(m7,m8)]

files = files[-c(1:17)]

save(files,file='files')
#### with batch effect ####
mouse_df_batch1 <- data.frame(query=files[c(1,2,2,3:12)],ref=files[c(14,13,15,16:25)])
mouse_df_batch2 <- data.frame(query=files[c(14,13,15,16:25)],ref=files[c(1,2,2,3:12)])
mouse_df_batch = rbind(mouse_df_batch1,mouse_df_batch2)

save(mouse_df_batch,file='mouse_df_batch')


############ Then we start easy mode !!!!! ############
############
scmap_cluster_res0.5_Mouse_batch = Main_compare_process_OtherTools(mouse_df_batch,method='scmap-cluster',folder=folder,mode='easy')
scmap_cluster_res0.5_Mouse_batch_p = Visualize_res_plot(scmap_cluster_res0.5_Mouse_batch,'mouse')
ggplot(scmap_cluster_res0.5_Mouse_batch_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("scmap_cluster_res0.5_Mouse_Batch_Easy.png",width=8,height=4)

##### we are here!!! red red red #######################
scmap_cluster_res0.5_Mouse_batch = Main_compare_process_OtherTools(mouse_df_batch,method='scmap-cluster',folder=folder,mode='hard')
scmap_cluster_res0.5_Mouse_batch_p = Visualize_res_plot(scmap_cluster_res0.5_Mouse_batch,'mouse')
ggplot(scmap_cluster_res0.5_Mouse_batch_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("scmap_cluster_res0.5_Mouse_Batch_Hard.png",width=8,height=4)


############
############ OK!!! we will check the easy mode !!!! ########
############ Then we start easy mode!!!! ###################
############

chetah_Mouse_batch = Main_compare_process_OtherTools(mouse_df_batch,method='chetah',folder=folder,mode='easy')
chetah_Mouse_batch_p = Visualize_res_plot(chetah_Mouse_batch,'mouse')
ggplot(chetah_Mouse_batch_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("chetah_Mouse_Batch_Easy.png",width=8,height=4)

chetah_Mouse_batch = Main_compare_process_OtherTools(mouse_df_batch,method='chetah',folder=folder,mode='hard')
chetah_Mouse_batch_p = Visualize_res_plot(chetah_Mouse_batch,'mouse')
ggplot(chetah_Mouse_batch_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("chetah_Mouse_Batch_Hard.png",width=8,height=4)


########### Next !!!! ###############################
seurat_Mouse_batch = Main_compare_process_OtherTools(mouse_df_batch,method='seurat',folder=folder,mode='easy')
seurat_Mouse_batch_p = Visualize_res_plot(seurat_Mouse_batch,'mouse')
ggplot(seurat_Mouse_batch_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("seurat_Mouse_Batch_Easy.png",width=8,height=4)

seurat_Mouse_batch = Main_compare_process_OtherTools(mouse_df_batch,method='seurat',folder=folder,mode='hard')
seurat_Mouse_batch_p = Visualize_res_plot(seurat_Mouse_batch,'mouse')
ggplot(seurat_Mouse_batch_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("seurat_Mouse_Batch_Hard.png",width=8,height=4)


########### Next !!!! ###############################
scpred_Mouse_batch = Main_compare_process_OtherTools(mouse_df_batch,method='scpred',folder=folder,mode='easy')
scpred_Mouse_batch_p = Visualize_res_plot(scpred_Mouse_batch,'mouse')
ggplot(scpred_Mouse_batch_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("scpred_Mouse_Batch_Easy.png",width=8,height=4)

scpred_Mouse_batch = Main_compare_process_OtherTools(mouse_df_batch,method='scpred',folder=folder,mode='hard')
scpred_Mouse_batch_p = Visualize_res_plot(scpred_Mouse_batch,'mouse')
ggplot(scpred_Mouse_batch_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("scpred_Mouse_Batch_Hard.png",width=8,height=4)



########### Next !!!! ###############################
scClassify_Mouse_batch = Main_compare_process_OtherTools(mouse_df_batch,method='scClassify',folder=folder,mode='easy')
scClassify_Mouse_batch_p = Visualize_res_plot(scClassify_Mouse_batch,'mouse')
ggplot(scClassify_Mouse_batch_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("scClassify_Mouse_Batch_Easy.png",width=8,height=4)

scClassify_Mouse_batch = Main_compare_process_OtherTools(mouse_df_batch,method='scClassify',folder=folder,mode='hard')
scClassify_Mouse_batch_p = Visualize_res_plot(scClassify_Mouse_batch,'mouse')
ggplot(scClassify_Mouse_batch_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("scClassify_Mouse_Batch_Hard.png",width=8,height=4)


############### OK !!!! ########
############### let us check the results !!!!!! #####################################

mouse_hard_batch1 <- data.frame(query=files[c(25)],ref=files[c(23,22,21,20,18,17)])
mouse_hard_batch2 <- data.frame(query=files[c(14)],ref=files[c(23,22,21,20,18,17)])
mouse_hard_batch3 <- data.frame(query=files[c(12)],ref=files[c(23,22,21,20,18,17)])
mouse_hard_batch4 <- data.frame(query=files[c(1)],ref=files[c(23,22,21,20,18,17)])
mouse_hard_batch = rbind(mouse_hard_batch1,mouse_hard_batch2,mouse_hard_batch3,mouse_hard_batch4)
save(mouse_hard_batch,file='mouse_hard_batch')



############### OK !!! #######################################################
###################################### let us check between datasets #########

scmap_cluster_res0.5_Mouse_batch = Main_compare_process_OtherTools(mouse_hard_batch,method='scmap-cluster',folder=folder,mode='NO')
scmap_cluster_res0.5_Mouse_batch_p = Visualize_res_plot(scmap_cluster_res0.5_Mouse_batch,'mouse')
ggplot(scmap_cluster_res0.5_Mouse_batch_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("scmap_cluster_res0.5_Mouse_Batch_VeryHard.png",width=8,height=4)

chetah_Mouse_batch = Main_compare_process_OtherTools(mouse_hard_batch,method='chetah',folder=folder,mode='NO')
chetah_Mouse_batch_p = Visualize_res_plot(chetah_Mouse_batch,'mouse')
ggplot(chetah_Mouse_batch_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("chetah_Mouse_Batch_VeryHard.png",width=8,height=4)

seurat_Mouse_batch = Main_compare_process_OtherTools(mouse_hard_batch,method='seurat',folder=folder,mode='NO')
seurat_Mouse_batch_p = Visualize_res_plot(seurat_Mouse_batch,'mouse')
ggplot(seurat_Mouse_batch_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("seurat_Mouse_Batch_VeryHard.png",width=8,height=4)

scpred_Mouse_batch = Main_compare_process_OtherTools(mouse_hard_batch,method='scpred',folder=folder,mode='NO')
scpred_Mouse_batch_p = Visualize_res_plot(scpred_Mouse_batch,'mouse')
ggplot(scpred_Mouse_batch_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("scpred_Mouse_Batch_VeryHard.png",width=8,height=4)

scClassify_Mouse_batch = Main_compare_process_OtherTools(mouse_hard_batch,method='scClassify',folder=folder,mode='NO')
scClassify_Mouse_batch_p = Visualize_res_plot(scClassify_Mouse_batch,'mouse')
ggplot(scClassify_Mouse_batch_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("scClassify_Mouse_Batch_VeryHard.png",width=8,height=4)











###### yellow ####
###### yellow ####
###### yellow ####
###### yellow ####

######
###### OK! let us prepare the Human datasets !!!! ##########
######



setwd('C:/Users/plyu3/Desktop/CellAnn_methods_test/Human_Pancreas_test2/')

files = list.files()
files = files[grep('pre3$',files)]


###### let us see these cell names !!!! #####################
####

_test_input_ClusterAvg

_test_input_GroundTruthLabel

_Ref_CTsubDEG_ORDER150

_Ref_CTDEG


#####
Ref 是 根据 UMAP 对每一个 Cell type 分类 #########
#####
Query 是 根据 cluster 分每一个 cell type ##########
#####


RNA_process_UMAP_Cluster <- function(x,res){
	#####
	DefaultAssay(x) = 'RNA'
	#####
	x <- NormalizeData(x)
    x <- FindVariableFeatures(x,selection.method ='vst',nfeatures = 5000)
    x <- ScaleData(x,  verbose = FALSE)
    x <- RunPCA(x, verbose = FALSE,npcs=100)
    x <- RunUMAP(x, reduction = "pca", dims = 1:100)
	x <- FindNeighbors(x, reduction = "pca", dims = 1:100)
	x <- FindClusters(x, resolution = res)
	#####
	return(x)
}

RNA_process_Cluster_to_CT <- function(x,tag){
	x$seurat_clusters = as.numeric(x$seurat_clusters)
	#### rm NA #####
	k = which(is.na(x$celltype) == T)
	if(length(k)>0){x = x[,-k]}
	####
	#### see the celltype and seurat_clusters table ####
	m = match(c('celltype','seurat_clusters'),colnames(x@meta.data))
	####
	tab = x@meta.data[,m]
	allcelltypes = levels(as.factor(tab$celltype))
	####
	align_table = list()
	for(i in 1:length(allcelltypes)){
		print(allcelltypes[i])
		tab_tmp = tab[which(tab$celltype == allcelltypes[i]),]
		tab_tmp_table = table(as.numeric(tab_tmp$seurat_clusters))
		print(tab_tmp_table)
		print(paste('total',sum(tab_tmp_table)))
		#######
		tab_tmp_res = tab_tmp
		rownames(tab_tmp_res) = NULL
		#######
		tab_tmp_res$index = paste(tab_tmp_res$celltype,as.numeric(tab_tmp_res$seurat_clusters),sep='@')
		tab_tmp_res = tab_tmp_res[!duplicated(tab_tmp_res$index),]
		m = match(tab_tmp_res$seurat_clusters,names(tab_tmp_table))
		tab_tmp_res$cellcount = as.numeric(tab_tmp_table)[m]
		#######
		####### added new clusters #######
		tab_tmp_res$new_clusters = paste(tab_tmp_res$seurat_clusters,i,sep='_')
		######## rm the cluster with too few cells #######
		k = which(tab_tmp_res$cellcount < 15)
		if(length(k) > 0 & sum(tab_tmp_table) > 50 & length(k) != length(tab_tmp_res$cellcount)){
			tab_tmp_res = tab_tmp_res[-k,]
		}
		if(length(k) > 0 & sum(tab_tmp_table) > 50 & length(k) == length(tab_tmp_res$cellcount)){
			tab_tmp_res$new_clusters = tab_tmp_res$new_clusters[1]
		}
		if(sum(tab_tmp_table) < 50){
			tab_tmp_res$new_clusters = tab_tmp_res$new_clusters[1]
		}
		########
		align_table = c(align_table,list(tab_tmp_res))
	}
	align_table = do.call('rbind',align_table)
	####
	#### Then replace the seurat clusters ######
	x$seurat_clusters_new = 'removed'
	for(i in 1:dim(align_table)[1]){
		tmp_celltype = align_table$celltype[i]
		tmp_seurat_clusters = align_table$seurat_clusters[i]
		k = which(x$celltype == tmp_celltype & x$seurat_clusters == tmp_seurat_clusters)
		x$seurat_clusters_new[k] = align_table$new_clusters[i]
	}
	x = subset(x,subset=seurat_clusters_new != 'removed')
	x$seurat_clusters_new = paste('C',x$seurat_clusters_new,sep='')
	####
	####
	x$seurat_clusters = x$seurat_clusters_new
	####
	#######
	png_file = paste(tag,'_test_input_cluster_test','.png',sep='')
	library(ggplot2)
	png(png_file,height=4000,width=5000,res=72*12)
	print(DimPlot(x, reduction = "umap",group.by='seurat_clusters_new',label = FALSE, label.size = 2.5, repel = TRUE))
	dev.off()
	png_file = paste(tag,'_test_input_celltype_author_test','.png',sep='')
	library(ggplot2)
	png(png_file,height=4000,width=5000,res=72*12)
	print(DimPlot(x, reduction = "umap",group.by='celltype',label = FALSE, label.size = 2.5, repel = TRUE))
	dev.off()
	#######
	return(x)
}



CellAnn_Avg_Mat <- function(data_mat,data_cluster,log='log',scale_factor=10000){
	######
	tag_cluster = unname(data_cluster)
	tag_cluster_level = levels(as.factor(tag_cluster))
	###### normalized back datasets ######
	if(log == 'log'){
		data_mat_exp = exp(data_mat)
		data_mat_exp = data_mat_exp-1
	}
	if(log == 'log2'){
		data_mat_exp = 2^(data_mat)
		data_mat_exp = data_mat_exp-1
	}
	print(paste('Sums:',head(colSums(data_mat_exp[,c(1:2)]))))
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
	merge_mat = round(log(merge_mat+1),5)
	return(merge_mat)
}


Avg_mat_to_df <- function(mat){
	df_1 = data.frame(GENE=rownames(mat))
	##
	df_2 = data.frame(mat)
	##
	out = cbind(df_1,df_2)
	##
	return(out)
}


df_to_mat <- function(df){
	rowN = df$GENE
	##
	mat = as.matrix(df[,-1])
	##
	rownames(mat) = rowN
	##
	return(mat)
}


#####
head(Seurat_Obj@meta.data)

library(Seurat)
index = files

########
######## Muraro and Xin need to replace!!!! ##########
######## 

setwd("/zp1/data/plyu3/M_G_expression_datasets")
Human_gtf = rtracklayer::import("gencode.v41.basic.annotation.gtf.gz")
Human_gtf = data.frame(Human_gtf)
Human_gtf = Human_gtf[,c(10,12)]
Human_gtf_index = paste(Human_gtf[,1],Human_gtf[,2],sep='@')
Human_gtf = Human_gtf[!duplicated(Human_gtf_index),]
Human_gtf[,1] = sapply(strsplit(Human_gtf[,1],split='.',fixed=T),function(x) x[[1]])
save(Human_gtf,file='Human_gtf')

######## replace genes #######
########

i = 3
tag = index[i]
Seurat_Obj = loadRData(tag)
Counts = Seurat_Obj[['RNA']]@counts
load('Human_gtf')

mat = Counts

Replaces_rownames_for_matrix <- function(mat,Human_gtf){
	k = which(rownames(mat) %in% Human_gtf[,1] == T)
	mat = mat[k,]
	m = match(rownames(mat),Human_gtf[,1])
	#####
	rownames(mat) = Human_gtf[m,2]
	#####
	return(mat)
}
head(colnames(Counts))
head(colnames(Counts_new))
Counts_new = Replaces_rownames_for_matrix(Counts,Human_gtf)
New_Seurat = CreateSeuratObject(Counts_new)
New_Seurat@meta.data = Seurat_Obj@meta.data

save(New_Seurat,file=tag)

i=6
tag = index[i]
Seurat_Obj = loadRData(tag)
Counts = Seurat_Obj[['RNA']]@counts
load('Human_gtf')
Counts_new = Replaces_rownames_for_matrix(Counts,Human_gtf)
New_Seurat = CreateSeuratObject(Counts_new)
New_Seurat@meta.data = Seurat_Obj@meta.data
save(New_Seurat,file=tag)



for(i in 3){
	#####
	tag = index[i]
	print(tag)
	Seurat_Obj = loadRData(tag)
	#####
	print(head(rownames(Seurat_Obj)))
	#####
	Seurat_Obj = RNA_process_UMAP_Cluster(Seurat_Obj,2)
	#####
	##### first prepare query #####
	#####
	Seurat_Obj <- RNA_process_Cluster_to_CT(Seurat_Obj,tag=tag)
	#####
	data = Seurat_Obj
	clusters = data$seurat_clusters
	celltypes = data$celltype
	####
	dat = data.frame(clusters,celltypes)
	index33 = paste(dat[,1],dat[,2])
	dat = dat[!duplicated(index33),]
	colnames(dat) = c('cluster','ground-truth')
	######
	new_File = paste(tag,'_test_input_GroundTruthLabel','.txt',sep='')
	write.table(dat,file=new_File,sep='\t',quote=F,row.names=F)
	######
	###### Next Query input ######
	######
	data_mat = Seurat_Obj[['RNA']]@data
	####
	data_cluster = Seurat_Obj$seurat_clusters
	####
	new_AvgFile = paste(tag,'_test_input_ClusterAvg.txt','.txt',sep='')
	####
	avg_mat = CellAnn_Avg_Mat(data_mat,data_cluster)
	####
	avg_df = Avg_mat_to_df(avg_mat)
	####
	write.table(avg_df,file=new_AvgFile,sep='\t',quote=F,row.names=F)
	####
	#### Next prepare ref ######
	####
	subset_matrix_res = subset_each_ct(data,resolution=0.3)
	subset_matrix = subset_matrix_res[[1]]
	subset_cell = subset_matrix_res[[2]]
	####
	new_File = paste(tag,'_Ref_SubCTAvg',sep='')
	####
	saveRDS(subset_matrix,file=new_File)
	#### Then calculate DEGs for sub ######
	m = match(colnames(data),subset_cell$cells)
	data$cell_id = subset_cell$cells[m]
	data$celltype_sub = subset_cell$new_sub_cluster[m]
	####
	print(tag)
	save(data,file=tag)
}


##########
##########

##########
##########

j = 7

subset_each_ct <- function(data,resolution=0.3){
	#####
	matrix_list = list()
	##### ct ######
	ct = levels(as.factor(data$celltype))
	#####
	cells_list = list()
	#####
	for(j in 1:length(ct)){
		print(j)
		print(ct[j])
		#######
		k = which(data$celltype == ct[j])
		#######
		if(length(k) > 1){
			print('large')
			#######
			sub_seurat = subset(data,subset = celltype == ct[j])
			sub_seurat = FindNeighbors(sub_seurat,reduction = "umap",dims = 1:2)
			sub_seurat = FindClusters(sub_seurat,resolution=resolution)
			sub_seurat_mat = sub_seurat[['RNA']]@data
			data_cluster = sub_seurat$seurat_clusters
			data_cluster = paste(ct[j],data_cluster,sep='@sub')
			sub_seurat_avg = CellAnn_Avg_Mat(data_mat=sub_seurat_mat,data_cluster)
			matrix_list = c(matrix_list,list(sub_seurat_avg))
			#######
			cells_table = data.frame(cells=colnames(sub_seurat_mat),new_sub_cluster = data_cluster)
			cells_list = c(cells_list,list(cells_table))
		}
	}
	matrix_list = do.call('cbind',matrix_list)
	cells_list = do.call('rbind',cells_list)
	####
	return(list(matrix_list,cells_list))
}


############
####### OK! then we perform self alignment !!! ###########
############

Main_compare_process_OtherTools2 <- function(compare_df,method='scmap-cluster',folder=folder,mode='easy'){
	library(Seurat)
	##########
	out_table_list = list()
	ref_label_list = list()
	##########
	for(i in 1:dim(compare_df)[1]){
		print(paste("num",i))
		print(paste('query:',compare_df[i,1],'  ','ref:',compare_df[i,2],sep=''))
		TAG = paste('query:',compare_df[i,1],'-->','ref:',compare_df[i,2],sep='')
		#######
		query = compare_df[i,1]
		ref = compare_df[i,2]
		query_index = paste0(query)
		ref_index = paste0(ref)
		####### load ref and query seurat object !!!! #######
		query_seurat = loadRData(query_index)
		ref_seurat = loadRData(ref_index)
		#######	OK!
		####### we next find how is the overlap between ref and query datasets ##########
		query_seurat_ct = levels(as.factor(query_seurat$celltype))
		ref_seurat_ct = levels(as.factor(ref_seurat$celltype))
		####### then we see how many overlap cell types between query and ref ###########
		overlap_ct = query_seurat_ct[which(query_seurat_ct %in% ref_seurat_ct == T)]
		#######
		if(length(overlap_ct) == 0){
			print('Warning: no overlaps between ref and query datasets !')
		}
		####### then we load query input matrix !!!! #####
		query_avg = paste(query,'_test_input_ClusterAvg.txt.txt',sep='')
		query_mat = read.table(query_avg,header=T)
		query_mat = df_to_mat(query_mat)
		####### ref input mat :: ref_mat #######
		####### ref input cell types:: ref_label
		ref_mat = ref_seurat[['RNA']]@data
		ref_label = unname(ref_seurat$celltype)
		####### then we will filter query datasets if the mode == 'easy' !!!! ####
		if(mode == 'easy'){
			#### we will filter the query matrix !!! #######
			#### we need to load the GroundTruth to get the cell types of the query matrix ########
			query_cluster_truth_index = paste(query,'_test_input_GroundTruthLabel.txt',sep='')
			query_cluster_truth = read.table(query_cluster_truth_index,sep='\t',header=T)
			#### clusters in query must be the overlapped ct #####
			query_seurat_cluster = query_cluster_truth$cluster[which(query_cluster_truth$ground.truth %in% overlap_ct == T)]
			#### Then we update the query matrix !!!! ####
			query_mat = query_mat[,which(colnames(query_mat) %in% query_seurat_cluster == T)]
		}
		########
		if(dim(query_mat)[2] < 2){
			print('Warning:: the query matrix must have at least 2 clusters !!!! Skip !!!!')
    		next
		}
		###### then Next we filter the ref datasets if the mode == 'hard' #######
		if(mode == 'hard'){
			### trim ref ####
			k = which(ref_seurat$celltype %in% overlap_ct == T)
			ref_cl = ref_seurat[,k]
			##### then we update ref_mat and ref label !!!! #######
			ref_mat = ref_cl[['RNA']]@data
			ref_label = unname(ref_cl$celltype)
		}
		######################################################
		####### prepare the datasets in the folder ###########
		####### scmap-cluster ################################
		if(method=='scmap-cluster'){
			library(scmap)
			library(SingleCellExperiment)
			#### train and test should be normalized counts #######
			query_mat_input = exp(query_mat)-1
			ref_mat_input = exp(ref_mat)-1
			ref_label_input = ref_label
			res = CellAnn_scmapcluster(ref_mat_input,query_mat_input,ref_label_input,threshold = 0.5,time = F)
			#####
			res_table = data.frame(cluster=colnames(query_mat_input),result=res)
			res_table = list(res_table)
			names(res_table) = TAG
			out_table_list = c(out_table_list,res_table)
			ref_label_list = c(ref_label_list,list(ref_label))
		}
		if(method=='chetah'){
			library(scmap)
			library(SingleCellExperiment)
			#### train and test should be normalized counts #######
			query_mat_input = exp(query_mat)-1
			ref_mat_input = exp(ref_mat)-1
			ref_label_input = ref_label
			res = CellAnn_chetah(ref_mat_input,query_mat_input,ref_label_input,time = F)
			#####
			res_table = data.frame(cluster=colnames(query_mat_input),result=res)
			res_table = list(res_table)
			names(res_table) = TAG
			out_table_list = c(out_table_list,res_table)
			ref_label_list = c(ref_label_list,list(ref_label))
		}
		if(method=='seurat'){
			library(Seurat)
			#### train and test should be normalized counts #######
			query_mat_input = exp(query_mat)-1
			ref_mat_input = exp(ref_mat)-1
			ref_label_input = ref_label
			res = CellAnn_seurat(ref_mat_input,query_mat_input,ref_label_input,time = F)
			#####
			res_table = data.frame(cluster=colnames(query_mat_input),result=res)
			res_table = list(res_table)
			names(res_table) = TAG
			out_table_list = c(out_table_list,res_table)
			ref_label_list = c(ref_label_list,list(ref_label))
		}
		if(method=='scpred'){
			library(Seurat)
			#### train and test should be normalized counts #######
			query_mat_input = exp(query_mat)-1
			ref_mat_input = exp(ref_mat)-1
			ref_label_input = ref_label
			res = CellAnn_scpred(ref_mat_input,query_mat_input,ref_label_input,time = F)
			#####
			res_table = data.frame(cluster=colnames(query_mat_input),result=res)
			res_table = list(res_table)
			names(res_table) = TAG
			out_table_list = c(out_table_list,res_table)
			ref_label_list = c(ref_label_list,list(ref_label))
		}
		if(method == 'scClassify'){
			library(scClassify)
			#### log-transformed (size-factor normalized) matrices as query datasets #####
			query_mat_input = query_mat
			ref_mat_input = ref_mat
			ref_label_input = ref_label
			####
			res = CellAnn_scClassify(query_mat_input,ref_mat_input,ref_label_input,time = F,prob_threshold=0.5)
			####
			res_table = data.frame(cluster=colnames(query_mat_input),result=res)
			res_table = list(res_table)
			names(res_table) = TAG
			out_table_list = c(out_table_list,res_table)
			ref_label_list = c(ref_label_list,list(ref_label))
		}

	}
	####### Then Next tools !!!! ##########################
	#######
	out_table_list_v = Visualize_res(out_table_list,ref_label_list)
	#######
	return(out_table_list_v)
}


setwd('C:/Users/plyu3/Desktop/CellAnn_methods_test/Human_Pancreas_test2/')

files = list.files()
files = files[grep('pre3$',files)]

Human_self_df = data.frame(query=files,ref=files)

library(ggplot2)
folder = "C:/Users/plyu3/Desktop/CellAnn_methods_test/Human_Pancreas_test2/"

compare_df = Human_self_df

scmap_cluster_res0.5_H_Self = Main_compare_process_OtherTools2(Human_self_df[c(1:4,6),],method='scmap-cluster',folder=folder,mode='easy')
scmap_cluster_res0.5_H_Self_p = Visualize_res_plot(scmap_cluster_res0.5_H_Self,'H')
ggplot(scmap_cluster_res0.5_H_Self_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("scmap_cluster_res0.5_H_Self.png",width=8,height=4)

############### Then we first run cellAnn for tools !!! ##############
library(ggplot2)
chetah_Mouse_Self = Main_compare_process_OtherTools2(Human_self_df[c(1:4,6),],method='chetah',folder=folder,mode='easy')
chetah_Mouse_Self_p = Visualize_res_plot(chetah_Mouse_Self,'H')
ggplot(chetah_Mouse_Self_p,aes(x=sample2,y=counts,fill=class)) + geom_bar(position="stack", stat="identity") + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('#006994','lightblue','lightgreen','grey','red','pink')) + ylab('Number of clusters') + xlab('')
ggsave("chetah_H_Self.png",width=8,height=4)




