####
#### this script test the memory usage in R #####
####

ssh plyu3@10.112.40.197
plyu3 njd$rft1


conda activate CellAnn_test
R

library(Seurat)
setwd("/zp1/data/plyu3/MCL_V3")
####
####
library("bench")
####
#### we will first get the compare_df !!! #######
####
library(Seurat)
setwd("/zp1/data/plyu3/MCL_V3")
####
####
#### 看来只能重新换datasets了 #####
#### using PBMC is the best solution !!! #######
####
#### line 2780 for prepare ###########
#### PBMC 没有那么多cell type ##########

cd /zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare

conda activate CellAnn_test
R

library(Seurat)

setwd("/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")

load('PBMC_Drop')
load('PBMC_inDrops')
load('PBMC_SeqWell')

load("PBMC_Chromium_V2")
load("PBMC_Chromium_V2A")
load("PBMC_Chromium_V2B")

load("PBMC_Chromium_V3")

PBMC_combine <- merge(PBMC_Drop,y = c(PBMC_inDrops,PBMC_SeqWell,PBMC_Chromium_V2,PBMC_Chromium_V2A,PBMC_Chromium_V2B,PBMC_Chromium_V3))

#### OK!!! ####

save(PBMC_combine,file="PBMC_combine")

#### OK!!! ####

size = c(100,500,1000,5000,10000,15000,20000,25000)


RNA_process_UMAP_Cluster <- function(x,res){
	#####
	DefaultAssay(x) = 'RNA'
	#####
	x <- NormalizeData(x)
    x <- FindVariableFeatures(x,selection.method ='vst',nfeatures = 2000)
    x <- ScaleData(x,  verbose = FALSE)
    x <- RunPCA(x, verbose = FALSE,npcs=50)
    x <- RunUMAP(x, reduction = "pca", dims = 1:50)
	x <- FindNeighbors(x, reduction = "pca", dims = 1:50)
	x <- FindClusters(x, resolution = res)
	#####
	return(x)
}

data = PBMC_combine

total_length = dim(data)[2]

for(i in 1:length(size)){
	#######
	tmp_size = size[i]
	#######
	tmp_random_list = sample(1:total_length, tmp_size)
	#######
	tmp_seurat = data[,tmp_random_list]
	####### then we process the tmp tmp_seurat ########
	tmp_seurat = RNA_process_UMAP_Cluster(tmp_seurat,res=5)
	print(table(tmp_seurat$seurat_clusters))
	#######
	tmp_file = paste0('PBMC_Ref_speed_test','_',tmp_size)
	#######
	save(tmp_seurat,file=tmp_file)
}

setwd("/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")

######### let us prepare the reference files for cellAnn ###################
#########


files = list.files()

files = files[grep("0$",files)]


CellAnn_Avg_Mat2 <- function(data_mat,data_cluster,log='log',scale_factor=10000){
	######
	tag_cluster = unname(data_cluster)
	tag_cluster_level = levels(as.factor(tag_cluster))
	###### normalized back datasets ######
	data_mat_exp = data_mat
	###### data_mat_exp is 1e5 normalize #######
	merge_mat = c()
	for(i in 1:length(tag_cluster_level)){
		index = which(data_cluster %in% tag_cluster_level[i] == T)
		index_mat = data_mat_exp[,index]
		print(dim(index_mat))
		######
		index_sum = Matrix::rowSums(index_mat)
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

for(i in 1:length(files)){
	######## "seurat_clusters_new" ########
	tmp_seurat = loadRData(files[i])
	########
	tmp_file_name_avg = paste0(files[i],'_SubCluster_Avg')
	########
	tmp_seurat_mat = tmp_seurat[['RNA']]@counts
	data_cluster = tmp_seurat$seurat_clusters
	ref_Avg_Mat <- CellAnn_Avg_Mat2(data_mat=tmp_seurat_mat,data_cluster)
	########
	save(ref_Avg_Mat,file=tmp_file_name_avg)
	######## Next we need to call cell markers for each main Cts #############
	tmp_file_name_CT = paste0(files[i],'_CT_Marker')
	tmp_ct_marker = runDEGs_Ref_sub(tmp_seurat,'COSG','seurat_clusters',100)
	save(tmp_ct_marker,file=tmp_file_name_CT)
	########
	tmp_file_name_CT = paste0(files[i],'_SubCT_Marker')
	tmp_ct_marker = runDEGs_Ref_sub(tmp_seurat,'COSG','seurat_clusters',100)
	save(tmp_ct_marker,file=tmp_file_name_CT)
	######## Next we need to call subCell markers for each sub Cts ###########
}

###
### OK!!! Next we will prepare the query datasets ####
###
### prepare the "_Cluster_Avg" ###
###

### line 1789 ####

tmp_query = "PBMC_inDrops"
tmp_seurat = loadRData(tmp_query)
tmp_file_name_avg = paste0(tmp_query,'_Cluster_Avg')
tmp_seurat_mat = tmp_seurat[['RNA']]@data
data_cluster = tmp_seurat$seurat_clusters
query_Avg_Mat <- CellAnn_Avg_Mat(data_mat=tmp_seurat_mat,data_cluster)

save(query_Avg_Mat,file=tmp_file_name_avg)


celltypes = tmp_seurat$celltype
dat = data.frame(data_cluster,celltypes)
#### Then we calculate the cell numbers for each !!!! #######
dat_res = table(data_cluster)
####
index = paste(dat[,1],dat[,2])
dat = dat[!duplicated(index),]
colnames(dat) = c('cluster','ground-truth')
m = match(dat$cluster,names(dat_res))
dat$number = as.numeric(dat_res)[m]
rownames(dat) = NULL
print(dat)

new_File = paste(tmp_query,'_GroundTruth_Cluster',sep='')
save(dat,file=new_File)


####
####




####
####
####

bnch <- bench::mark(
  data.frame(x = runif(4000, 1, 1000), y=runif(4000, 1, 1000))
)
print(max(bnch$mem_alloc))
print(max(bnch$total_time))

####
####
#### 2884 ####
####
####



files = list.files()
Ref_new = paste0('PBMC_Ref_speed_test_',as.character(c(100,500,1000,5000,10000,15000,20000,25000)))
compare_df = data.frame(Query_new="PBMC_inDrops",Ref_new=Ref_new)
folder = "/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare"



loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}


Main_compare_process_OtherTools_Single_Cell_times_MEM <- function(compare_df,method='scmap-cluster',folder=folder){
	library(Seurat)
	##########
	times_list = list()
	mem_list = list()
	##########
	for(i in 1:dim(compare_df)[1]){
		#######
		TAG = paste('query:',compare_df$Query_new[i],'-->','ref:',compare_df$Ref_new[i],sep='')
		print(paste("NONONONONONONONONO",i))
		print(TAG)
		#######
		query = compare_df$Query_new[i]
		ref = compare_df$Ref_new[i]
		####### load ref and query seurat object !!!! #######
		setwd(folder)
		#######
		query_seurat = loadRData(query)
		ref_seurat = loadRData(ref)
		#######
		### gene1 = rownames(query_seurat)
		### gene2 = rownames(ref_seurat)
		#######	OK! ############
		query_mat = query_seurat[['RNA']]@data
		query_label = unname(query_seurat$seurat_clusters)
		#######
		ref_mat = ref_seurat[['RNA']]@data
		ref_label = unname(ref_seurat$celltype)
		####### then we will filter query datasets if the mode == 'easy' !!!! ####
		####### prepare the datasets in the folder #####
		####### scmap-cluster First we don't need to try scmap-cluster !!!!! ################################
		if(method=='scmap-cluster'){
			library(scmap)
			library(SingleCellExperiment)
			#### we should load the input average expression data as the input !!!!! #########
			query_avg = paste0(query,'_Cluster_Avg') ### red ####
			query_mat = loadRData(query_avg)
			#### we will use the threshold to 0.5 #####
			#### then we load the query labels !!!! ###
			####
			query_label_index = paste0(query,'_GroundTruth_Cluster')
			query_label_tab = loadRData(query_label_index)
			m1 = match(colnames(query_mat),query_label_tab$cluster)
			query_label_tab = query_label_tab[m1,]
			####
			####
			#### train and test should be normalized counts #######
			query_mat_input = exp(query_mat)-1
			ref_mat_input = exp(ref_mat)-1
			ref_label_input = ref_label
			t1<-Sys.time()
			res = bench::mark(CellAnn_scmapcluster(ref_mat_input,query_mat_input,ref_label_input,threshold = 0.5,time = F))
			t2<-Sys.time()
			#diff_time = difftime(t2,t1,units="secs")
			######
			res1 = list(res$mem_alloc)
			names(res1) = TAG
			mem_list = c(mem_list,res1)
			#########
			res2 = list(res$total_time)
			names(res2) = TAG
			times_list = c(times_list,res2)
			#########
		}
		if(method=='chetah'){
			library(CHETAH)
			library(SingleCellExperiment)
			#### train and test should be normalized counts #######
			query_mat_input = exp(query_mat)-1
			ref_mat_input = exp(ref_mat)-1
			ref_label_input = ref_label
			t1<-Sys.time()
			res = bench::mark(CellAnn_chetah(ref_mat_input,query_mat_input,ref_label_input,time = F))
			t2<-Sys.time()
			######
			res1 = list(res$mem_alloc)
			names(res1) = TAG
			mem_list = c(mem_list,res1)
			#########
			res2 = list(res$total_time)
			names(res2) = TAG
			times_list = c(times_list,res2)
			#########
		}
		if(method=='seurat'){
			library(Seurat)
			#### train and test should be normalized counts #######
			query_mat_input = exp(query_mat)-1
			ref_mat_input = exp(ref_mat)-1
			ref_label_input = ref_label
			t1<-Sys.time()
			res = bench::mark(CellAnn_seurat(ref_mat_input,query_mat_input,ref_label_input,time = F))
			t2<-Sys.time()
			######
			res1 = list(res$mem_alloc)
			names(res1) = TAG
			mem_list = c(mem_list,res1)
			#########
			res2 = list(res$total_time)
			names(res2) = TAG
			times_list = c(times_list,res2)
			#########
		}
		if(method=='scpred'){
			#### train and test should be normalized counts #######
			query_mat_input = exp(query_mat)-1
			ref_mat_input = exp(ref_mat)-1
			ref_label_input = ref_label
			ref_label_input = as.character(ref_label_input)
			ref_label_input = paste0('C',ref_label_input)
			#########
			print(table(ref_label_input))
			#########
			Overlap_gene = rownames(query_mat_input)[which(rownames(query_mat_input) %in% rownames(ref_mat_input) == T)]
			query_mat_input = query_mat_input[which(rownames(query_mat_input) %in% Overlap_gene == T),]
			ref_mat_input = ref_mat_input[which(rownames(ref_mat_input) %in% Overlap_gene == T),]
			print(dim(query_mat_input))
			print(dim(ref_mat_input))
			#########
			t1<-Sys.time()
			res = bench::mark(CellAnn_scpred(ref_mat_input,query_mat_input,ref_label_input,time = F))
			t2<-Sys.time()
			######
			res1 = list(res$mem_alloc)
			names(res1) = TAG
			mem_list = c(mem_list,res1)
			#########
			res2 = list(res$total_time)
			names(res2) = TAG
			times_list = c(times_list,res2)
			#########
		}
		if(method == 'scClassify'){
			library(scClassify)
			#### log-transformed (size-factor normalized) matrices as query datasets #####
			query_mat_input = query_mat
			ref_mat_input = ref_mat
			ref_label_input = ref_label
			####
			t1<-Sys.time()
			print(head(rownames(query_mat_input)))
			print(head(rownames(ref_mat_input)))
			k = which(rownames(query_mat_input) %in% rownames(ref_mat_input) == T)
			print(length(k))
			res = bench::mark(CellAnn_scClassify(query_mat_input,ref_mat_input,ref_label_input,time = F,prob_threshold=0.5))
			t2<-Sys.time()
			######
			res1 = list(res$mem_alloc)
			names(res1) = TAG
			mem_list = c(mem_list,res1)
			#########
			res2 = list(res$total_time)
			names(res2) = TAG
			times_list = c(times_list,res2)
			#########
		}
	}
	####### Then Next tools !!!! ##########################
	####### names(times_list) = TAG
	#######
	return(list(times_list,mem_list))
}


#### scmap !!! #################

scmap_cluster_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[1:8,],method='scmap-cluster',folder="/zp1/data/plyu3/MCL_V3")

save(scmap_cluster_MEM,file="scmap_cluster_MEM")

cp MEM_test_other_tools.R scmap_test_other_tools.R

conda activate CellAnn_test
R

nohup Rscript scmap_test_other_tools.R &

#### OK!!! Next tools !!!!! ####

chetah_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[1:7,],method='chetah',folder="/zp1/data/plyu3/MCL_V3")

save(chetah_MEM,file="chetah_MEM")

cp scmap_test_other_tools.R chetah_test_other_tools.R
nohup Rscript chetah_test_other_tools.R &



#### OK!!! Next tools !!!! #####

seurat_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df,method='seurat',folder="/zp1/data/plyu3/MCL_V3")

save(seurat_MEM,file="seurat_MEM")

cp scmap_test_other_tools.R seurat_test_other_tools.R

nohup Rscript seurat_test_other_tools.R &

Rscript seurat_test_other_tools.R

#### OK!!! Next tools !!!!! #####

scPred_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df,method='scpred',folder="/zp1/data/plyu3/MCL_V3")

save(scPred_MEM,file="scPred_MEM")

cd /zp1/data/plyu3/MCL_V3

cp MEM_test_other_tools.R scPred1_MEM_test_other_tools.R

vi scPred1_MEM_test_other_tools.R



cp MEM_test_other_tools.R scPred2_MEM_test_other_tools.R
cp MEM_test_other_tools.R scPred3_MEM_test_other_tools.R
cp MEM_test_other_tools.R scPred4_MEM_test_other_tools.R
cp MEM_test_other_tools.R scPred5_MEM_test_other_tools.R
cp MEM_test_other_tools.R scPred6_MEM_test_other_tools.R
cp MEM_test_other_tools.R scPred7_MEM_test_other_tools.R

#### 转移到新的服务器上面 #######

ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

#### 先安装scclassify #########
BiocManager::install("scClassify")

conda activate seurat4

conda install -c bioconda r-scpred







##### OK!!! Next tools !!!! #####

scClassify_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[1:7,],method='scClassify',folder="/zp1/data/plyu3/MCL_V3")

save(scClassify_MEM,file="scClassify_MEM")


cp MEM_test_other_tools.R scClassify1_MEM_test_other_tools.R

scClassify_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[1:3,],method='scClassify',folder="/zp1/data/plyu3/MCL_V3")

nohup Rscript scClassify1_MEM_test_other_tools.R &

cp scClassify1_MEM_test_other_tools.R scClassify4_MEM_test_other_tools.R

nohup Rscript scClassify4_MEM_test_other_tools.R &

cp scClassify1_MEM_test_other_tools.R scClassify5_MEM_test_other_tools.R

nohup Rscript scClassify5_MEM_test_other_tools.R &

cp scClassify1_MEM_test_other_tools.R scClassify6_MEM_test_other_tools.R

nohup Rscript scClassify6_MEM_test_other_tools.R &

cp scClassify1_MEM_test_other_tools.R scClassify7_MEM_test_other_tools.R

nohup Rscript scClassify7_MEM_test_other_tools.R &


#####
##### OK!!! we first need to check scClassify_MEM ######
#####




##### OK!!! Next is cellAnn !!! #######
##### OK!!! Next is cellAnn !!! #######
##### OK!!! Next is cellAnn !!! #######
#####


#####
#####
conda activate CellAnn_test
R
library(Seurat)
setwd("/zp1/data/plyu3/MCL_V3")
library("bench")

load("compare_df")
####
#### we will first get the compare_df !!! #######
####

folder = "/zp1/data/plyu3/MCL_V3"


loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}



equal_matrix <- function(query_mat,ref_mat){
	##########
	genes_overlap = rownames(query_mat)[which(rownames(query_mat) %in% rownames(ref_mat) == T)]
	########## head(rownames(query_mat2))
	m = match(genes_overlap,rownames(query_mat))
	query_mat2 = query_mat[m,]
	########## head(rownames(ref_mat2))
	m = match(genes_overlap,rownames(ref_mat))
	ref_mat2 = ref_mat[m,]
	##########
	########## OK !!!! Then normalize !!!! ########
	##########
	ref_mat2 = round(ref_mat2,3)
	query_mat2 = round(query_mat2,3)
	##########
	ref_mat_input3 = exp(ref_mat2)-1
	#print(head(colSums(ref_mat_input3)))
	query_mat_input3 = exp(query_mat2)-1
	#print(colSums(query_mat_input3))
	##########
	combined_mat = cbind(query_mat_input3,ref_mat_input3)
	#######
	combined_mat_norm = limma::normalizeBetweenArrays(combined_mat,method="quantile")
	combined_mat_norm = log(combined_mat_norm+1)
	combined_mat_norm = round(combined_mat_norm,2)
	#######
	query_mat_input4 = combined_mat_norm[,1:dim(query_mat_input3)[2]]
	ref_mat_input4 = combined_mat_norm[,-c(1:dim(query_mat_input3)[2])]
	#######
	outputlist = list(query_mat_input4,ref_mat_input4)
	return(outputlist)
}

selection_DEGs <- function(all_used_genes,ref_marker,Top=50){
	#### class(ref_marker) ######
	res = apply(ref_marker,2,function(x) length(which(x %in% all_used_genes == T)))
	num = min(res)
	####
	num_new = min(num,Top)
	#### OK!! Then see the Top genes !!!! #######
	ref_marker_cl = ref_marker[1:num_new,]
	####
	ref_marker_cl = reshape2::melt(as.matrix(ref_marker_cl))
	####
	all_markers = ref_marker_cl$value
	all_markers = all_markers[!duplicated(all_markers)]
	####
	return(all_markers)
}


###### we use cor.fk to calculate the correlations ######
calculate_Cor <- function(query_mat_input,ref_mat_input,DEGs_overlap){
	k1 = which(DEGs_overlap %in% rownames(query_mat_input) == T)
	k2 = which(DEGs_overlap %in% rownames(ref_mat_input) == T)
	k3 = k1[which(k1 %in% k2 == T)]
	DEGs_overlap = DEGs_overlap[k3]
	######
	query_mat_input_cl = query_mat_input[which(rownames(query_mat_input) %in% DEGs_overlap == T),]
	ref_mat_input_cl = ref_mat_input[which(rownames(ref_mat_input) %in% DEGs_overlap == T),]
	######
	######
	m1 = match(DEGs_overlap,rownames(query_mat_input_cl))
	m2 = match(DEGs_overlap,rownames(ref_mat_input_cl))
	query_mat_input_cl= query_mat_input_cl[m1,]
	ref_mat_input_cl= ref_mat_input_cl[m2,]
	######
	merge_mat = cbind(query_mat_input_cl,ref_mat_input_cl)
	Cor_res <- pcaPP::cor.fk(merge_mat)
	###### split the Cor_res #######
	query_dim = dim(query_mat_input_cl)[2]
	ref_dim = dim(ref_mat_input_cl)[2]
	######
	Cor_res = Cor_res[,-c(1:query_dim)]
	Cor_res = Cor_res[c(1:query_dim),]
	###### Then we output the most largest clusters ########
	return(Cor_res)
}

######
###### Get hightest correlated cells !!!! #################
######
Res_mat_highest_celltype <- function(res_mat,cutoff){
	res_list = list()
	for(i in 1:dim(res_mat)[1]){
		res_mat_tmp = res_mat[i,]
		k = which(res_mat_tmp == max(res_mat_tmp))
		max_cor = res_mat_tmp[k]
		##### #######
		k2 = which(res_mat_tmp <= max_cor & res_mat_tmp >= cutoff)
		if(length(k2) == 0){
			res_mat_tmp_k2 = "Unassigned"
			res_list = c(res_list,list("Unassigned"))
		}
		if(length(k2) > 0){
		##### #######
			res_mat_tmp_k2 = res_mat_tmp[k2]
			res_mat_tmp_k2 = sort(res_mat_tmp_k2,decreasing=T)
			#####
			if(length(res_mat_tmp_k2) >3){
				res_mat_tmp_k2 = res_mat_tmp_k2[1:3]
			}
			res_list = c(res_list,list(names(res_mat_tmp_k2)))
		}
		######
		######
	}
	names(res_list) = rownames(res_mat)
	return(res_list)
}


########
Analysis_cor <- function(cor_res,lower_cutoff = 0.4){
	####
	# lower_cutoff = 0.35
	###### print the max of cor_res ######
	print(apply(cor_res,1,max))
	######
	cor_resv = as.vector(cor_res)
	cor_resv = sort(cor_resv,decreasing=T)
	######
	model <- mclust::densityMclust(cor_resv,G=1:3)
	###### First we need to know how many models !!!!#########
	number_model = length(levels(as.factor(model$classification)))
	######
	###### Then we get the parameters for each model !!!! ####
	######
	model_mean_total = model$parameters$mean
	model_sd_total = model$parameters$variance$sigmasq
	######
	if(length(model_sd_total) == 1){
		model_sd_total = rep(model_sd_total,number_model)
	}
	###### OK!!! Next we find the cutoffs ########
	if(number_model == 3){
		#### we selected to 2!!! #####
		#### we will find the sencond clusters ####
		tmp_mean = model_mean_total[2]
		tmp_sd = model_sd_total[2]
		cutoff = qnorm(0.75,mean=tmp_mean,sd=sqrt(tmp_sd))
	}
	#######
	if(number_model == 2){
		#### we selected to 2!!! #####
		#### we use the cutoff between the 2 peaks !!!! ##########
		model_classification = as.numeric(model$classification)
		k = which(model$classification %in% model$classification[1] == F)
		index = max(cor_resv[k])
		cutoff = index
	}
	if(number_model == 1){
		#### we selected to 2!!! #####
		#### we use the cutoff between the 2 peaks !!!! ##########
		tmp_mean = model_mean_total[1]
		tmp_sd = model_sd_total[1]
		cutoff = qnorm(0.75,mean=tmp_mean,sd=sqrt(tmp_sd))
	}
	if(cutoff < lower_cutoff){
		cutoff = lower_cutoff
	}
	#####
	return(cutoff)
}

compared_stat <- function(candidate_align,ref_marker_sub,query_mat_1){
	######
	res_tab = as.character(sapply(candidate_align,function(x) x[[1]]))
	res_tab = gsub('@(.+)','',res_tab)
			for(j in 1:length(candidate_align)){
				#print(j)
				tmp = candidate_align[[j]]
				tmp_index = sapply(strsplit(tmp,split='@'),function(x) x[[1]])
				if(length(levels(as.factor(tmp_index))) == 1){
					next
				}else{
					print('compare!')
					ct_list = list()
					for(ii in 1:length(tmp)){
						#print(paste0("ii=",ii))
						sub_ct = tmp[ii]
						subDEGs = ref_marker_sub
						sub_ct_DEGs_index = which(colnames(subDEGs) == sub_ct)
						#####
						sub_ct_DEGs = subDEGs[,sub_ct_DEGs_index]
						#####
						index_1 = which(colnames(query_mat_1) == names(candidate_align)[j])
						#####
						sub_query_mat_input = query_mat_1[,index_1]
						query_sub_ct = sub_query_mat_input[which(names(sub_query_mat_input) %in% sub_ct_DEGs == T)]
						ct_list = c(ct_list,list(query_sub_ct))
					}
					#### then get DEGs #####
					indexJ = Compare3_corr(ct_list,tmp_index)
					res_tab[j] = indexJ
				}
			}
	return(res_tab)
}

Compare3_corr <- function(ct_list,tmp_index){
	######
	empty_matrix = matrix(0,nrow=length(ct_list),ncol=length(ct_list))
	######
	for(i in 1:length(ct_list)){
		for(j in 1:length(ct_list)){
			v_i = ct_list[[i]]
			v_j = ct_list[[j]]
			#### i vs j #########
			res1 = wilcox.test(v_i, v_j, alternative = "greater")
			res2 = wilcox.test(v_i, v_j, alternative = "less")
			if(res1$p.value < 0.05 & res2$p.value > 0.05){
				empty_matrix[i,j] = 1
			}
			if(res1$p.value > 0.05 & res2$p.value < 0.05){
				empty_matrix[i,j] = -1
			}
		}
	}
	matrix_row_max = apply(empty_matrix,1,sum)
	k = which(matrix_row_max == max(matrix_row_max))
	if(length(k) > 1){
		tmp_index1 = tmp_index[k]
		tmp_index1 = tmp_index1[!duplicated(tmp_index1)]
		tmp_index1 = paste(tmp_index1,collapse=' & ')
	}else{
		tmp_index1 = tmp_index
	}
	######
	return(tmp_index1)
}

########
Main_compare_process_CellAnn_time <- function(compare_df,folder=folder,time=T){
	times_list = list()
	mem_list = list()
	##########
	query_label_list = list()
	ref_label_list = list()
	out_table_list = list()
	##########
	for(i in 1:dim(compare_df)[1]){
		t1<-Sys.time()
		print(paste('NOOOOOOO',i))
		print(paste('query:',compare_df$Query_new[i],'  ','ref:',compare_df$Ref_new[i],sep=''))
		TAG = paste('query:',compare_df$Query_new[i],'-->','ref:',compare_df$Ref_new[i],sep='')
		#######
		query = compare_df$Query_new[i]
		ref = compare_df$Ref_new[i]
		#######
		setwd(folder)
		#######
		tmp_functions <- function(query=query,ref=ref,TAG=TAG){
			####### first load the query_Avg expression matrix !!!!! ########
			query_avg = paste0(query,'_Cluster_Avg') ### red ####
			query_mat = loadRData(query_avg)
			print("loadRData(query_avg)")
			#print(bench::mark(loadRData(query_avg))$mem_alloc)
			#### we will use the threshold to 0.7 #####
			#### then we load the query labels !!!! ###
			#query_label_index = paste0(query,'_GroundTruth_Cluster')
			#query_label_tab = loadRData(query_label_index)
			#m1 = match(colnames(query_mat),query_label_tab$cluster)
			#query_label_tab = query_label_tab[m1,]
			####### OK!!! Next !!!! ###################
			####### then we load Refenece input matrix !!!! #####
			ref_avg = paste0(ref,'_SubCluster_Avg')
			ref_mat = loadRData(ref_avg)
			print(paste('Ref_clusters::',dim(ref_mat)))
			#print(bench::mark(loadRData(ref_avg))$mem_alloc)
			#######
			ref_marker_index = paste0(ref,'_CT_Marker')
			ref_marker = loadRData(ref_marker_index)
			#print(bench::mark(loadRData(ref_marker_index))$mem_alloc)
			#######
			ref_marker_subindex = paste0(ref,'_SubCT_Marker')
			ref_marker_sub = loadRData(ref_marker_subindex)
			#######
			ref_label = colnames(ref_mat)
			ref_label = sapply(strsplit(ref_label,split="@"),function(x) x[[1]])
			####### OK!!!! #######
			#######
			####### Let us calculate the correlations !!!!! ################
			####### first we need to make the 2 matrix equal !!!! ##########
			Mat_list = equal_matrix(query_mat,ref_mat)
			#print(bench::mark(equal_matrix(query_mat,ref_mat))$mem_alloc)
			#######
			query_mat_1 = Mat_list[[1]]
			ref_mat_1 = Mat_list[[2]]
			#######
			all_used_genes <- rownames(query_mat_1)
			all_used_DEGs <- selection_DEGs(all_used_genes,ref_marker,Top=100)
			####### Next we need to calculate correlations !!!!! ###########
			####### let us try to use combined markers !!! #################
			##########
			#all_used_DEGs_total = c(all_used_DEGs,all_used_DEGs_sub)
			#all_used_DEGs_total = all_used_DEGs_total[!duplicated(all_used_DEGs_total)]
			##########
			cor_res = calculate_Cor(query_mat_input=query_mat_1,ref_mat_input=ref_mat_1,DEGs_overlap=all_used_DEGs)
			########## Next we need the cutoffs !!!! #######################
			cutoff = Analysis_cor(cor_res,lower_cutoff = 0.4)
			########## Next we get the highest correlated cells !!!!! ######
			candidate_align = Res_mat_highest_celltype(cor_res,cutoff)
			##########
			##########
			########## "red" ####################
			########## load DEGs #######
			##########
			res_max = apply(cor_res,1,max)
			Unassigned_index = which(res_max < cutoff)
			#####
			##### Then we compare the DEGs !!! ###########
			#####
			res_tab = compared_stat(candidate_align,ref_marker_sub,query_mat_1)
			#####
			#####
			if(length(Unassigned_index) > 0){
				res_tab[Unassigned_index] = 'Unassigned'
			}
			res_table = data.frame(cluster=colnames(query_mat_1),result=res_tab)
			res_table = list(res_table)
			names(res_table) = TAG
		}
		#######
		###### we need match the order of the avg matrix #######
 		#####
 		##### t2<-Sys.time()
		res = bench::mark(tmp_functions(query=query,ref=ref,TAG=TAG))
		res1 = list(res$mem_alloc)
		names(res1) = TAG
		mem_list = c(mem_list,res1)
		#########
		res2 = list(res$total_time)
		names(res2) = TAG
		times_list = c(times_list,res2)
	}
	return(list(times_list,mem_list))
}

########## 

CellAnn_MEM = Main_compare_process_CellAnn_time(compare_df,folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare",time=T)
CellAnn_MEM_1 = CellAnn_MEM[[1]]
CellAnn_MEM_2 = CellAnn_MEM[[2]]
save(CellAnn_MEM_1,file="CellAnn_time")
save(CellAnn_MEM_2,file="CellAnn_mem")


library("bench")

setwd("/zp1/data/plyu3/MCL_V3")

load("scmap_cluster_MEM")

load("CellAnn_MEM")

load("seurat_MEM")

load("chetah_MEM")
######### OK!!! Next we need to integrate these 2 functions !!! ########
#########
#########
save(scPred_MEM[[1]],file="scPred_time")
save(scPred_MEM[[2]],file="scPred_mem")



tmp_MEM = CellAnn_MEM

##########
##########

CellAnn_mems_res = Convert_MEM_list_to_data_frame(CellAnn_MEM,'CellAnn')

scmap_cluster_res = Convert_MEM_list_to_data_frame(scmap_cluster_MEM,'scmap_cluster')

chetah_mems_res = Convert_MEM_list_to_data_frame(chetah_MEM,'chetah')

seurat_mems_res = Convert_MEM_list_to_data_frame(seurat_MEM,'Seurat')

###############
###############

combined_plot_res = rbind(CellAnn_mems_res,scmap_cluster_res,seurat_mems_res,chetah_mems_res)

combined_plot_res$cells = as.numeric(combined_plot_res$cells)

library(ggplot2)

ggplot(combined_plot_res,aes(x=cells,y=log10(mems+1))) + geom_point(aes(color=tag)) + geom_smooth(method = "lm",formula = y ~ poly(log10(x+1), 2),se = FALSE,aes(color=tag)) + theme_classic() + scale_x_continuous(breaks=c(1000,5000,10000,20000,40000,60000,80000,100000)) + geom_hline(yintercept=c(log10(1+1),log10(10+1),log10(100+1)),color="grey",linetype='dashed')+ theme(axis.text.x=element_text(angle = 45, hjust = 0.5, vjust = -0.1))

ggsave("test2.png",width=6,height=3)

##############
##############
############## Next we will integrate all the results !!! ######
##############
##############


cd /zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare

ls -lh

conda activate CellAnn_test
R

library(Seurat)


#### first we will load the 


conda activate CellAnn_test
R

setwd("/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")

library(bench)

loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}


CellAnn_mem <- loadRData("CellAnn_mem")
CellAnn_time <- loadRData("CellAnn_time")

##### next ########

Seurat_mem <- loadRData("seurat_mem")
Seurat_time <- loadRData("seurat_time")


##### next ########

Scmap_cluster_time <- loadRData("scmap_cluster_time")
Scmap_cluster_mem <- loadRData("scmap_cluster_mem")

##### next ########
chetah_MEM_12 <- loadRData("chetah_mem_12")
chetah_MEM_3 <- loadRData("chetah_mem_3")
chetah_MEM_4 <- loadRData("chetah_mem_4")
chetah_MEM_5 <- loadRData("chetah_mem_5")
chetah_MEM_6 <- loadRData("chetah_mem_6")
chetah_MEM_7 <- loadRData("chetah_mem_7")
chetah_MEM_8 <- loadRData("chetah_mem_8")

chetah_MEM <- c(chetah_MEM_12,chetah_MEM_3,chetah_MEM_4,chetah_MEM_5,chetah_MEM_6,chetah_MEM_7,chetah_MEM_8)

save(chetah_MEM,file='chetah_MEM')

##### 

chetah_time_12 <- loadRData("chetah_time_12")
chetah_time_3 <- loadRData("chetah_time_3")
chetah_time_4 <- loadRData("chetah_time_4")
chetah_time_5 <- loadRData("chetah_time_5")
chetah_time_6 <- loadRData("chetah_time_6")
chetah_time_7 <- loadRData("chetah_time_7")
chetah_time_8 <- loadRData("chetah_time_8")

chetah_time <- c(chetah_time_12,chetah_time_3,chetah_time_4,chetah_time_5,chetah_time_6,chetah_time_7,chetah_time_8)

save(chetah_time,file='chetah_time')



##### Next !!! ####
scClassify_mem_13 <- loadRData("scClassify_mem_13")
scClassify_mem_4 <- loadRData("scClassify_mem_4")
scClassify_mem_5 <- loadRData("scClassify_mem_5")
scClassify_mem_6 <- loadRData("scClassify_mem_6")
scClassify_mem_7 <- loadRData("scClassify_mem_7")
scClassify_mem_8 <- loadRData("scClassify_mem_8")

scClassify_mem <- c(scClassify_mem_13,scClassify_mem_4,scClassify_mem_5,scClassify_mem_6,scClassify_mem_7,scClassify_mem_8)

save(scClassify_mem,file='scClassify_mem')

##### Next !!! #####
scClassify_time_13 <- loadRData("scClassify_time_13")
scClassify_time_4 <- loadRData("scClassify_time_4")
scClassify_time_5 <- loadRData("scClassify_time_5")
scClassify_time_6 <- loadRData("scClassify_time_6")
scClassify_time_7 <- loadRData("scClassify_time_7")
scClassify_time_8 <- loadRData("scClassify_time_8")

scClassify_time <- c(scClassify_time_13,scClassify_time_4,scClassify_time_5,scClassify_time_6,scClassify_time_7,scClassify_time_8)

save(scClassify_time,file='scClassify_time')

##### Next !!! #####

scPred_mem_13 <- loadRData("scPred_mem_13")
scPred_mem_4 <- loadRData("scPred_mem_4")
scPred_mem_5 <- loadRData("scPred_mem_5")
scPred_mem_6 <- loadRData("scPred_mem_6")

scPred_mem <- c(scPred_mem_13,scPred_mem_4,scPred_mem_5,scPred_mem_6)
##### scPred 6 and 7 !!!! #########
##### Next !!! #####

scPred_time_13 <- loadRData("scPred_time_13")
scPred_time_4 <- loadRData("scPred_time_4")
scPred_time_5 <- loadRData("scPred_time_5")
scPred_time_6 <- loadRData("scPred_time_6")

scPred_time <- c(scPred_time_13,scPred_time_4,scPred_time_5,scPred_time_6)


#####
##### Next we will merge all the samples !!! #######
#####
Convert_MEM_list_to_data_frame <- function(tmp_MEM,tag='cellann'){
	#########
	mems = c()
	#########
	for(i in 1:length(tmp_MEM)){
		tmptmp = tmp_MEM[[i]]
		tmptmp = as.character(tmptmp)
		if(length(grep("MB",tmptmp)) == 1){
			print("convert MB to GB !!!")
			tmptmp_number = gsub("MB","",tmptmp)
			tmptmp_number = as.numeric(tmptmp_number)
			tmptmp_number = tmptmp_number/1024
		}
		if(length(grep("GB",tmptmp)) == 1){
			#print("convert MB to GB !!!")
			tmptmp_number = gsub("GB","",tmptmp)
			tmptmp_number = as.numeric(tmptmp_number)
		}
		mems = c(mems,tmptmp_number)
	}
	options(scipen=10)
	#########
	index = c(100,500,1000,5000,10000,15000,20000,25000)
	index = as.character(index)
	#########
	index_cl = index[1:length(tmp_MEM)]
	#########
	dat = data.frame(cells=index_cl,mems=mems,tag=tag)
	#########
	#########
	return(dat)
}


CellAnn_mems_res = Convert_MEM_list_to_data_frame(CellAnn_mem,'CellAnn')
Scmap_cluster_mems_res = Convert_MEM_list_to_data_frame(Scmap_cluster_mem,'Scmap_cluster')
chetah_mems_res = Convert_MEM_list_to_data_frame(chetah_MEM,'Chetah')
Seurat_mems_res = Convert_MEM_list_to_data_frame(Seurat_mem,'Seurat')
scPred_mems_res = Convert_MEM_list_to_data_frame(scPred_mem,'scPred')
scClassify_mems_res = Convert_MEM_list_to_data_frame(scClassify_mem,'scClassify')

####### #####
#######

combined_plot_res = rbind(CellAnn_mems_res,Scmap_cluster_mems_res,chetah_mems_res,Seurat_mems_res,scPred_mems_res,scClassify_mems_res)

combined_plot_res$cells = as.numeric(combined_plot_res$cells)

library(ggplot2)

ggplot(combined_plot_res,aes(x=cells,y=log10(mems+1))) + geom_point(aes(color=tag)) + geom_smooth(method = "lm",formula = y ~ poly(log10(x+1), 2),se = FALSE,aes(color=tag)) + theme_classic() + scale_x_continuous(breaks=c(1000,5000,10000,15000,20000,25000)) + geom_hline(yintercept=c(log10(1+1),log10(10+1),log10(100+1)),color="grey",linetype='dashed')+ theme(axis.text.x=element_text(angle = 45, hjust = 0.5, vjust = -0.1))

ggsave("test2.png",width=6,height=3)


####### Next is the time !!! ##########
#######
#######

tmp_time = CellAnn_time

Convert_time_list_to_data_frame <- function(tmp_time,tag='cellann'){
	#########
	mems = c()
	#########
	for(i in 1:length(tmp_time)){
		print(i)
		tmptmp = tmp_time[[i]]
		tmptmp = as.character(tmptmp)
		if(length(grep("m$",tmptmp)) == 1){
			print("convert m to s !!!")
			tmptmp_number = gsub("m","",tmptmp)
			tmptmp_number = as.numeric(tmptmp_number)
			tmptmp_number = tmptmp_number*60
		}
		if(length(grep("ms$",tmptmp)) == 1){
			print("convert ms to s !!!")
			tmptmp_number = gsub("ms","",tmptmp)
			tmptmp_number = as.numeric(tmptmp_number)
			tmptmp_number = round(tmptmp_number/1000,3)
		}
		if(length(grep("s",tmptmp)) == 1 & length(grep("ms",tmptmp)) == 0){
			#print("convert MB to GB !!!")
			tmptmp_number = gsub("s","",tmptmp)
			tmptmp_number = as.numeric(tmptmp_number)
		}
		mems = c(mems,tmptmp_number)
	}
	options(scipen=10)
	#########
	index = c(100,500,1000,5000,10000,15000,20000,25000)
	index = as.character(index)
	#########
	index_cl = index[1:length(tmp_time)]
	#########
	dat = data.frame(cells=index_cl,mems=mems,tag=tag)
	#########
	#########
	return(dat)
}

############### #################

CellAnn_time_res = Convert_time_list_to_data_frame(CellAnn_time,'CellAnn')
Scmap_cluster_time_res = Convert_time_list_to_data_frame(Scmap_cluster_time,'Scmap_cluster')
chetah_time_res = Convert_time_list_to_data_frame(chetah_time,'Chetah')
Seurat_time_res = Convert_time_list_to_data_frame(Seurat_time,'Seurat')
scPred_time_res = Convert_time_list_to_data_frame(scPred_time,'scPred')
scClassify_time_res = Convert_time_list_to_data_frame(scClassify_time,'scClassify')


##### ####
#####


combined_plot_res = rbind(CellAnn_time_res,Scmap_cluster_time_res,chetah_time_res,Seurat_time_res,scPred_time_res,scClassify_time_res)

combined_plot_res$cells = as.numeric(combined_plot_res$cells)

library(ggplot2)

ggplot(combined_plot_res,aes(x=cells,y=log10(mems+1))) + geom_point(aes(color=tag)) + geom_smooth(method = "lm",formula = y ~ poly(log10(x+1), 2),se = FALSE,aes(color=tag)) + theme_classic() + scale_x_continuous(breaks=c(1000,5000,10000,15000,20000,25000)) + geom_hline(yintercept=c(log10(10+1),log10(60+1),log10(3600+1)),color="grey",linetype='dashed')+ theme(axis.text.x=element_text(angle = 45, hjust = 0.5, vjust = -0.1))

ggsave("test3.png",width=6,height=3)

#####
##### ####
##### count the clusters ########
#####













