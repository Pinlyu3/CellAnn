
#### load the library ###############
library(Seurat)
setwd("/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
library("bench")
load("compare_df")

library(CHETAH)
library(scPred)
library(scmap)
library(scClassify)

####


loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}

####
#### OK!!!! #####
CellAnn_chetah = function(train,test,label_train,time = F){
  ##########
  library(SingleCellExperiment)
  library(CHETAH)
  ###########
  start_time = Sys.time()
  ###########
  sce = SingleCellExperiment(assays = list(counts = train),colData = data.frame(celltypes = label_train))
  sce_test = SingleCellExperiment(assays = list(counts = test))
  ###########
  sce_test = CHETAHclassifier(input = sce_test, ref_cells = sce)
  ###########
  tmp_tab = sce_test@colData
  ###########
  predict_label = unname(tmp_tab$celltype_CHETAH)
  ########### Then we may found the node !!! ##########
  ### PlotCHETAH(input = sce_test, interm = TRUE)
  nodes = sce_test@int_metadata$CHETAH$nodetypes
  #### 
  k = grep("Node",predict_label)
  #### we change node to cell types ##########
  if(length(k) > 0){
  		for(ki in 1:length(k)){
  			tmp_node = predict_label[k[ki]]
  			tmp_node_index = as.numeric(gsub('Node','',tmp_node)) + 1
  			tmp_ct = names(nodes[[tmp_node_index]])
  			tmp_ct = paste(tmp_ct,collapse=' & ')
  			predict_label[k[ki]] = tmp_ct
  		}
  }
  ############
  end_time = Sys.time()
  ###########
  times = as.numeric(difftime(end_time,start_time,units = 'secs'))
  ###########
  if(time){
    return(list(predict_label = predict_label,times = times))
  }
  return(predict_label)
}
#### Then Next tools !!!! #################
CellAnn_seurat = function(train,
                  test,
                  label_train,
                  k.filter = NA,
                  time = T,
                  selection.method = 'vst',
                  nfeatures = 2000,
                  mean.cutoff = c(0.1, 8),
                  dispersion.cutoff = c(1, Inf),
                  prediction.score.max = 0.5){
  #######
  start_time = Sys.time()
  ########
  reference = CreateSeuratObject(train)
  reference = NormalizeData(reference,verbose = F)
  reference = FindVariableFeatures(reference,mean.cutoff = mean.cutoff,dispersion.cutoff = dispersion.cutoff,verbose = F)
  reference$celltype = label_train
  query = CreateSeuratObject(test)
  query = NormalizeData(query,verbose = F)
  query = FindVariableFeatures(query,mean.cutoff = mean.cutoff,dispersion.cutoff = dispersion.cutoff,verbose = F)
  #######
  anchors = FindTransferAnchors(reference,query,k.filter = NA)
  k.weight = dim(anchors@anchors)[1]
  predictions = try(TransferData(anchors,as.character(reference$celltype)))
  while(inherits(predictions,'try-error')){
  	if(k.weight > 100){
  		k.weight = 100
  	}
  	k.weight=k.weight-1
  	print(k.weight)
  	predictions = try(TransferData(anchors,as.character(reference$celltype),k.weight=k.weight))
  }
  #######
  query = AddMetaData(query, metadata = predictions)
  ######
  print(summary(query$prediction.score.max))
  query$predicted.id[which(query$prediction.score.max < prediction.score.max)] <- 'Unassigned'
  #######
  predict_label = unname(query$predicted.id)
  ######
  end_time = Sys.time()
  ######
  times = as.numeric(difftime(end_time,start_time,units = 'secs'))
  ######
  if(time){
    return(list(predict_label = predict_label,times = times))
  }
  return(predict_label)
}
#### OK!!!!! ### The Next tools !!!! ######
####
CellAnn_scpred = function(train,
                  test,
                  label_train,
                  model = 'svmRadial',
                  reclassify = NULL,
                  time = F,
                  threshold = 0.55){

  library(scPred)
  library(dplyr)
  start_time = Sys.time()
  #####
  print(summary(rowSums(train)))
  print(summary(rowSums(test)))
  ######
  reference = CreateSeuratObject(train,min.cells = 10,min.features = 10)
  query = CreateSeuratObject(test,min.cells = 10,min.features = 10)
  reference = reference %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30)
  ###########
  query = NormalizeData(query)
  reference$cell_type = label_train
  #####
  k = which(table(reference$cell_type) < 20)
  if(length(k) > 0){
  	print(names(table(reference$cell_type))[k])
  	k2 = which(reference$cell_type %in% names(table(reference$cell_type))[k] == T)
  	reference = reference[,-k2]
  	reference$cell_type = as.character(reference$cell_type)
  	print(table(reference$cell_type))
  }
  ##### ####### ############## #################################### red ####### change to seurat_clusters should be celltype #####
  reference = getFeatureSpace(reference, "cell_type")
  reference = trainModel(reference)
  query = scPredict(query,
                    reference,
                    threshold = threshold)
  predict_label = unname(query$scpred_prediction)
  end_time = Sys.time()
  times = as.numeric(difftime(end_time,start_time,units = 'secs'))
  if(time){
    return(list(predict_label = predict_label,times = times))
  }
  return(predict_label)
}
#### OK!!! #### Next !!!! #################
CellAnn_scClassify <- function(test = test,
                  train = train,
                  label_train = label_train,
                  time=T,
                  prob_threshold=0.5
                  ){
	#### first train the model ####
	####
	scClassify_res_ensemble <- scClassify(exprsMat_train = train,
                                      cellTypes_train = label_train,
                                      exprsMat_test = test,
                                      tree = "HC",
                                      algorithm = "WKNN",
                                      selectFeatures = c("limma"),
                                      similarity = c("pearson", "cosine"),
                                      weighted_ensemble = FALSE,
                                      returnList = FALSE,
                                      verbose = FALSE)
	####
	####
	start_time = Sys.time()
	pred_res <- scClassify_res_ensemble$testRes$test$ensembleRes$cellTypes
	#####
	end_time = Sys.time()
	times = as.numeric(difftime(end_time,start_time,units = 'secs'))
	if(time){
    	return(list(predict_label = pred_res$ensembleRes,times = times))
  	}
	return(pred_res)
}
####
CellAnn_scmapcluster = function(train,
                        test,
                        label_train,
                        threshold = 0.7,
                        time = T){
  start_time = Sys.time()
  ######
  sce = SingleCellExperiment(list(counts = train),colData = data.frame(cell_type1 = label_train))
  logcounts(sce) = log2(counts(sce) + 1)
  rowData(sce)$feature_symbol = rownames(sce)
  sce = selectFeatures(sce)
  ######
  sce_test = SingleCellExperiment(list(counts = test))
  logcounts(sce_test) = log2(counts(sce_test) + 1)
  rowData(sce_test)$feature_symbol = rownames(sce_test)
  ######
  sce = indexCluster(sce)
  scmapCluster_results = scmapCluster(projection = sce_test,index_list = list(sce@metadata$scmap_cluster_index),threshold = threshold)
  predict_label = scmapCluster_results$combined_labs
  #######
  end_time = Sys.time()
  #######
  times = as.numeric(difftime(end_time,start_time,units = 'secs'))
  #######
  if(time){
    return(list(predict_label = predict_label,times = times))
  }
  ########
  return(predict_label)
}


####

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

#######
####### we need find the reason ######
#######


#### scmap !!! #################
print("scmap_cluster_MEM")
scmap_cluster_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df,method='scmap-cluster',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
scmap_cluster_MEM_1 = scmap_cluster_MEM[[1]]
scmap_cluster_MEM_2 = scmap_cluster_MEM[[2]]
save(scmap_cluster_MEM_1,file="scmap_cluster_time")
save(scmap_cluster_MEM_2,file="scmap_cluster_mem")

#### OK!!! Next tools !!!!! ####
print("chetah_MEM")
chetah_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[1:2,],method='chetah',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
chetah_time_12 = chetah_MEM[[1]]
chetah_mem_12 = chetah_MEM[[2]]
save(chetah_time_12,file="chetah_time_12")
save(chetah_mem_12,file="chetah_mem_12")



chetah_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[3,],method='chetah',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
chetah_time_3 = chetah_MEM[[1]]
chetah_mem_3 = chetah_MEM[[2]]

save(chetah_time_3,file='chetah_time_3')
save(chetah_mem_3,file='chetah_mem_3')


chetah_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[4,],method='chetah',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
chetah_time_4 = chetah_MEM[[1]]
chetah_mem_4 = chetah_MEM[[2]]
save(chetah_time_4,file='chetah_time_4')
save(chetah_mem_4,file='chetah_mem_4')


chetah_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[5,],method='chetah',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")

chetah_time_5 = chetah_MEM[[1]]
chetah_mem_5 = chetah_MEM[[2]]

save(chetah_time_5,file='chetah_time_5')
save(chetah_mem_5,file='chetah_mem_5')



chetah_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[6,],method='chetah',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
chetah_time_6 = chetah_MEM[[1]]
chetah_mem_6 = chetah_MEM[[2]]

save(chetah_time_6,file='chetah_time_6')
save(chetah_mem_6,file='chetah_mem_6')

chetah_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[8,],method='chetah',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")

chetah_time_7 = chetah_MEM[[1]]
chetah_mem_7 = chetah_MEM[[2]]

save(chetah_time_7,file='chetah_time_7')
save(chetah_mem_7,file='chetah_mem_7')


chetah_time_8 = chetah_MEM[[1]]
chetah_mem_8 = chetah_MEM[[2]]

save(chetah_time_8,file='chetah_time_8')
save(chetah_mem_8,file='chetah_mem_8')




#### OK!!! Next tools !!!! #####
print("seurat_MEM")
seurat_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df,method='seurat',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
seurat_MEM_1 = seurat_MEM[[1]]
seurat_MEM_2 = seurat_MEM[[2]]
save(seurat_MEM_1,file="seurat_time")
save(seurat_MEM_2,file="seurat_mem")


#### OK!!! Next tools !!!!! #####
print("scPred_MEM")
scPred_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[1:3,],method='scpred',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
scPred_MEM_1 = scPred_MEM[[1]]
scPred_MEM_2 = scPred_MEM[[2]]
save(scPred_MEM_1,file="scPred_time_13")
save(scPred_MEM_2,file="scPred_mem_13")


scPred_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[4,],method='scpred',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
scPred_MEM_1 = scPred_MEM[[1]]
scPred_MEM_2 = scPred_MEM[[2]]
save(scPred_MEM_1,file="scPred_time_4")
save(scPred_MEM_2,file="scPred_mem_4")



scPred_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[5,],method='scpred',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
scPred_MEM_1 = scPred_MEM[[1]]
scPred_MEM_2 = scPred_MEM[[2]]
save(scPred_MEM_1,file="scPred_time_5")
save(scPred_MEM_2,file="scPred_mem_5")


scPred_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[6,],method='scpred',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
scPred_MEM_1 = scPred_MEM[[1]]
scPred_MEM_2 = scPred_MEM[[2]]
save(scPred_MEM_1,file="scPred_time_6")
save(scPred_MEM_2,file="scPred_mem_6")


#####
#####


scPred_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[7,],method='scpred',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
scPred_MEM_1 = scPred_MEM[[1]]
scPred_MEM_2 = scPred_MEM[[2]]
save(scPred_MEM_1,file="scPred_time_7")
save(scPred_MEM_2,file="scPred_mem_7")


scPred_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[8,],method='scpred',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
scPred_MEM_1 = scPred_MEM[[1]]
scPred_MEM_2 = scPred_MEM[[2]]
save(scPred_MEM_1,file="scPred_time_8")
save(scPred_MEM_2,file="scPred_mem_8")



save(scPred_MEM[[1]],file="scPred_time")
save(scPred_MEM[[2]],file="scPred_mem")


##### OK!!! Next tools !!!! #####
print("scClassify_MEM")
scClassify_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[2:8,],method='scClassify',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
save(scClassify_MEM[[1]],file="scClassify_time")
save(scClassify_MEM[[2]],file="scClassify_mem")



scClassify_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[1:3,],method='scClassify',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
scClassify_MEM_1 = scClassify_MEM[[1]]
scClassify_MEM_2 = scClassify_MEM[[2]]
save(scClassify_MEM_1,file="scClassify_time_13")
save(scClassify_MEM_2,file="scClassify_mem_13")

scClassify_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[4,],method='scClassify',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
scClassify_MEM_1 = scClassify_MEM[[1]]
scClassify_MEM_2 = scClassify_MEM[[2]]
save(scClassify_MEM_1,file="scClassify_time_4")
save(scClassify_MEM_2,file="scClassify_mem_4")


scClassify_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[5,],method='scClassify',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
scClassify_MEM_1 = scClassify_MEM[[1]]
scClassify_MEM_2 = scClassify_MEM[[2]]
save(scClassify_MEM_1,file="scClassify_time_5")
save(scClassify_MEM_2,file="scClassify_mem_5")


scClassify_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[6,],method='scClassify',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
scClassify_MEM_1 = scClassify_MEM[[1]]
scClassify_MEM_2 = scClassify_MEM[[2]]
save(scClassify_MEM_1,file="scClassify_time_6")
save(scClassify_MEM_2,file="scClassify_mem_6")


scClassify_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[7,],method='scClassify',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
scClassify_MEM_1 = scClassify_MEM[[1]]
scClassify_MEM_2 = scClassify_MEM[[2]]
save(scClassify_MEM_1,file="scClassify_time_7")
save(scClassify_MEM_2,file="scClassify_mem_7")


scClassify_MEM = Main_compare_process_OtherTools_Single_Cell_times_MEM(compare_df[8,],method='scClassify',folder="/zp1/data/plyu3/CellAnn_test_AUC/PMBC_prepare")
scClassify_MEM_1 = scClassify_MEM[[1]]
scClassify_MEM_2 = scClassify_MEM[[2]]
save(scClassify_MEM_1,file="scClassify_time_8")
save(scClassify_MEM_2,file="scClassify_mem_8")




##### OK!!! #####################
#####
print("Done!!!!")






