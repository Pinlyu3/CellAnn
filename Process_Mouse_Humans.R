########
########

### on the server ###

########
########

### load the datasets ###

########
########

setwd("/zp1/data/plyu3/M_G_expression_datasets")

########
########

ssh [plyu3@omb2.onc.jhmi.edu](mailto:plyu3@omb2.onc.jhmi.edu)

plyu3 U[9C20&&

ssh plyu3@omb2.onc.jhmi.edu

########
########

conda activate seurat4

R

setwd('/zp1/data/plyu3/M_G_expression_datasets')

library("rtracklayer")

Human_annotations = import('gencode.v41.basic.annotation.gtf.gz')

Human_gene_type = Human_annotations$gene_type

Human_gene_name = Human_annotations$gene_name

Human_gene_tab = data.frame(Human_gene_type,Human_gene_name)

Human_gene_tab$index = paste(Human_gene_tab$Human_gene_type,Human_gene_tab$Human_gene_name)

Human_gene_tab_cl = Human_gene_tab[!duplicated(Human_gene_tab$index),]

table(Human_gene_tab_cl$Human_gene_type)

#### we only use the protein_coding genes #########

####

Human_gene_tab_clcl = Human_gene_tab_cl[which(Human_gene_tab_cl$Human_gene_type == 'protein_coding'),]

####

Human_gene_list = Human_gene_tab_clcl$Human_gene_name

#### rm the MT and rm the ENSG genes #######

k1 = grep("^MT",Human_gene_list)
k2 = grep("^ENSG",Human_gene_list)

####
Human_gene_list = Human_gene_list[c(-k1,-k2)]
length(Human_gene_list)

####
saveRDS(Human_gene_list,file='Human_gene_list')



#### OK next is the Mouse list ###########


setwd('/zp1/data/plyu3/M_G_expression_datasets')

library("rtracklayer")

Mouse_annotations = import('gencode.vM30.basic.annotation.gtf.gz')

Mouse_gene_type = Mouse_annotations$gene_type

Mouse_gene_name = Mouse_annotations$gene_name

Mouse_gene_tab = data.frame(Mouse_gene_type,Mouse_gene_name)

Mouse_gene_tab$index = paste(Mouse_gene_tab$Mouse_gene_type,Mouse_gene_tab$Mouse_gene_name)

Mouse_gene_tab_cl = Mouse_gene_tab[!duplicated(Mouse_gene_tab$index),]

table(Mouse_gene_tab_cl$Mouse_gene_type)

##### OK !!! Then we filtered the ############
#### we only use the protein_coding genes #########

####

Mouse_gene_tab_clcl = Mouse_gene_tab_cl[which(Mouse_gene_tab_cl$Mouse_gene_type == 'protein_coding'),]

####
Mouse_gene_list = Mouse_gene_tab_clcl$Mouse_gene_name

#### rm the MT and rm the ENSG genes #######

k1 = grep("^mt",Mouse_gene_list)

####

Mouse_gene_list = Mouse_gene_list[c(-k1)]

length(Mouse_gene_list)

####

saveRDS(Mouse_gene_list,file='Mouse_gene_list')

#### #### ####






















