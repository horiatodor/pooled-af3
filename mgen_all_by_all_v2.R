#########
#analysis of mgen all by all (partial), using the new analysis pipeline. 
#########

########
#reading in the relevant stuff
########

#read in all of the description files and concetenate them
description_files <- paste0("/Volumes/AF3SSD/job_info/", list.files("/Volumes/AF3SSD/job_info/"))
description_files <- description_files[-grep("matrix_of_pairs", description_files)]
list_of_descriptions <- list()

for (i in seq_along(description_files)){
  list_of_descriptions[[i]] <- read.csv(description_files[i], row.names = 1, stringsAsFactors = FALSE)
}

all_descriptions <- do.call("rbind", list_of_descriptions)  
all_descriptions[,1] <- apply(do.call("rbind", strsplit(all_descriptions[,1], split = "__")),1,paste, collapse = "_")

rm(list_of_descriptions)

#size of jobs in amino acids
job_sizes <- rep(NA, dim(all_descriptions)[1])
for (i in seq_along(job_sizes)){
  of_int <- match(strsplit(all_descriptions[i,2], split = "::::")[[1]],names(mgen_proteins))
  job_sizes[i] <- sum(sapply(mgen_proteins[of_int], nchar))
}


#read in the list of proteins
mgen_proteins <- seqinr::read.fasta("mgen_proteins.txt", forceDNAtolower = FALSE, as.string = TRUE, whole.header = TRUE)
names(mgen_proteins) <- do.call("rbind", strsplit(do.call("rbind", strsplit(names(mgen_proteins), split = "locus_tag="))[,2], split = "] "))[,1]

########
#read in all of the IPTMS  
########
iptms_and_paes <- get_list_of_iptm(afolder_of_zip = "/Volumes/AF3SSD/AF3_completed_zips/", 
                                   description = all_descriptions, return_paes = FALSE, return_sd = TRUE)

list_of_iptm <- iptms_and_paes[[1]]
list_of_paes <- iptms_and_paes[[2]]
names(list_of_paes) <- names(list_of_iptm)

temp2 <- resolve_duplicates(list_of_iptm) 
list_of_paes <- resolve_duplicates(list_of_paes)[[1]]

########
#do the size correction 
########

list_of_iptm_size_corr <- size_correct(temp2[[1]], mgen_proteins)[[1]]
list_of_iptm_size_corr_ratio <- size_correct(temp2[[1]], mgen_proteins, type_of_corr = "ratio")[[1]]
list_of_iptm_size_corr_variable <- size_correct(temp2[[1]], mgen_proteins, type_of_corr = "variable")[[1]]
list_of_iptm_size_corr_eff <- size_correct(temp2[[1]], mgen_proteins, type_of_corr = "subtract",protein_sizes = mol_mass)[[1]]


########
#get the replicates and the matrix
########

temp <- make_matrix_all_by_all(list_of_iptm_size_corr, mgen_proteins, resolve_reps = "mean")
temp_ratio <- make_matrix_all_by_all(list_of_iptm_size_corr_ratio, mgen_proteins, resolve_reps = "mean")
temp_variable <- make_matrix_all_by_all(list_of_iptm_size_corr_variable, mgen_proteins, resolve_reps = "mean")
temp_eff <- make_matrix_all_by_all(list_of_iptm_size_corr_eff, mgen_proteins, resolve_reps = "mean")

temp_paes <- make_matrix_all_by_all(list_of_paes, mgen_proteins, resolve_reps = "mean")
temp_uncorrected <- make_matrix_all_by_all(temp2[[1]], mgen_proteins, resolve_reps = "mean")

matrix_of_iptm <- temp[[1]]
matrix_of_iptm_cor <- cor(matrix_of_iptm, use = "pair")
diag(matrix_of_iptm_cor) <- NA

#temp_paes <- make_matrix_all_by_all(list_of_paes, mgen_proteins, resolve_reps = "mean")
matrix_of_iptm[which(temp[[3]] == 1, arr.ind = TRUE)] <- NA

########
#replicate analysis
########
replicate_matrix <- expand_list_of_pairs(temp[[2]])
replicate_matrix_sds <- expand_list_of_pairs(temp_paes[[2]])

plot(replicate_matrix, pch = 16, col = scales::alpha("black", 0.1),
     ylim = c(-0.05,0.8), xlim = c(-0.05,0.8))

#plotting uncertainty
segments(replicate_matrix[,1], replicate_matrix[,2]-replicate_matrix_sds[,2],
         replicate_matrix[,1], replicate_matrix[,2]+replicate_matrix_sds[,2], col = scales::alpha("black", 0.1))

segments(replicate_matrix[,1]-replicate_matrix_sds[,1], replicate_matrix[,2],
         replicate_matrix[,1]+replicate_matrix_sds[,1], replicate_matrix[,2], col = scales::alpha("black", 0.1))

#plotting STRING
replicate_list <- expand_list_of_pairs(temp[[2]] ,return_list = TRUE)
for (i in seq_along(replicate_list)){
  the_ind <- match(strsplit(names(replicate_list)[i], split = "::::")[[1]], rownames(string_mat[[2]]))
  if (string_mat[[2]][the_ind[1], the_ind[2]] > 800){
    points(replicate_list[[i]], col = scales::alpha("red", 0.5))
  }
}

cor(replicate_matrix[,1], replicate_matrix[,2])

########
#benchmark against STRING
########

#######
#string
string_file <- read.delim("243273.protein.links.detailed.v12.0.txt", sep = " ")
string_file[,1] <- do.call("rbind", strsplit(string_file[,1], split = "\\."))[,2]
string_file[,2] <- do.call("rbind", strsplit(string_file[,2], split = "\\."))[,2]
string_names <- apply(string_file[,1:2],1,paste, collapse = "::::")

#string mat
source("~/Documents/BMK analysis/make_matrix_of_string.R")
string_mat <- make_matrix_of_string(rownames(matrix_of_iptm), rownames(matrix_of_iptm), string_file[,c(1,2,10,7)])
string_mat[[1]][which(is.na(string_mat[[1]]), arr.ind = TRUE)] <- 0
string_mat[[2]][which(is.na(string_mat[[2]]), arr.ind = TRUE)] <- 0

#some plots
source("~/Documents/AF3/assess_string_v3.R")
string_assess <- assess_string_v2(matrix_of_iptm, string_mat[[1]], to_test = seq(-0.1,0.7, by = 0.05), string_threshold = 800)
string_assess_exp <- assess_string_v2(matrix_of_iptm, string_mat[[2]], to_test = seq(-0.1,0.7, by = 0.05), string_threshold = 800)

#ROC plots (all interactions and experimental)
plot(assess_string_v2(matrix_of_iptm, string_mat[[1]], to_test = seq(-0.1,0.7, by = 0.001), string_threshold = 800)[,c(6,3)], type = 'l')
points(assess_string_v2(matrix_of_iptm, string_mat[[2]], to_test = seq(-0.1,0.7, by = 0.001), string_threshold = 800)[,c(6,3)], lty = 2, type = 'l')
legend("topleft", lty = c(1,2), legend = c("STRING all > 800","STRING exp > 800"))
abline(a=0,b=1,col = 'red')

#
pROC::auc(string_mat[[1]] > 800, as.numeric(matrix_of_iptm))
pROC::auc(string_mat[[2]] > 800, as.numeric(matrix_of_iptm))

pROC::auc(string_mat[[1]] > 800, as.numeric(matrix_of_iptm_cor))
pROC::auc(string_mat[[2]] > 800, as.numeric(matrix_of_iptm_cor))

#
poop <- apply(simplify2array(list(matrix_of_iptm, matrix_of_iptm_cor)), c(1, 2), mean)

#
combined_score <- matrix_of_iptm + matrix_of_iptm_cor*0.2

########
#benchmark against PPI data
########


#######
#read in the PPIs
all_interactions <- data.frame(readxl::read_excel("abb3758_table_s1.xlsx", sheet = 4))

#read in all of the sensed proteins
mpn_proteins1 <- data.frame(readxl::read_excel("abb3758_table_s1.xlsx", sheet = 2))[,c(1,3)]
mpn_proteins2 <- data.frame(readxl::read_excel("abb3758_table_s1.xlsx", sheet = 3))[,c(1,3)]
all_uniprots <- unique(c(unlist(mpn_proteins1),unlist(mpn_proteins2)))

#read in the translation from unpirot to mpn
uni_to_pmn <- read.delim("uniprotkb_taxonomy_id_272634_2025_04_07.tsv")
all_included_mpn <- rep(NA, length(all_uniprots))
for (i in seq_along(all_uniprots)){
  the_split <- strsplit(uni_to_pmn[match(all_uniprots[i],uni_to_pmn[,1]),5], split = " ")[[1]]
  of_int <- grep("MPN_", the_split)
  if (length(of_int) == 1){all_included_mpn[i] <- the_split[grep("MPN_", the_split)]}
}

#OMA translation from mgen to mpn
mgen_to_mpn2 <- read.delim("OMA-Pairs_MYCPN-MYCGE_2025-03-19_202844.tsv")[,c(16,8)]
mgen_to_mpn2[,1] <- paste("MG_" ,do.call("rbind", strsplit(mgen_to_mpn2[,1], split = "MG"))[,2], sep = "")

#run this code to get mgen to mpn rough translation
#source("~/Documents/AF3/mpn_to_mgen_mapping.R")
all_interactions$Gene_ID1 <- mgen_to_mpn2[match(all_interactions$Gene_ID1, mgen_to_mpn2[,2]),1]
all_interactions$Gene_ID2 <- mgen_to_mpn2[match(all_interactions$Gene_ID2, mgen_to_mpn2[,2]),1]
all_interactions <- all_interactions[which(!is.na(all_interactions$Gene_ID1) & !is.na(all_interactions$Gene_ID2)),]

#fix the "." genes
all_interactions$Gene_ID1[grep("\\.", all_interactions$Gene_ID1)] <- c("MG_480","MG_491","MG_480","MG_480","MG_515",
                                                                        "","","MG_480","MG_480","MG_478","MG_515", "MG_491", "MG_480","")
all_interactions$Gene_ID2[grep("\\.", all_interactions$Gene_ID2)] <- c("MG_480","MG_480","MG_515","MG_480","MG_522","MG_515")

#get matching iptm
matching_iptms <- rep(NA, 492)
for (i in seq_along(matching_iptms)){
  of_int1 <- grep(all_interactions$Gene_ID1[i], rownames(matrix_of_iptm))
  of_int2 <- grep(all_interactions$Gene_ID2[i], rownames(matrix_of_iptm))
  if (length(of_int1) == 1 & length(of_int2) == 1){
    matching_iptms[i] <- matrix_of_iptm[of_int1, of_int2]
  }
}

#
write.csv(cbind(all_interactions, matching_iptms),"poopsicle.csv")

temp_file <- cbind(all_interactions$Gene_ID1, all_interactions$Gene_ID2, all_interactions[,14:16])

#mgen not hit
mgen_expressed <- mgen_to_mpn2[match(all_included_mpn, mgen_to_mpn2[,2]),1]

#string mat
source("~/Documents/BMK analysis/make_matrix_of_string.R")
ppi_mat <- make_matrix_of_string(rownames(matrix_of_iptm), rownames(matrix_of_iptm), temp_file[,1:3])[[1]]
ppi_mat2 <- make_matrix_of_string(rownames(matrix_of_iptm), rownames(matrix_of_iptm), temp_file[,c(2,1,3)])[[1]]

ppi_mat2[which(is.na(ppi_mat2), arr.ind = TRUE)] <- 0
ppi_mat[which(is.na(ppi_mat), arr.ind = TRUE)] <- 0

ppi_mat <- ppi_mat +  ppi_mat2

#NA any mgens that are not in the mpn mass spec data (since they wouldnt )
not_expressed <- which(!(rownames(ppi_mat) %in% mgen_expressed))
#ppi_mat[not_expressed,not_expressed] <- NA

#some plots
source("~/Documents/Msmeg_Chem_CRISPRi/assess_string_v2.R")
ppi_assess <- assess_string_v2(matrix_of_iptm, ppi_mat, to_test = seq(0,0.7, by = 0.05), string_threshold = 0.1)
plot(ppi_assess[,c(4,3)], type = 'l', ylim = c(0,1))

#ROC plots (PPI from Mpn interactions)
plot(assess_string_v2(matrix_of_iptm, ppi_mat, to_test = seq(-0.1,0.7, by = 0.01), string_threshold = 0.1)[,c(6,3)], type = 'l')
abline(a=0,b=1,col = 'red')

pROC::auc(as.numeric(ppi_mat > 0.1),as.numeric(matrix_of_iptm))

#

########################################################
#write the lists and STRING to csv for Pedro
string_names <- apply(string_file[,1:2],1,paste0, collapse = "::::")
match(names(temp[[2]])[1:10], string_names)

#
poop <- cbind(string_file[match(names(temp[[2]]), string_names),], 
              write.list.to.matrix(temp[[2]]))
rownames(poop) <- names(temp[[2]])
write.csv(poop, "mgen_pairs_3.26.csv")


#
write.csv(matrix_of_iptm, "matrix_of_iptm_4.1.csv")
write.csv(matrix_of_iptm_cor, "matrix_of_iptm_cor.csv")

######################
#write the current list of hits
#with eggnog annotations
mgen_annot <- read.csv("mgen_annot2.csv")[,-1]

#
poop <- cbind(string_file[match(names(temp[[2]]), string_names),], 
              names(temp[[2]]), sapply(temp[[2]],mean))
poop <- poop[order(poop[,12], decreasing = TRUE)[1:1000],c(11,7,10,12)]

da_genes <- do.call("rbind", strsplit(poop[,1], split = "::::"))
poop <- cbind(poop, mgen_annot[match(da_genes[,1], mgen_annot[,2]), c(6,9,20),], 
                    mgen_annot[match(da_genes[,2], mgen_annot[,2]), c(6,9,20),])

write.csv(poop,"top1000_iptm.csv")

######################
#write the current list of hits
#with eggnog annotations
mgen_annot <- read.csv("mgen_annot2.csv")[,-1]

#
poop <- cbind(string_file[match(names(temp_ratio[[2]]), string_names),], 
              names(temp_ratio[[2]]), sapply(temp_ratio[[2]],mean))
poop <- poop[order(poop[,12], decreasing = TRUE)[1:1000],c(11,7,10,12)]

da_genes <- do.call("rbind", strsplit(poop[,1], split = "::::"))
poop <- cbind(poop, mgen_annot[match(da_genes[,1], mgen_annot[,2]), c(6,9,20),], 
              mgen_annot[match(da_genes[,2], mgen_annot[,2]), c(6,9,20),])

write.csv(poop,"top1000_iptm_ratio.csv")

#################
poop <- matrix_of_iptm_cor
poop[upper.tri(poop)] <- NA
of_int <- which(poop > 0.4, arr.ind = TRUE)

da_genes <- cbind(rownames(matrix_of_iptm)[of_int[,1]], rownames(matrix_of_iptm)[of_int[,2]])
da_genes_bind <- apply(da_genes,1,paste, collapse = "::::")

poop <- cbind(da_genes, string_file[match(da_genes_bind, string_names),c(7,10)], poop[of_int],
              mgen_annot[match(da_genes[,1], mgen_annot[,2]), c(6,9,20),],
              mgen_annot[match(da_genes[,2], mgen_annot[,2]), c(6,9,20),])

write.csv(poop,"top1000_iptm_cor.csv")



