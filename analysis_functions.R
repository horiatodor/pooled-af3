
############
#set of functions to analyse AF3 genome-wide results
############

############
#function to get a list of iptm matrixes
############
#folder_level 1 for downloads, folder_level 2 for macosx compression

get_list_of_iptm <- function(afolder=NULL, azip=NULL, afolder_of_zip = NULL, 
                             description, description_delimiter = "::::", #folder_level = 2, 
                             return_paes = FALSE, return_sd = FALSE){
  
  #this function is set up to read in file in 3 different formats, all of which I'll probably use at some point
  list_of_iptm <- list()
  if (return_paes){list_of_pae <- list()}
  if (return_sd){list_of_sds <- list()}
  
  #1 (unzipped) folder of folds
  if (!is.null(afolder)){
    
    #get the files
    folders_of_interest <- list.dirs(afolder,recursive = FALSE)
    
    for (i in seq_along(folders_of_interest)){
      
      #get the files we care about
      all_files <- list.files(folders_of_interest[i])
      files_we_want <- all_files[grep("summary_confidences", all_files)]
      
      #read the iptms in
      list_of_runs <- list()
      if (return_paes){list_of_mins <- list()}
      
      for (j in seq_along(files_we_want)){
        list_of_runs[[j]] <- jsonlite::fromJSON(paste0(folders_of_interest[i], "/",files_we_want[j]))[[2]]
        if (return_paes){list_of_mins[[j]] <- jsonlite::fromJSON(paste0(folders_of_interest[i], "/",files_we_want[j]))[[3]]}
      }
      
      #average them, 
      iptms <- apply(simplify2array(list_of_runs), c(1, 2), mean)
      if (return_paes){paes <- apply(simplify2array(list_of_mins), c(1, 2), mean)}
      if (return_sd){sds <- apply(simplify2array(list_of_runs), c(1, 2), sd)}
        
      #get the gene names
      the_name <- strsplit(strsplit(files_we_want[1], split = "_summary_confidences_0.json")[[1]], split = "fold_")[[1]][[2]]
      
      #if its in the descriptions
      if (the_name %in% description[,1]){
        index <- which(description[,1] == the_name)
        gene_names <- unlist(strsplit(description[index,2], split = description_delimiter))
  
        #prepare the output
        rownames(iptms) <- gene_names
        colnames(iptms) <- gene_names
        
        if (return_paes){
          rownames(paes) <- gene_names
          colnames(paes) <- gene_names
        }
        
        if (return_sd){
          rownames(sds) <- gene_names
          colnames(sds) <- gene_names
        }
        
        #join the output
        list_of_iptm[[length(list_of_iptm)+1]] <- iptms
        names(list_of_iptm)[length(list_of_iptm)] <- the_name
        
        #paes
        if (return_paes){
          list_of_pae[[length(list_of_pae)+1]] <- paes
          names(list_of_pae)[length(list_of_pae)] <- the_name
        }
        
        #sds
        if (return_sd){
          list_of_sds[[length(list_of_sds)+1]] <- sds
          names(list_of_sds)[length(list_of_sds)] <- the_name
        }
        
      }
    }
  }
  
  #2 zipped folder of folds
  if (!is.null(azip)){
    
    #get all of the filenames (with summary_confidences)
    filenames <- unzip(azip, list=TRUE)[,1]
    macosx_stuff <- grep("__MACOSX", filenames)
    folder_level <- 1
    if (length(macosx_stuff) > 0){
      filenames <- filenames[-macosx_stuff]
      folder_level <- 2
    } 
    
    #folders in there
    folders <- do.call("rbind", strsplit(filenames,"/"))[,folder_level]
    unique_folders <- unique(folders)
    unique_folders <- unique_folders[which(unique_folders %in% description[,1])]
    
    #for each folder
    for (i in seq_along(unique_folders)){
      
      #get the files we care about
      files_we_want <- filenames[intersect(which(folders == unique_folders[i]),grep("summary_confidences", filenames))]
      
      #read the iptms in
      list_of_runs <- list()
      if (return_paes){list_of_mins <- list()}
      
      for (j in seq_along(files_we_want)){
        list_of_runs[[j]] <- jsonlite::fromJSON(unz(azip,files_we_want[j]))[[2]]
        if (return_paes){list_of_mins[[j]] <- jsonlite::fromJSON(unz(azip,files_we_want[j]))[[3]]}
      }
      
      #average them
      iptms <- apply(simplify2array(list_of_runs), c(1, 2), mean)
      if (return_paes){paes <- apply(simplify2array(list_of_mins), c(1, 2), mean)}
      if (return_sd){sds <- apply(simplify2array(list_of_runs), c(1, 2), sd)}
      
      index <- which(description[,1] == unique_folders[i])
      gene_names <- unlist(strsplit(description[index,2], split = description_delimiter))
      
      #prepare the output
      rownames(iptms) <- gene_names
      colnames(iptms) <- gene_names
      
      if (return_paes){
        rownames(paes) <- gene_names
        colnames(paes) <- gene_names
      }
      
      if (return_sd){
        rownames(sds) <- gene_names
        colnames(sds) <- gene_names
      }
      
      #join the output
      list_of_iptm[[length(list_of_iptm)+1]] <- iptms
      names(list_of_iptm)[length(list_of_iptm)] <- unique_folders[i]
      
      if (return_paes){
        list_of_pae[[length(list_of_pae)+1]] <- paes
        names(list_of_pae)[length(list_of_pae)] <- unique_folders[i]
      }
      
      if (return_sd){
        list_of_sds[[length(list_of_sds)+1]] <- sds
        names(list_of_sds)[length(list_of_sds)] <- unique_folders[i]
      }
    }

  }
  
  #3 zipped folder of folds
  if (!is.null(afolder_of_zip)){
    
    #initialize and recursively call self?
    zip_filesnames <- paste0(afolder_of_zip, list.files(afolder_of_zip))
    
    #recursively call the read zip?
    for (i in seq_along(zip_filesnames)){
      print(i)
      the_temp <- get_list_of_iptm(azip = zip_filesnames[i], description = description, #folder_level = folder_level, 
                                   return_paes=return_paes, return_sd = return_sd)
      
      if (return_paes){
        list_of_iptm[[i]]<- the_temp[[1]]
        list_of_pae[[i]]<- the_temp[[2]]
      } 
      
      if (return_sd){
        list_of_iptm[[i]]<- the_temp[[1]]
        list_of_sds[[i]]<- the_temp[[2]]
      }
      
      if (!return_paes & !return_sd){
        list_of_iptm[[i]] <- the_temp
      }
    }
    
    #join
    list_of_iptm <- do.call("c", list_of_iptm)
    if (return_paes){list_of_pae <- do.call("c", list_of_pae)}
    if (return_sd){list_of_sds <- do.call("c", list_of_sds)}
  }
  
  #return statement
  if (return_paes){return(list(list_of_iptm, list_of_pae))}
  if (return_sd){return(list(list_of_iptm, list_of_sds))}
  if (!return_paes & !return_sd){return(list_of_iptm)}
  
}

############
#function to size correct
############

size_correct <- function(list_of_iptm, list_of_proteins, coefs = NULL, root_to_use = 1/2, type_of_corr = "subtract", protein_sizes=NULL){
  
  #innitialize
  list_of_iptm_size_corr <- list_of_iptm
  
  #do the thing
  #initialize
  list_of_size_matrices <- list()
  if (is.null(protein_sizes)){
    protein_sizes <- sapply(list_of_proteins,nchar)
  }
  
  #for loop
  for (i in seq_along(list_of_iptm)){
    
    gene_lengths <- protein_sizes[match(rownames(list_of_iptm[[i]]),names(protein_sizes))]
    temp <- list_of_iptm[[i]]
    diag(temp) <- NA
    
    list_of_size_matrices[[i]] <- data.frame(rep(rownames(list_of_iptm[[i]]), each = length(gene_lengths)), 
                                             rep(rownames(list_of_iptm[[i]]), length(gene_lengths)),
                                             rep(gene_lengths, each = length(gene_lengths)), 
                                             rep(gene_lengths, length(gene_lengths)), 
                                             as.vector(temp),
                                             rep(sum(gene_lengths), length(as.vector(temp))),
                                             rep(length(gene_lengths), length(as.vector(temp))))
    
  }
  
  list_of_size_matrices <- do.call("rbind", list_of_size_matrices)
  list_of_size_matrices <- list_of_size_matrices[which(!is.na(list_of_size_matrices[,5])),]
  
  #get rid of duplicates
  list_of_size_matrices <- list_of_size_matrices[!duplicated(list_of_size_matrices[,1:2]),]
  
  #run the regression
  if (is.null(coefs)){
    coefs <- summary(robustbase::lmrob(list_of_size_matrices[,5] ~ as.numeric((list_of_size_matrices[,3] + list_of_size_matrices[,4])^(root_to_use))))$coef[,1]
    #coefs <- summary(robustbase::lmrob(list_of_size_matrices[,5] ~ as.numeric((list_of_size_matrices[,3] + list_of_size_matrices[,4])^(root_to_use)) + list_of_size_matrices[,7]))$coef[,1]
  }

  #do the correction
  #for loop
  for (i in seq_along(list_of_iptm)){
    
    gene_lengths <- protein_sizes[match(rownames(list_of_iptm[[i]]),names(protein_sizes))]
    mat1 <- matrix(rep(gene_lengths, times = length(gene_lengths)), ncol = length(gene_lengths))
    mat2 <- matrix(rep(gene_lengths, each = length(gene_lengths)), ncol = length(gene_lengths))
    
    if (type_of_corr == "subtract"){
      list_of_iptm_size_corr[[i]] <- list_of_iptm[[i]]-(coefs[1] + coefs[2]*((mat1+mat2)^(root_to_use)))
    }

    if (type_of_corr == "ratio"){
      list_of_iptm_size_corr[[i]] <- list_of_iptm[[i]]/(coefs[1] + coefs[2]*((mat1+mat2)^(root_to_use)))
    }
    
    if (type_of_corr == "variable"){
      coef1 <- -4.326e-02 + 5.753e-04*dim(list_of_iptm[[i]])[1]
      coef2 <-  4.810e-03 + -2.756e-05*dim(list_of_iptm[[i]])[1]
      list_of_iptm_size_corr[[i]] <- list_of_iptm[[i]]-(coef1 + coef2*((mat1+mat2)^(root_to_use)))
    }
  }
  
  #print
  print(coefs)
  print("normal int= -0.036, coef = 0.0044")
  
  #return
  return(list(list_of_iptm_size_corr, coefs, list_of_size_matrices))
  
}

############
#function to process into few by all
############
make_matrix_few_by_all <- function(list_of_iptm, constant_index = c(1,2,3,4,5)){
  
  #initialize the lists
  list_of_results <- list() 
  list_of_repeats <- list()
  
  #go through the list
  for (i in seq_along(list_of_iptm)){
    list_of_results[[i]] <- list_of_iptm[[i]][-constant_index,constant_index]
    list_of_repeats[[i]] <- list_of_iptm[[i]][constant_index,constant_index]
  }
  
  #make and prettify a matrix
  temp <- apply(simplify2array(list_of_repeats), c(1, 2), mean)
  diag(temp) <- NA
  matrix_of_results <- rbind(temp, do.call("rbind", list_of_results))
  
  #return
  return(matrix_of_results)
  
}


############
#function to process into all by all
############
make_matrix_all_by_all <- function(list_of_iptm, list_of_proteins, resolve_reps = "max"){
  
  #initialize
  matrix_of_iptm <- matrix(NA, length(list_of_proteins), length(list_of_proteins))
  rownames(matrix_of_iptm) <- names(list_of_proteins)
  colnames(matrix_of_iptm) <- names(list_of_proteins)
  matrix_of_replicate_number <- matrix_of_iptm
  
  #######################
  #get out all of the individual pairs
  list_of_pairs <- list()
  
  #
  for (i in seq_along(list_of_iptm)){
    
    temp <- list_of_iptm[[i]]
    
    #just upper tri
    mat1 <- matrix(rep(rownames(temp), each = dim(temp)[1]), ncol = dim(temp)[1])
    mat2 <- matrix(rep(rownames(temp), times = dim(temp)[1]), ncol = dim(temp)[1])
    mat3 <- matrix(rep(names(list_of_iptm)[i], dim(temp)[1]*dim(temp)[1]), ncol = dim(temp)[1])
    mat4 <- matrix(rep(diag(temp), each = dim(temp)[1]), ncol = dim(temp)[1]) +
      matrix(rep(diag(temp), times = dim(temp)[1]), ncol = dim(temp)[1])
    
    list_of_pairs[[i]] <- data.frame(mat1[upper.tri(mat1)], 
                                     mat2[upper.tri(mat2)],
                                     mat3[upper.tri(mat3)],
                                     temp[upper.tri(temp)],
                                     mat4[upper.tri(mat4)])
    
  }
  
  matrix_of_pairs <- do.call("rbind", list_of_pairs)
  matrix_of_pairs <- matrix_of_pairs[which(!is.na(matrix_of_pairs[,4])),]
  
  #######################
  #put them into the matrix and 
  #list of replicates
  list_of_pairs_in_set <- vector("list", length(list_of_proteins)^2) 
  matching_list_of_runs <- vector("list", length(list_of_proteins)^2) 
  the_names <- rep("", length(list_of_proteins)^2)
  
  #for loop
  for (i in 1:(length(list_of_proteins)-1)){
    
    of_int <- which(matrix_of_pairs[,1] %in% names(list_of_proteins)[i] |
                      matrix_of_pairs[,2] %in% names(list_of_proteins)[i])
    one_pairs <- matrix_of_pairs[of_int, ]
    
    for (j in (i+1):length(list_of_proteins)){
      
      #what do we have
      of_int2 <- which(one_pairs[,1] %in% names(list_of_proteins)[j] |
                         one_pairs[,2] %in% names(list_of_proteins)[j])
      
      #if we have one value put it in and put the number into the other matrix
      if (length(of_int2) == 1){
        
        #iptm matrix
        matrix_of_iptm[i,j] <- one_pairs[of_int2,4]
        matrix_of_iptm[j,i] <- one_pairs[of_int2,4]
        
        #replicates
        matrix_of_replicate_number[i,j] <- 1
        matrix_of_replicate_number[j,i] <- 1
        
      }
      
      #if we have more than one value, 
      if (length(of_int2) > 1){
        
        #put in the average in the matrix
        if (resolve_reps == "max"){
          matrix_of_iptm[i,j] <- max(one_pairs[of_int2,4])
          matrix_of_iptm[j,i] <- max(one_pairs[of_int2,4])
        }
        
        if (resolve_reps == "mean"){
          matrix_of_iptm[i,j] <- mean(one_pairs[of_int2,4])
          matrix_of_iptm[j,i] <- mean(one_pairs[of_int2,4])
        }
        
        if (resolve_reps == "min"){
          matrix_of_iptm[i,j] <- min(one_pairs[of_int2,4])
          matrix_of_iptm[j,i] <- min(one_pairs[of_int2,4])
        }
        
        if (resolve_reps == "median"){
          matrix_of_iptm[i,j] <- median(one_pairs[of_int2,4])
          matrix_of_iptm[j,i] <- median(one_pairs[of_int2,4])
        }
        
        if (resolve_reps == "sum_ptm"){
          matrix_of_iptm[i,j] <- one_pairs[of_int2[which.max(one_pairs[of_int2,5])],4]
          matrix_of_iptm[j,i] <- one_pairs[of_int2[which.max(one_pairs[of_int2,5])],4]
        }
        
        if (resolve_reps == "geo_mean"){
          matrix_of_iptm[i,j] <- exp(mean(log(one_pairs[of_int2,4])))
          matrix_of_iptm[j,i] <- exp(mean(log(one_pairs[of_int2,4])))
        }
        
        if (resolve_reps == "special"){
          poop <- ifelse(min(one_pairs[of_int2,4]) < 0, min(one_pairs[of_int2,4]), mean(one_pairs[of_int2,4]))
          matrix_of_iptm[i,j] <- poop
          matrix_of_iptm[j,i] <- poop
        }
        
        if (resolve_reps == "first"){
          matrix_of_iptm[i,j] <- one_pairs[of_int2[1],4]
          matrix_of_iptm[j,i] <- one_pairs[of_int2[1],4]
        }
        
        #put the number in the matrix
        matrix_of_replicate_number[i,j] <- length(of_int2)
        matrix_of_replicate_number[j,i] <- length(of_int2)
        
      }
      
      #list of pairs
      if (length(of_int2) > 0){
        
        the_index <- ((i-1)*length(list_of_proteins))+j
        list_of_pairs_in_set[[the_index]] <- one_pairs[of_int2,4]
        matching_list_of_runs[[the_index]] <- one_pairs[of_int2,3]
        the_names[the_index] <- paste0(names(list_of_proteins)[i],"::::",names(list_of_proteins)[j])
        
      }
    }
  }
  
  #do the names
  names(list_of_pairs_in_set) <- the_names
  names(matching_list_of_runs) <- the_names
  
  #remove empties
  not_empties <- which(sapply(list_of_pairs_in_set, length) > 0)
  list_of_pairs_in_set <- list_of_pairs_in_set[not_empties]
  matching_list_of_runs <- matching_list_of_runs[not_empties]
  
  #print some summary stats
  print(paste0("number_of_runs = ", length(list_of_iptm)))
  print(paste0("number_of_pairs = ", dim(matrix_of_pairs)[1]))
  
  #return
  return(list(matrix_of_iptm, list_of_pairs_in_set, matrix_of_replicate_number, matching_list_of_runs, matrix_of_pairs))
  
}

############
#function to assess replicates from all by all
############

expand_list_of_pairs <- function(list_of_pairs_in_set, return_list = FALSE){
  
  #for each pair that has more than one
  multiple_ind <- which(sapply(list_of_pairs_in_set, length) > 1)
  replicate_matrix_list <- list()
  
  for (i in seq_along(multiple_ind)){
    
    number_of_reps <- length(list_of_pairs_in_set[[multiple_ind[i]]])
    
    poop <- matrix(1:(number_of_reps^2),ncol=number_of_reps)
    replicate_matrix_list[[i]] <- expand.grid(list_of_pairs_in_set[[multiple_ind[i]]], 
                                              list_of_pairs_in_set[[multiple_ind[i]]])[poop[upper.tri(poop)],, drop = FALSE]
    

  }

  names(replicate_matrix_list) <- names(list_of_pairs_in_set)[multiple_ind]
  
  if (return_list){return(replicate_matrix_list)}  
  return(do.call("rbind", replicate_matrix_list))
  
}


############
#write list to matrix
############

#
write.list.to.matrix <- function(alist, blank_to_use = ""){
  list_lengths <- sapply(alist, length)
  amatrix <- matrix(blank_to_use, length(alist), max(list_lengths, na.rm = TRUE))
  rownames(amatrix) <- names(alist)
  for (i in seq_along(alist)){
    amatrix[i,1:list_lengths[i]] <- alist[[i]]
  }
  return(amatrix)
}

############
#function to resolve multiple runs of the same protein set and provide their replicability 
############

resolve_duplicates <- function(list_of_iptm){
  
  #initialize
  list_averages <- list()
  list_of_dupes <- list()
  unique_names <- unique(names(list_of_iptm))
  
  #if there are no duplicates, just quit
  if (length(unique_names) == length(list_of_iptm)){
    return(list(list_of_iptm, list_of_dupes))
  }

  #for loop
  for (i in seq_along(unique_names)){
    
    #which one is unique
    of_int <- which(names(list_of_iptm) == unique_names[i])
    
    #
    if (length(of_int) == 1){
      #if 1 just put it in
      list_averages[[i]] <- list_of_iptm[[of_int]]
    } else {
      list_averages[[i]] <- apply(simplify2array(list_of_iptm[of_int]), c(1, 2), mean)
      list_of_dupes[[length(list_of_dupes)+1]] <- list_of_iptm[of_int]
    }
  }
    
  names(list_averages) <- unique_names
  
  #return averaged iptm list and list of duplicates
  return(list(list_averages, list_of_dupes))
  
}

#######
#generate unique triples
#######

get_triples <- function(amatrix){
  
  #sort the names so they'll always be in the same order
  the_order <- order(rownames(amatrix))
  amatrix <- amatrix[the_order,the_order]
  
  #initialize
  v1 <- 1:length(the_order)
  combinations_to_use <- lapply(seq_along(v1), function(i) combn(v1, i, FUN = list))[[3]]
  results <- matrix(NA, length(combinations_to_use), 2)
  rownames(results) <- rep("", length(combinations_to_use))
    
  #nested for loops
  for (i in seq_along(combinations_to_use)){
    poop <- expand.grid(combinations_to_use[[i]], combinations_to_use[[i]])
    poop <- poop[which(poop[,1] < poop[,2]),]
    results[i,1] <- min(amatrix[as.matrix(poop)])
    rownames(results)[i] <- paste0(rownames(amatrix)[combinations_to_use[[i]]], collapse = "::::")  
  }
  
  #return
  return(results)
  
}