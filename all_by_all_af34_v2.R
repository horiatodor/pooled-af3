#
source("~/Documents/AF3/make_list_of_af3.R")

#read
#list_of_proteins <- seqinr::read.fasta("mgen_proteins.txt", forceDNAtolower = FALSE, as.string = TRUE, whole.header = TRUE)
#names(list_of_proteins) <- do.call("rbind", strsplit(do.call("rbind", strsplit(names(list_of_proteins), split = "locus_tag="))[,2], split = "] "))[,1]

#################
#make all by all lists
#################

all_by_all <- function(list_of_proteins, max_protein_size = 2500){
  
  #remove proteins that are too long
  prot_length <- sapply(list_of_proteins, nchar)
  list_of_proteins <- list_of_proteins[which(prot_length <= max_protein_size)]

  #shuffle the list of proteins
  list_of_proteins <- list_of_proteins[sample(length(list_of_proteins))]
  
  #initialize results output and matrix to keep track kof pairs
  list_of_results <- list()
  matrix_of_pairs <- matrix (0L, length(list_of_proteins), length(list_of_proteins))
  rownames(matrix_of_pairs) <- names(list_of_proteins)
  colnames(matrix_of_pairs) <- names(list_of_proteins)
  diag(matrix_of_pairs) <- 1L
  
  #for loop
  while (sum(matrix_of_pairs == 0) > 0){
    
    #pick a protein that
    blanks_per_gene <- apply(matrix_of_pairs == 0,1,sum)
    proteins_of_int <- sample_real(which(blanks_per_gene == max(blanks_per_gene)),1)
    
    #make a list
    poop <- make_one_list(list_of_proteins, proteins_of_int, max_run_size = 5000, matrix_of_pairs)
    list_of_results[[length(list_of_results)+1]] <- poop[[1]]
    matrix_of_pairs <- make_redundancy_matrix2(list_of_results, list_of_proteins)
    
    #where we are
    if (length(list_of_results) %% 100 == 0){print(length(list_of_results))}
  }
  
  #the return statement
  return(list(list_of_results, matrix_of_pairs))
  
} 

#
#sum(sapply(list_of_results, length))
#hist(matrix_of_pairs, 100)


#################
#code to assess the redundancy of a set - or quantify the number of new interactiosn
#################
#matrix_of_current_things <- matrix_of_pairs
#list_of_possible_splits <- the_lists
assess_redundancy <- function(list_of_possible_splits, matrix_of_current_things, list_of_proteins){
  
  new_interactions <- rep(NA, length(list_of_possible_splits))
  
  #for each thing in the 
  for (k in seq_along(list_of_possible_splits)){
    
    #get all of the pairs
    pairs_in_set <- expand.grid(names(list_of_possible_splits[[k]]), names(list_of_possible_splits[[k]]))
    pairs_in_set <- cbind(match(pairs_in_set[,1], names(list_of_proteins)), match(pairs_in_set[,2], names(list_of_proteins)))
    pairs_in_set <- pairs_in_set[which(pairs_in_set[,1] != pairs_in_set[,2]),]
    
    #plug into the matrix
    new_interactions[k] <- length(which(matrix_of_current_things[pairs_in_set] == 0))
    
  }
  
  #return
  return(sum(new_interactions)/length(list_of_possible_splits))
  
}

#
#list_of_jobs <- all_by_all(list_of_proteins)

#################
#code to assess the redundancy of a set - or quantify the number of new interactiosn
#################
make_redundancy_matrix2 <- function(a_run_list, list_of_proteins){
  
  #initialize results output and matrix to keep track kof pairs
  matrix_of_pairs <- matrix (0L, length(list_of_proteins), length(list_of_proteins))
  rownames(matrix_of_pairs) <- names(list_of_proteins)
  colnames(matrix_of_pairs) <- names(list_of_proteins)
  diag(matrix_of_pairs) <- 1L
  
  #get all pairs from within the lists
  for (k in seq_along(a_run_list)){
    
    #get all of the pairs
    matched_names <- match(names(a_run_list[[k]]), names(list_of_proteins))
    pairs_in_set <- cbind(rep(matched_names, times = length(matched_names)),
                          rep(matched_names, each = length(matched_names)))
    pairs_in_set <- pairs_in_set[which(pairs_in_set[,1] != pairs_in_set[,2]),]
    
    #plug into the matrix
    matrix_of_pairs[pairs_in_set] <- matrix_of_pairs[pairs_in_set] + 1L
    
  }
  
  return(matrix_of_pairs)
  
} 

######################
#first lets get the number of tokens in each maboober
get_number_of_aa <- function(listicle){
  number_of_aa <- rep(NA, length(listicle))
  for (i in seq_along(number_of_aa)){number_of_aa[i] <- sum(sapply(listicle[[i]], nchar))}
  return(number_of_aa)
}

######
#
update_matrix <- function(new_list, matrix_of_pairs, list_of_proteins){
  
  #get all pairs from within the lists
  for (k in seq_along(new_list)){
    
    #get all of the pairs
    if (length(new_list[[k]]))
      pairs_in_set <- expand.grid(names(new_list[[k]]), names(new_list[[k]]))
    pairs_in_set <- cbind(match(pairs_in_set[,1], names(list_of_proteins)), match(pairs_in_set[,2], names(list_of_proteins)))
    pairs_in_set <- pairs_in_set[which(pairs_in_set[,1] != pairs_in_set[,2]),]
    
    #plug into the matrix
    matrix_of_pairs[pairs_in_set] <- matrix_of_pairs[pairs_in_set] + 1
  }
  
  return(matrix_of_pairs)
  
}

#
make_one_list <- function(list_of_proteins, first_protein_index, max_run_size = 5000, matrix_of_pairs){
  
  #initialize a temp group
  temp_group <- list_of_proteins[[first_protein_index]]
  names(temp_group)[1] <- names(list_of_proteins)[first_protein_index]
  too_long <- rep(FALSE,10)
  group_index <- first_protein_index

  for (i in 1:10){
    
  #while the total length of the group is less than the max
  while (!too_long[i]){
    
    #for each gene that is not calculate how many new interactions adding each gene 
    available_genes <- which(!names(list_of_proteins) %in% names(temp_group))
    
    number_added <- apply(matrix_of_pairs[group_index,available_genes, drop = FALSE] == 0,2,sum)
    possible_genes <- available_genes[which(number_added == max(number_added))]
    
    #if there's only one possible gene, add a few options
    #if (length(possible_genes) == 1 & sum(sapply(temp_group, nchar)) > 4500){
     # possible_genes <- available_genes[order(number_added, decreasing = TRUE)[1:3]]}
    
    random_protein_index <- sample_real(possible_genes,1)
    random_protein <- list_of_proteins[[random_protein_index]]
    
    #if that protein doesnt make the length go over max, add it to the list and set the things
    current_group_length <- sum(sapply(temp_group, nchar)) + nchar(random_protein)
    
    if (current_group_length < max_run_size){
      
      #addit to the group
      group_index <- c(group_index, random_protein_index)
      temp_group[[length(temp_group)+1]] <- random_protein
      names(temp_group)[length(temp_group)] <- names(list_of_proteins[random_protein_index])
      
      #update the matrix
      matrix_of_pairs <- update_matrix(list(temp_group), matrix_of_pairs, list_of_proteins) 
      
    } else {
      too_long[i] <- TRUE
    }
  }
  }
  
  return(list(temp_group, matrix_of_pairs))
  
}
