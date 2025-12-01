#read
#list_of_proteins <- seqinr::read.fasta("subtilis_proteins.txt", forceDNAtolower = FALSE, as.string = TRUE, whole.header = TRUE)
#names(list_of_proteins) <- do.call("rbind", strsplit(do.call("rbind", strsplit(names(list_of_proteins), split = "locus_tag="))[,2], split = "] "))[,1]
#proteins_of_interest <- list_of_proteins[1:5]

#################
#make a list of list from a protein list
#################

make_list_of_af3_runs2 <- function(list_of_proteins, proteins_of_interest, max_run_size = 5000){
  
  #get rid of proteins that are just too long
  if (length(proteins_of_interest) > 1){
    max_protein_size <- max_run_size - sum(sapply(proteins_of_interest, nchar)) 
  } else {
    max_protein_size <- max_run_size
  }
  
  list_of_proteins <- list_of_proteins[which(sapply(list_of_proteins, nchar) <= max_protein_size)]
  
  #initialize
  list_of_af3_runs <- list()
  gene_included <- rep(FALSE, length(list_of_proteins))
  gene_touched <- rep(0, length(list_of_proteins))
  
  #while loop
  while(sum(gene_included) < length(list_of_proteins)){
    
    #initialize a temp group
    temp_group <- proteins_of_interest
    too_long <- rep(FALSE,10)
    available_genes <- which(gene_touched == 0)
    
    for (i in 1:10){
      #while the total length of the group is less than the max
      while (!too_long[i] & sum(gene_included) < length(list_of_proteins)){
        
        #pick a protein
        random_protein_index <- sample_real(available_genes,1)
        random_protein <- list_of_proteins[[random_protein_index]]
        
        #if that protein doesnt make the length go over max, add it to the list and set the things
        if (length(temp_group) > 0){
          current_group_length <- sum(sapply(temp_group, nchar)) + nchar(random_protein)
        } else {
          current_group_length <- nchar(random_protein)
        }
        
        if (current_group_length < max_run_size){
          gene_included[random_protein_index] <- TRUE
          gene_touched[random_protein_index] <- gene_touched[random_protein_index] + 1
          temp_group[[length(temp_group)+1]] <- random_protein
          names(temp_group)[length(temp_group)] <- names(list_of_proteins[random_protein_index])
        } else {
          too_long[i] <- TRUE
        }
      }
    }
    
    #save the temp group
    list_of_af3_runs[[length(list_of_af3_runs)+1]] <- temp_group
    
  }
  
  #return
  return(list_of_af3_runs)
  
}


#
sample_real <- function(x, size, replace = F, prob = NULL) {
  if (length(x) == 1) return(x)
  base::sample(x, size = size, replace = replace, prob = prob)
}


#######
