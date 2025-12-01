
require(jsonlite)
#jsonlite guide https://arxiv.org/pdf/1403.2805

#' 
#' each .json file contains [batch size] of jobs - right now set to 30, since that's the daily use per day.
#' first, setup_run calls Horia's make_list_of_af3.R script.
#' This produces a huge list of lists of proteins for every job. From there, they are broken up into batches and written:
#' 1. setup_run calls job_parser, which breaks up the list into smaller lists of [batch size] 
#' 2. job_parser makes an information file (which proteins are in which job) for each batched list of jobs, 
#'    and calls job_writer to actually write the .json format for each job. It takes all the .json formats for each 
#'    batch and puts them into one .json file, so you can drag and drop one file instead of 30 into the web server 
#'    (google lets you put up to 100 jobs in a .json, but 30 is easier for bookkeeping)
#' 4. finally, setup_run calls write_files to write everything to the specified directory.
#' 


#functions ################

#takes in all the info of few and all proteins ('tokens' - eventually we will generalize this to include RNA), call's Horia's code to
#make all the runs, then breaks them up into .json files. Finally, writes the .json files and informational .csv files to directory
setup_run <- function(all_tokens, tokens_of_int, name_prefix = "alphafold_runs", output_dir = "", batch_size = 30) {
  #make the list of runs
  runs <- make_list_of_af3_runs(all_tokens, tokens_of_int)

  #call to get all the runs  written into json files
  parsed_jobs <- job_parser(runs, name_prefix, batch_size = 30) #batch size indicates how many runs to put in each json
  
  #write all batches to separate files - and information files after. 
  write_files(parsed_jobs, name_prefix, output_dir)
}

#this writes all the parsed jobs into whatever output directory you give it. For now, it dumps both the .json files and informational
# .csv files into the same folder, but it would be easy enough to change this.
write_files <- function(parsed_jobs, name_prefix, output_dir) {
  for(i in seq_along(parsed_jobs)) {
    file_name <- paste(name_prefix,i,".json",sep="")
    write(jsonlite::toJSON(parsed_jobs[[i]][[2]],auto_unbox = TRUE, pretty=FALSE),paste(output_dir,file_name,sep=""))
  }
  for(i in seq_along(parsed_jobs)) {
    file_name <- paste(name_prefix,i,".csv",sep="")
    write.csv(parsed_jobs[[i]][[1]],paste(output_dir,file_name,sep=""))
  }
  print("done writing")
}


#returns list of lists. Each outer list is up to [batch size] jobs, containing the table of info and the list which can be converted into a json file
#name prefix is passed to functions for naming each job and bookkeeping in the information file
#this function doles out the lists in the specified batch size
job_parser <- function(list_of_af3s, name_prefix ="", batch_size = 30) {
  
  batch_list <- list()
  num_batches <- ceiling(length(list_of_af3s) / batch_size) #only accepts up to 100 jobs at a time
  #print(num_batches)
  start <- 1
  
  for(i in seq_along(1:num_batches)) {
    end <- i * batch_size
    #print(start)
    #print(end)
    if (end > length(list_of_af3s)) {
      end <- length(list_of_af3s)
    }
    batch_list <- append(batch_list, list(job_filer(list_of_af3s[start:end], (start - 1), name_prefix)))
    start <- (i) * batch_size + 1
    #print("done batching")
     
  }
 return(batch_list) 
}

#given set of protein lists, make 1. an information file of each job including the name, proteins, and number of proteins and 2. a list to be converted to json
job_filer <- function(list_of_af3s,num_to_start = 1, name_to_append = "") {
  
  num_jobs <- length(list_of_af3s)
  info_table <- matrix(NA, nrow = num_jobs, ncol = 3)
  colnames(info_table) <- c("job_number","proteins_tested","num_proteins_tested")
  jsonlist <- list()

  for(i in seq_along(list_of_af3s)) {
    #delimiter is :::: to try and avoid anything that could possibly be used in annotation
    info_table[i,] <- c(paste(name_to_append,num_to_start + i,sep="_"), paste(names(list_of_af3s[[i]]),collapse="::::"),length(list_of_af3s[[i]]))
    
    jsonlist <- append(jsonlist,job_writer(list_of_af3s[[i]],paste(name_to_append,num_to_start + i,sep="_")))
    
  }
 
  return(list(info_table,jsonlist)) 
}

#given single list of proteins, make a job 
job_writer <- function(set_of_proteins,job_num) {
  
  #make list of all protein chains and sequences for formatting
  sequencesList <- list()
  for(j in seq_along(set_of_proteins)) {
    sequencesList = append(sequencesList,list(list(proteinChain = list(sequence=set_of_proteins[[j]][1],count=1))))
  }
  
  #add random extra info alphafold needs in the .json
  currlist <- list(list(
    name=toString(job_num),
    modelSeeds=list(),
    sequences=(sequencesList),
    dialect="alphafoldserver",
    version=1
  ))
  
  jobjson <- (currlist)
  return(jobjson)
  
}
