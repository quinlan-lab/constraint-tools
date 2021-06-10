##################################################################################
##### Function to determine cds/transcript length of genes from the MAF file #####
##################################################################################
## Load in necessary packages
library(data.table)
library(dplyr)
library(biomaRt)
library(pbmcapply)

## Define the function which uses biomaRt to get gene cds lengths
get_length <- function(input_filename, biomart_log_dir, output_filename) {
  ## Read in the variant data
  df <- fread(input_filename, header=FALSE, sep="\t") %>% `colnames<-`(c("chromosome", "start", "stop", "gene", "variant_class", "cancer_type", "sample_id", "ref", "alt", "working_group"))
  
  ## Define mart and ensembl database for biomaRt
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "http://dec2016.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
  #ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  ###########################################################################
  ##### For every gene in the df --> get the length (CDS or transcript) #####
  ###########################################################################
  ## Get the unique gene names in the variant dataset
  genes <- df$gene %>% unique
  
  ## Define the attributes (output of query)
  attributes <- c("external_gene_name", "ensembl_gene_id", "ensembl_transcript_id", 
                  "transcript_count", "cds_length", "transcript_length")
  
  ## Define the filters (inputs of query)
  filters <- "external_gene_name"
  
  ## Run the biomaRt query
  output <- getBM(attributes = attributes, filters = filters,
                  values = genes, mart = ensembl, useCache = FALSE)
  output <- data.table(output)

  ########################################## 
  ##### Get genes without biomaRt data #####
  ##########################################
  ## Initiate file to store genes without biomaRt data
  cat("", file = paste0(biomart_log_dir, "/genes_no_biomaRt_data.txt"), append = FALSE)
  warning_message <- paste0("The biomaRt query of some genes yielded no data... Skipping these genes but please refer to the file: ",
                            biomart_log_dir, "/genes_no_biomaRt_data.txt for these gene names...")
  no_biomart_genes <- genes[-which(genes %in% output$external_gene_name)]
  write(no_biomart_genes, 
        file = paste0(biomart_log_dir, "/genes_no_biomaRt_data.txt"), 
        append = TRUE)
  print(warning_message)
  
  #####################################################################################################################
  ##### Get the cds_length or transcript_length (includes CDS and 5'/3' UTRs) and other information for each gene #####
  #####################################################################################################################
  genes <- output$external_gene_name %>% unique
  gene_name <- genes[2]
  
  ## Initiate file to store genes without CDS length
  cat("", file = paste0(biomart_log_dir, "/genes_no_cds_length.txt"), append = FALSE)

  ## Iterate through each gene from the biomaRt query
  gene_length <- rbindlist(pbmclapply(mc.cores = 2, genes, function(gene_name, output){
    ## Subset biomaRt output for gene of interest
    gene_df <- output[external_gene_name == gene_name]
    
    ## Check to see if the gene has mutliple gene_ids and/or does not have a reported CDS length
    if(all(is.na(gene_df$cds_length)) | (length(unique(gene_df$ensembl_gene_id)) > 1)) {
      warning_message <- paste0("The gene: ", gene_name, " does not have a CDS length and/or has multiple gene ids... Skipping this gene...")
      print(warning_message)
      write(gene_name, file = paste0(biomart_log_dir, "/genes_no_cds_length.txt"), append = TRUE)
      return(data.table())
    }
    
    ## Get the maximum cds_length and other relevant information
    if (!all(is.na(gene_df$cds_length))) {
      cds_length <- max(na.omit(gene_df$cds_length))
      transcript_num <- gene_df$transcript_count[1]
      return(data.table("gene" = gene_name,
                        "length" = cds_length,
                        "transcript_num" = transcript_num))
    }
  }, output=output))
  
  #################################
  ##### Write the output file #####
  #################################
  write.table(x = gene_length, file = output_filename, quote = FALSE, sep = "\t", row.names = FALSE)
}

## Define arguments
#input_filename <- "/uufs/chpc.utah.edu/common/HIPAA/u1240855/git/somccr/data/output/somccr_v3/tcga_icgc.sorted.filtered.snv.CDS.bed"
#biomart_log_dir <- "/uufs/chpc.utah.edu/common/HIPAA/u1240855/git/somccr/data/output/biomart_log"
#output_filename <- "/uufs/chpc.utah.edu/common/HIPAA/u1240855/git/somccr/data/output/gene_cds_lengths.txt"
args <- commandArgs(trailingOnly = TRUE)
get_length(input_filename = args[1], biomart_log_dir = args[2], output_filename = args[3])





