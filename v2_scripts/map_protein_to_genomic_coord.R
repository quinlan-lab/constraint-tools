library(data.table)
library(dplyr)
library(ensembldb)
library(AnnotationHub)
library(EnsDb.Hsapiens.v86)
library(biomaRt)

## NOTE: Using ensembl v86 to map protein to genomic coordinates
## Ensembl archive: http://oct2016.archive.ensembl.org/index.html
## Mapping strategy with ensembldb: https://www.bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/coordinate-mapping.html

#######################################
##### Define the ensembl database #####
#######################################
edb <- EnsDb.Hsapiens.v86

###############################################
##### Define the gene and interval region #####
###############################################
## Using cBioPortal example for HSP90: https://www.cbioportal.org/results/mutations?tab_index=tab_visualize&Action=Submit&session_id=60f64961e4b015b63e9f61f5&plots_horz_selection=%7B%7D&plots_vert_selection=%7B%7D&plots_coloring_selection=%7B%7D

gene <- "HSP90AA1"
p_start <- 253
p_end <- 268
transcript <- "ENST00000216281"

##########################################################
##### Use biomaRt to map transcript id to protein id #####
##########################################################

## Define the biomaRt dataset
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "oct2016.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")

## Define the attributes
attributes <- c("ensembl_peptide_id")

## Define the filters
filters <- "ensembl_transcript_id"

## Run the biomaRt query
p_id <- getBM(attributes = attributes, filters = filters, values = transcript, mart = mart, useCache = FALSE)$ensembl_peptide_id

#########################################################################
##### Map protein coordinates to genomic coordinates with ensembldb #####
#########################################################################

## Define the IRanges object
p_ranges <- IRanges(start = p_start, end = p_end, names = p_id)

## Get the genomic coordinates
g_coord <- proteinToGenome(p_ranges, edb)

## Extract genomic coordinate information
chr <- paste0("chr", g_coord[[1]]@seqnames %>% as.character())
g_start <- g_coord[[1]]@ranges@start
g_end <- g_coord[[1]]@ranges@start + g_coord[[1]]@ranges@width - 1

print(paste0("Genomic coordinates for ", gene, " are ", chr, ":", g_start, "-", g_end))
