#######################################################################
##### Function to make mutation lolliplots with additional tracks #####
#######################################################################

## Load in necessary packages
library(data.table)
library(dplyr)
library(Gviz)
library(rtracklayer)
library(trackViewer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(biomaRt)

## Use this link: https://bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/trackViewer.html#Plot_the_lollipop_plot_with_the_coverage_and_annotation_tracks

## Read in mutation data
df <- fread("~/Google Drive/Quinlan Lab - PhD/Projects/somatic_ccr/data/output/mc3_pcawg_v3/mc3_pcawg.sorted.filtered.CDS.bed")

gene <- "PIK3CA"
mut_df <- df[gene == "PIK3CA" & cancer_type == "Breast invasive carcinoma"]

## Steps: 
# 1) Use biomaRt to get coding sequence exons (defined as exons with known CDS start/end in biomaRt)
# 2) Plot mutations whose start positions lies within a an exonic coding sequence 

## Define the mart --> Use GRCh37 
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)

## Get gene information
gene_info <- getBM(attributes = c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "rank", "start_position", "end_position",
                                  "cds_start", "cds_end", "ensembl_transcript_id", "cds_length"),
                   filters = "external_gene_name", values = gene, mart = ensembl) %>% data.table %>% na.omit

## Use transcript with the longest coding sequence length
gene_df <- gene_info[which(cds_length == max(as.numeric(gene_info$cds_length)))]

## Setup transcript biomaRt information for bedtools intersect 
gene_df$chromosome_name <- paste0("chr", gene_df$chromosome_name)
gene_df$exon_chrom_start <- gene_df$exon_chrom_start - 1

## Write out gene mutation and biomart data for bedtools intersect
setwd("~/Google Drive/Quinlan Lab - PhD/Projects/somatic_ccr/data/output/mc3_pcawg_v3/intersect/")
output_file <- paste0(gene, "_biomart.bed")
write.table(x = gene_df, file = output_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
mut_file <- paste0(gene, "_mut.bed")
write.table(x = mut_df, file = mut_file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

## Run bedtools bash command to remove BED variants whose start/stop coordinates lie outside coding exons
cmd <- "bedtools intersect -a *_mut.bed -b *_biomart.bed > cds_mut.bed"
system(cmd)
cds_mut_df <- fread("cds_mut.bed") %>% `colnames<-`(colnames(mut_df))


x








names(id) <- gene
genes <- geneTrack(ids = as.character(id), txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, symbols = gene, asList=FALSE)

id <- getBM(attributes = c("entrezgene_id", "cds_start", "cds_end", "ensembl_exon_id"), filters = "external_gene_name", values = gene, mart = ensembl)

gene <- geneTrack(get(gene, org.Hs.egSYMBOL2EG), TxDb.Hsapiens.UCSC.hg19.knownGene)[[1]]


chr <- df$chromosome %>% unique()
positions <- df$start
data.table(table(positions))
mut <- GRanges(chr, IRanges(start = unique(positions), width = 1))
mut$score <- data.table(table(positions))$N


positions <- df$start
data.table(table(positions))
mut <- GRanges(chr, IRanges(start = unique(positions), width = 1))
mut$score <- data.table(table(positions))$N

