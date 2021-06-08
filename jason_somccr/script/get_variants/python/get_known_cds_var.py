#!/usr/bin/env python3

import pandas as pd
import sys

## Define input variables
input_filename = sys.argv[1]
known_cds_genes_file = sys.argv[2]
output_filename = sys.argv[3]

## Read in the sorted and filtered variant file
df = pd.read_csv(input_filename, sep="\t", header=None)
df.columns = ['chromosome', 'start', 'stop', 'gene', 'variant_class', 'cancer_type', 'sample_id', 'ref', 'alt', 'working_group']

############################################################
##### Get variants within genes with known CDS lengths #####
############################################################

## Get genes with known CDS length
cds_gene_df = pd.read_csv(known_cds_genes_file, sep="\t")

## Get unique list of CDS genes
cds_genes = cds_gene_df['gene'].tolist()
cds_genes = list(set(cds_genes))

## Get rows from the combined MAF dataset with cds genes
cds_maf = df[df["gene"].isin(cds_genes)]

#################################
##### Write the output file #####
#################################
cds_maf.to_csv(output_filename, sep="\t", index=False)