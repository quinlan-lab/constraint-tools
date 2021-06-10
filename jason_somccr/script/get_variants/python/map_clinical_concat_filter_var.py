#!/usr/bin/env python3

import pandas as pd
import sys

## Define input variables
tcga_file = sys.argv[1]
icgc_file = sys.argv[2]
output_filename = sys.argv[3]

##################################################################################
##### Function to map TCGA and iCGC cancer codes to a consensus disease name #####
##################################################################################

 ## Read in the clinical information mapping file for icgc and tcga
file = "~/git/somccr/data/clinical/map_icgc_tcga_cancer_codes.xlsx"

## Only select for columns of interest
tcga = ['TCGA_Project_Code', 'Cancer_Type_Simplified']
icgc = ['ICGC_Project_Code', 'Cancer_Type_Simplified']
map_disease_tcga = pd.read_excel(file, 'master', usecols=tcga)
map_disease_icgc = pd.read_excel(file, 'master', usecols=icgc)

## Read in TCGA and ICGC data
tcga_df = pd.read_csv(tcga_file, sep="\t", compression="gzip")
icgc_df = pd.read_csv(icgc_file, sep="\t", compression="gzip")

## Map the project code to the cancer type name
tcga_df = pd.merge(tcga_df, map_disease_tcga, on="TCGA_Project_Code")
icgc_df = pd.merge(icgc_df, map_disease_icgc, on="ICGC_Project_Code")

## Change the column name of the column used to map project codes to TCGA and ICGC samples --> used to plot # of samples for each cancer type
tcga_df.rename(columns={'Tumor_Sample_Barcode':'Sample_ID'}, inplace=True)
icgc_df.rename(columns={'Donor_ID':'Sample_ID'}, inplace=True)

###################################################################################
##### Concatenate TCGA and ICGC variant data with mapped clinical information #####
###################################################################################

## Define columns of interest
fields = ['Chromosome', 'Start_Position', 'End_Position', 'Hugo_Symbol', 
          'Variant_Classification', 'Cancer_Type_Simplified', 'Sample_ID', 
          'Reference_Allele', 'Tumor_Seq_Allele2']

tcga_df = tcga_df[fields]
icgc_df = icgc_df[fields]

## Specify which study the variants came from 
tcga_df['study'] = "mc3_wes"
icgc_df['study'] = "pcawg_wgs"

## Concatenate mc3 and pcawg variants 
maf = pd.concat([tcga_df, icgc_df])

## Convert start_position to 0-based positions
maf.loc[:, 'Start_Position'] = maf['Start_Position'].apply(lambda x: x - 1)

###############################################
##### Select for specific variant classes #####
###############################################

## Remove non-cds/exonic variants 
variants_to_keep = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation',
                   'Nonsense_Mutation', 'Nonstop_Mutation', 'RNA', 'Silent', 'Splice_Site', 'Translation_Start_Site',
                   'Targeted_Region', 'Start_Codon_Del', 'Start_Codon_Ins', 'Start_Codon_SNP', 'Stop_Codon_Del', 
                   'Stop_Codon_Ins', 'De_novo_Start_InFrame', 'De_novo_Start_OutOfFrame']
keep_variants = maf.loc[maf["Variant_Classification"].isin(variants_to_keep)]

## Add "chr" to the chromosome name 
keep_variants['Chromosome'] = "chr" + keep_variants['Chromosome'].astype(str)

########################
##### Write output #####
########################
keep_variants.to_csv(output_filename, sep="\t", header=False, index=False)

