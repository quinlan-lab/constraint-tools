#!/usr/bin/env python

import pandas as pd
import sys

#####################################
##### Read in the ICGC MAF file #####
#####################################

## Define dtypes
# icgc_dtypes = pd.read_csv("https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz",
#                       compression='gzip', sep="\t", nrows=10)
# dtypes_dict = icgc_dtypes.dtypes.apply(lambda x: x.name).to_dict()

## Specify dtypes and read in ICGC MAF file
print("Reading in ICGC MAF file")
pcawg_df = pd.read_csv("https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz",
                      compression='gzip', sep="\t", header=0)

## Rename columns Start/End_position to Start/End_Position
pcawg_df.rename(columns={'Start_position':'Start_Position', 'End_position':'End_Position'}, inplace=True)

#####################################################
##### Map clinical information to ICGC varaints #####
#####################################################

## Define columns to read in
fields = ['project_code', 'icgc_donor_id']

## Read in the MC3 clinical data
clinical_pcawg = pd.read_csv("~/git/somccr/data/clinical/pcawg_mapping.tsv", 
                             sep="\t", usecols=fields)

## Obtain dataset with the TCGA barcode and cancer type
clinical_pcawg = clinical_pcawg.drop_duplicates()
clinical_pcawg.columns = ['ICGC_Project_Code', 'Donor_ID']

# ## Remove the country designation in the cancer type code
# cancer_type = clinical_pcawg['cancer_type_pcawg']
# cancer_type_split = cancer_type.str.rsplit(pat="-", n=2, expand=True)[0]

# ## Add the split ICGC barcodes to the pcawg clinical dataset
# clinical_pcawg['cancer_type_pcawg'] = cancer_type_split

## Map the barcodes to the cancer type 
pcawg_df = pd.merge(pcawg_df, clinical_pcawg, on = "Donor_ID")

###########################
##### Output the data #####
###########################

## Define the output file 
output_filename = sys.argv[1]

## Write the output file
print("Writing ICGC raw varaint file")
pcawg_df.to_csv(output_filename, sep="\t", header=True, index=False, compression="gzip")


