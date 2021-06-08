#!usr/bin/env python

import pandas as pd
import sys

#####################################
##### Read in the TCGA MAF file #####
#####################################

## Define dtypes
# tcga_dtypes = pd.read_csv("https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc", 
#                      sep="\t", compression="gzip", header=0, nrows=10)
# dtypes_dict = tcga_dtypes.dtypes.apply(lambda x: x.name).to_dict()

## Specify dtypes and read in TCGA MAF file
print("Reading in TCGA MAF file")
mc3_df = pd.read_csv("https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc", 
                     sep="\t", compression="gzip", header=0)


#####################################################
##### Map Clinical information to TCGA variants #####
#####################################################

## Define columns to read in
fields = ['case_submitter_id', 'project_id']

## Read in the MC3 clinical data
clinical_mc3 = pd.read_csv("~/git/somccr/data/clinical/mc3_mapping.tsv", sep="\t", usecols=fields)

## Obtain dataset with the TCGA barcode and cancer type
clinical_mc3 = clinical_mc3.drop_duplicates()
clinical_mc3.columns = ['Tumor_Sample_Barcode_split', 'TCGA_Project_Code']

## Split the Tumor_Sample_Barcode column so that they match the barcodes in the 'clinical' dataset
barcodes = mc3_df['Tumor_Sample_Barcode']
barcodes_split = barcodes.str.rsplit(pat="-", n=4, expand=True)[0]

## Add the split TCGA barcodes to the mc3 variant dataset
mc3_df['Tumor_Sample_Barcode_split'] = barcodes_split

## Map the barcodes to the cancer type 
mc3_df = pd.merge(mc3_df, clinical_mc3, on = "Tumor_Sample_Barcode_split")

###########################
##### Output the data #####
###########################

## Define the output file 
output_filename = sys.argv[1]

## Write the output file
print("Writing TCGA raw variant file")
mc3_df.to_csv(output_filename, sep="\t", header=True, index=False, compression="gzip")
