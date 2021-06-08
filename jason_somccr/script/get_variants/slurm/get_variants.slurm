#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --account=ucgd-rw
#SBATCH --partition=ucgd-rw
#SBATCH -o get_variants.%j.log

## Define variables
dir=$HOME/git/somccr
reference_dir=/scratch/ucgd/lustre-work/quinlan/u1240855/somccr/data/reference
output_dir=/scratch/ucgd/lustre-work/quinlan/u1240855/somccr/data/output

## Load modules
module use $HOME/MyModules
module load miniconda3/latest
module load bedtools

#-----------------------------------------------------------------------------#

## Process TCGA variants
printf "\n\n##########################################################################\n"
printf "Processing TCGA variants... \n\n"
file=$output_dir/tcga.bed.gz
if [ ! -f "$file" ]; then ## Check if TCGA bed file does not exists
    echo "$file does not exist, getting TCGA variants..."
    python $dir/script/get_variants/python/get_tcga.py $output_dir/tcga.bed.gz
else
    echo "$file already exists, moving on..."
fi

#-----------------------------------------------------------------------------#

## Process ICGC variants
printf "\n##########################################################################\n"
printf "Processing ICGC variants... \n\n"
file=$output_dir/icgc.bed.gz
if [ ! -f "$file" ]; then ## Check if ICGC bed file does not exists
    echo "$file does not exist, getting ICGC variants..."
    python $dir/script/get_variants/python/get_icgc.py $output_dir/icgc.bed.gz
else
    echo "$file already exists, moving on..."
fi

#-----------------------------------------------------------------------------#

## Concatenate TCGA and ICGC varaints
printf "\n##########################################################################\n"
printf "Concatenating raw TCGA and ICGC variants... \n\n"
file=$output_dir/tcga_icgc.concat.bed
if [ ! -f "$file" ]; then
    echo "$file does not exist, concatenating TCGA and ICGC variants..."
    python $dir/script/get_variants/python/map_clinical_concat_filter_var.py \
    $output_dir/tcga.bed.gz \
    $output_dir/icgc.bed.gz \
    $output_dir/tcga_icgc.concat.bed
else
    echo "$file already exists, moving on..."
fi

#-----------------------------------------------------------------------------#

########################################################################################
######### Get a sorted/merged bed file of cds coordinates (hg19 from UCSC) #############
##### Refer to images in ~/git/somccr/presentations/img/*.png for more information #####
########################################################################################

## Check if the cds_coordinate file already exists
printf "\n##########################################################################\n"
printf "Checking sorted and merged hg19 UCSC CDS bed file... \n\n"
file=$reference_dir/cds/cds_ucsc_gene_coord.sorted.merged.bed
if [ ! -f "$file" ]; then ## Check if sorted & merged cds bed file does not exists
    echo "$file does not exist, sorting and merging with bedtools"
    ## Sort and merge cds interavls
    cat $reference_dir/cds/cds_ucsc_gene_coord.bed | sort -k1,1 -k2,2n | bedtools merge > $reference_dir/cds/cds_ucsc_gene_coord.sorted.merged.bed
    
    ## Remove 'chr' prefix
    cat cds_ucsc_gene_coord.sorted.merged.bed | sed 's/chr//g' > cds_ucsc_gene_coord.sorted.merged.noprefix.bed
    
    ## Rename sex chromosomes
    cat cds_ucsc_gene_coord.sorted.merged.noprefix.bed | \
    awk 'BEGIN {FS=OFS="\t"}; $1=="X" {$1=23}; {print $0}' | \
    awk 'BEGIN {FS=OFS="\t"}; $1=="Y" {$1=24}; {print $0}' > cds_ucsc_gene_coord.sorted.merged.noprefix.renameXY.bed
    
    ## Remove scaffold chromosomes
    cat cds_ucsc_gene_coord.sorted.merged.noprefix.renameXY.bed | \
    awk '/^[0-9]*\t/ {print $0}' > cds_ucsc_gene_coord.sorted.merged.noprefix.renameXY.onlyChr.bed

else
    echo "$file already exists, moving on..."
fi

#-----------------------------------------------------------------------------#

## Sort and filter variants (with bedtools)
printf "\n##########################################################################\n"
printf "Sorting and filtering TCGA and ICGC variants... \n\n"

file=$output_dir/tcga_icgc.sorted.filtered.snv.CDS.bed
if [ ! -f "$file" ]; then
    echo "$file does not exist, sorting and filtering variants..."
    cd $output_dir
    
    echo "Sorting variants..."
    cat tcga_icgc.concat.bed | sort -k1,1 -k2,2n > tcga_icgc.sorted.bed
    
    echo "Removing varaints in segdup/selfchain intervals..."
    bedtools intersect -a tcga_icgc.sorted.bed \
    -b $reference_dir/segdup_selfchain/merged_segdups_selfchain.hg19.txt -v \
    > tcga_icgc.sorted.filtered.bed
    
    echo "Selecting for SNVs..."
    cat tcga_icgc.sorted.filtered.bed | awk '{if ($3-$2==1) print $0}' > tcga_icgc.sorted.filtered.snv.bed
    
    echo "Selecting for variants in CDS regions of UCSC genes..."
    bedtools intersect -a tcga_icgc.sorted.filtered.snv.bed -b $reference_dir/cds/cds_ucsc_gene_coord.sorted.merged.bed \
    > tcga_icgc.sorted.filtered.snv.CDS.bed 
else
    echo "$file already exists, moving on..."
fi

#-----------------------------------------------------------------------------#

## Get genes with known CDS lengths (from biomaRt)
printf "\n##########################################################################\n"
printf "Obtaining genes within known CDS lengths... \n\n"
input_file=$output_dir/tcga_icgc.sorted.filtered.snv.CDS.bed
biomart_log_dir=$output_dir/biomart_log
output_file=$output_dir/gene_cds_lengths.txt
if [ ! -f $output_file ]; then
    echo "$output_file does not exist, running biomaRt to identify genes with known CDS lengths..."
    Rscript $dir/script/get_variants/Rscript/get-length.R $input_file $biomart_log_dir $output_file
else 
    echo "$output_file already exists, moving on..."
fi

#-----------------------------------------------------------------------------#

## Subset for variants in genes with known CDS lengths 
printf "\n##########################################################################\n"
printf "Subsetting variants in genes with known CDS lengths... \n\n"
input_file=$output_dir/tcga_icgc.sorted.filtered.snv.CDS.bed
known_cds_genes=$output_dir/gene_cds_lengths.txt
output_file=$output_dir/tcga_icgc.sorted.filtered.snv.CDS.knownCDSlength.bed
python $dir/script/get_variants/python/get_known_cds_var.py $input_file $known_cds_genes $output_file

#-----------------------------------------------------------------------------#

## Remove 'chr' prefix
cat $output_dir/tcga_icgc.sorted.filtered.snv.CDS.knownCDSlength.bed | \
sed 's/^chr//g' | tail -n+2 > $output_dir/tcga_icgc.sorted.filtered.snv.CDS.knownCDSlength.no_chr_prefix.bed

## Store the original header
head $output_dir/tcga_icgc.sorted.filtered.snv.CDS.knownCDSlength.bed | tail -n+2  > $output_dir/header.txt

## Replace the header of the sed-ed file 
cat $output_dir/header.txt $output_dir/tcga_icgc.sorted.filtered.snv.CDS.knownCDSlength.no_chr_prefix.bed \
> $output_dir/tcga_icgc.sorted.filtered.snv.CDS.knownCDSlength.no_chr_prefix.fix_header.bed

## Change name of sex chromosomes
cat $output_dir/tcga_icgc.sorted.filtered.snv.CDS.knownCDSlength.no_chr_prefix.fix_header.bed | \
awk 'BEGIN {FS=OFS="\t"}; $1=="X" {$1=23}; {print $0}' | \
awk 'BEGIN {FS=OFS="\t"}; $1=="Y" {$1=24}; {print $0}' \
> $output_dir/tcga_icgc.sorted.filtered.snv.CDS.knownCDSlength.no_chr_prefix.fix_header.renameXY.bed


