import sys

## Meant to take the output from bedtools intersects
## Between chen-mchale.kmerSizes.trainSets.enhancer-exon.bed and the CCDG, gnomAD, and 1000G DEL beds
## And add the following for each region:
## total number of DELs that overlap the region
## total number of singleton DELs that overlap the region
## the Max_AF of the DELs that overlap
## fraction of the region that is found to have overlaps

def read_bedout(bedout, source):
    for line in bedout:
        fields = line.strip('\n').split('\t')
        
        ## Create an ID for each region
        # TODO: Tom, index "3" below will no longer correspond to "position", and needs to changed 
        region_id = str(fields[0]) + '.' + str(fields[3])
        
        ## Save original data from chen-mchale.kmerSizes.trainSets.enhancer-exon.bed
        if region_id not in original_data:
            # TODO: Tom, presumably the index "10" below needs to be changed?
            original_data[region_id] = fields[0:10]

        ## Initialize region_id in each of the dicts
        if region_id not in n_DELs[source]:
            n_DELs[source][region_id] = 0
        if region_id not in n_singletons[source]:
            n_singletons[source][region_id] = 0
        if region_id not in max_AFs[source]:
            max_AFs[source][region_id] = []
        if region_id not in bp_overlaps[source]:
            bp_overlaps[source][region_id] = []

        ## Identify if a region has overlaps
        ## Get and save info accordingly
        # TODO: Tom, we'lll need to edit all the indices in this code block, as necessary 
        if fields[11] == '-1': 
            continue
        elif fields[11] != '-1':
            region_size = int(fields[2]) - int(fields[1])
            af = fields[17]
            hets = fields[19]
            hom_alts = fields[20]
            bps = int(fields[21])
            fraction_bp = bps/region_size
            
            n_DELs[source][region_id] += 1
            
            if source == 'CCDG':
                ccdg_id = fields[16]
                if ccdg_id in ccdg_singletons:
                    n_singletons[source][region_id] += 1
            elif source != 'CCDG' and (hets + hom_alts == 1):
                n_singletons[source][region_id] += 1

            max_AFs[source][region_id].append(af)
            bp_overlaps[source][region_id].append(fraction_bp)

## Print out a new header
# TODO: Tom, I don't know what you'll want to change existing_header to, 
#       but the header of chen-mchale.kmerSizes.trainSets.enhancer-exon.bed is: 
#       chromosome      start   end     N_observed      N_bar_3_noncoding       N_bar_3_coding  N_bar_3_chenWindows     N_bar_5_noncoding       N_bar_5_coding  N_bar_5_chenWindows     N_bar_7_noncoding       N_bar_7_coding     N_bar_7_chenWindows     enhancer overlap        merged_exon overlap     window overlaps enhancer        window overlaps merged_exon     window overlaps (enhancer, merged_exon) new chen zscore
existing_header = ["#chromosome", "start", "end", "position", "N_bar", "N_observed", "K_bar", "K_observed", "M", "chen_zscore"]
ccdg_to_header = ["CCDG_DELs", "CCDG_singletons", "CCDG_Max_AF", "CCDG_bp_overlap"] 
gnomad_to_header = ["gnomAD_DELs", "gnomAD_singletons", "gnomAD_Max_AF", "gnomAD_bp_overlap"] 
thousandg_to_header = ["1000G_DELs", "1000G_singletons", "1000G_Max_AF", "1000G_bp_overlap"] 
header = existing_header + ccdg_to_header + gnomad_to_header + thousandg_to_header
print('\t'.join(header))

## Open the bedtools output
out_CCDG = open(sys.argv[1])
out_gnomAD = open(sys.argv[2])
out_1000G = open(sys.argv[3])

## Initialize dicts
original_data = {}            
n_DELs = {}
n_singletons = {}
max_AFs = {}
bp_overlaps = {}
        
sources = ['CCDG', 'gnomAD', '1000G']
for i in sources:
    if i not in n_DELs:
        n_DELs[i] = {}
    if i not in n_singletons:
        n_singletons[i] = {}
    if i not in max_AFs:
        max_AFs[i] = {}
    if i not in bp_overlaps:
        bp_overlaps[i] = {}

## CCDG singletons are weird so read in 
## singleton IDs from a source file
## Save them to an array
singletons = open('/scratch/ucgd/lustre-work/quinlan/u0055382/SV_constraint/SV_data/CCDG_DEL_singleton_ids.txt', 'r')
ccdg_singletons = []
for line in singletons:
        singleton_id = line.strip('\n')
        if singleton_id not in ccdg_singletons:
            ccdg_singletons.append(singleton_id)

## Actually process the data
read_bedout(out_CCDG, 'CCDG')
read_bedout(out_gnomAD, 'gnomAD')
read_bedout(out_1000G, '1000G')

## Get final numbers for Max_AFs and bp_overlaps
## Print out the data
for i in original_data:
    out = original_data[i]
    for j in sources:
        max_af = 0
        if len(max_AFs[j][i]) > 0:
            max_af = max(max_AFs[j][i])
        bp_overlap = 0
        if len(bp_overlaps[j][i]) > 0:
            bp_overlap = max(bp_overlaps[j][i])
        
        add_to_out = [str(n_DELs[j][i]), str(n_singletons[j][i]), str(max_af), str(bp_overlap)]
        out = out + add_to_out
    
    print('\t'.join(list(map(str,out))))
