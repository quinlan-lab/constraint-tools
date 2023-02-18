#!/bin/sh

START=1
n=$1
DIR='/scratch/ucgd/lustre-work/quinlan/u0055382/SV_constraint'
CCDG='/scratch/ucgd/lustre-work/quinlan/u0055382/SV_constraint/SV_data/CCDG_DELs.bed.gz'
gnomAD='/scratch/ucgd/lustre-work/quinlan/u0055382/SV_constraint/SV_data/gnomAD_DELs.bed.gz'
thousG='/scratch/ucgd/lustre-work/quinlan/u0055382/SV_constraint/SV_data/1000G_DELs.bed.gz'
gsort='gsort /dev/stdin tmp.genome'

for (( i=$START; i<=$n; i++ ))
do 
    printf '#!/bin/bash'"\n" > $DIR/sbatch_scripts/sbatch_$i.sh
    printf '#SBATCH --job-name=sbatch_process_bedtools.'$i"\n" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf '#SBATCH --account quinlan-rw'"\n" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf '#SBATCH --mem=40000'"\n" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf '#SBATCH --time=01:00:00'"\n" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf '#SBATCH --nodes=1'"\n" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf '#SBATCH --cpus-per-task=1'"\n" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf '#SBATCH -p quinlan-shared-rw'"\n" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf '#SBATCH -o '$DIR'/sbatch_logs/sbatch_process_bedtools.'$i'.%%j-%%N.out'"\n" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf '#SBATCH -e '$DIR'/sbatch_logs/sbatch_process_bedtools.'$i'.%%j-%%N.err'"\n\n" >> $DIR/sbatch_scripts/sbatch_$i.sh

    printf 'set -euo pipefail'"\n\n" >> $DIR/sbatch_scripts/sbatch_$i.sh

    printf 'python3 simulate_dataset.py '$CCDG' CCDG_'$i' | gsort /dev/stdin '$DIR'/GRCh38.autosomes_only.genome > simulated_inputs_CCDG/sim_'$i'.bed'"\n" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf 'bedtools intersect \' >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf "\n\t" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf -- '-a <(python /scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/experiments/germline-model/chen-et-al-2022/get_noncoding_windows.py) \' >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf "\n\t" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf -- '-b simulated_inputs_CCDG/sim_'$i'.bed \' >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf "\n\t" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf -- '-wao > simulated_bedtools_outputs_CCDG/intersect_'$i'.bed' >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf "\n\n" >> $DIR/sbatch_scripts/sbatch_$i.sh

    printf 'rm tmp_input.CCDG_'$i'.bed'"\n" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf 'rm tmp_shuffled.CCDG_'$i'.bed'"\n\n" >> $DIR/sbatch_scripts/sbatch_$i.sh

    printf 'python3 simulate_dataset.py '$gnomAD' gnomAD_'$i' | gsort /dev/stdin '$DIR'/GRCh38.autosomes_only.genome > simulated_inputs_gnomAD/sim_'$i'.bed'"\n" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf 'bedtools intersect \' >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf "\n\t" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf -- '-a <(python /scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/experiments/germline-model/chen-et-al-2022/get_noncoding_windows.py) \' >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf "\n\t" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf -- '-b simulated_inputs_gnomAD/sim_'$i'.bed \' >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf -- '-wao > simulated_bedtools_outputs_gnomAD/intersect_'$i'.bed' >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf "\n\n" >> $DIR/sbatch_scripts/sbatch_$i.sh

    printf 'rm tmp_input.gnomAD_'$i'.bed'"\n" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf 'rm tmp_shuffled.gnomAD_'$i'.bed'"\n\n" >> $DIR/sbatch_scripts/sbatch_$i.sh

    printf 'python3 simulate_dataset.py '$thousG' 1000G_'$i' | gsort /dev/stdin '$DIR'/GRCh38.autosomes_only.genome > simulated_inputs_1000G/sim_'$i'.bed'"\n" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf 'bedtools intersect \' >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf "\n\t" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf -- '-a <(python /scratch/ucgd/lustre-work/quinlan/u6018199/constraint-tools/experiments/germline-model/chen-et-al-2022/get_noncoding_windows.py) \' >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf "\n\t" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf -- '-b simulated_inputs_1000G/sim_'$i'.bed \' >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf "\n\t" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf -- '-wao > simulated_bedtools_outputs_1000G/intersect_'$i'.bed' >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf "\n\n" >> $DIR/sbatch_scripts/sbatch_$i.sh

    printf 'rm tmp_input.1000G_'$i'.bed'"\n" >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf 'rm tmp_shuffled.1000G_'$i'.bed'"\n\n" >> $DIR/sbatch_scripts/sbatch_$i.sh

    printf 'python3 '$DIR'/process_bedtools_intersect_wao_out.py \' >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf "\n\t" >> /scratch/ucgd/lustre-work/quinlan/u0055382/SV_constraint/sbatch_scripts/sbatch_$i.sh
    printf $DIR'/simulated_bedtools_outputs_CCDG/intersect_'$i'.bed \' >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf "\n\t" >> /scratch/ucgd/lustre-work/quinlan/u0055382/SV_constraint/sbatch_scripts/sbatch_$i.sh
    printf $DIR'/simulated_bedtools_outputs_gnomAD/intersect_'$i'.bed \' >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf "\n\t" >> /scratch/ucgd/lustre-work/quinlan/u0055382/SV_constraint/sbatch_scripts/sbatch_$i.sh
    printf $DIR'/simulated_bedtools_outputs_1000G/intersect_'$i'.bed | \' >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf "\n\t" >> /scratch/ucgd/lustre-work/quinlan/u0055382/SV_constraint/sbatch_scripts/sbatch_$i.sh
    printf 'gsort /dev/stdin/ '$DIR'/GRCh38.autosomes_only.genome > '$DIR'/simulated_final_beds/final_'$i'.bed' >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf "\n\n" >> $DIR/sbatch_scripts/sbatch_$i.sh

    printf 'python3 '$DIR'/save_fraction_windows_that_overlap_DEL__grouped_by_zscore.tom.py '$DIR'/simulated_final_beds/final_'$i'.bed '$i >> $DIR/sbatch_scripts/sbatch_$i.sh
    printf "\n" >> $DIR/sbatch_scripts/sbatch_$i.sh
done
