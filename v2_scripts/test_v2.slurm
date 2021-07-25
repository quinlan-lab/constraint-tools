#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw
#SBATCH -o v2_scripts.output

root="/scratch/ucgd/lustre-work/quinlan/u1240855/constraint-tools/v2_scripts"

cd root

bash test_v2.sh
