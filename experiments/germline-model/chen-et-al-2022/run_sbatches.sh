#!/bin/sh

START=1
END=1000
DIR='/scratch/ucgd/lustre-work/quinlan/u0055382/SV_constraint/sbatch_scripts'

for (( i=$START; i<=$END; i++ ))
do 
    sbatch $DIR/sbatch_$i.sh
done
