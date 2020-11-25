#!/bin/bash
#SBATCH -t 15 -p debug -q wildfire


module purge
module load python/3.6.4

stepup="population"
#script="$1"

python -u VarInc_db.py \
--species chimp \
--data local_counts \
--pval 0.01
