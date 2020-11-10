#!/bin/bash
#SBATCH -n 8
#SBATCH -t 06:30:00
#SBATCH --mem=34GB

module purge
module load python/3.6.4

python -u power_analysis.py \
--species chimp \
--data local_counts \
--collapsed \
--norm
