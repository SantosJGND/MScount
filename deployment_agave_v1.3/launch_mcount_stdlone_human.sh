#!/bin/bash
#SBATCH -n 2
#SBATCH -t 00:30:00
#SBATCH --mem=8GB

module purge
module load python/3.6.4

stepup="increment" 

## if stepup == population:
# sampling is done from args.samp (below) to N (total N in vcf) in steps of args.steps (below).
# e.g. : samp= 2, steps= 1 goes from 2 to N covers all sample sizes.
## if stepup == increment:
# sampling is done from 2 to args.samp in args.steps.
# eg. samp= 2000, steps= 1998 goes from 2 to N covers all sample sizes.

## simsN == 0 uses everything. simsN > 0 randomly selects some directories.

## npacks: the number of vcfs processed by bash job. 
# e.g. npacks == 1 deploys as many jobs as there are vcf dirs.

python -u mcount_sdtlone_dploy.py \
--species human \
--data FP \
--db_dir db_dir \
-t 1:00:00 \
--nodes 6 \
--mem 20GB \
--stepup $stepup \
--simsN 0 \
--samp 2000 \
--steps 1998 \
--reps 5 \
--haps \
--npacks 1 \
--deployment agave \
--log_dir logs
