
## Paper_name Methods

> deployment_agave_v1.3

Repository for the analysis of mutation count variance and comparison. 

It is assumed that the a large number of directories of simulated data (each cotaining : vcf, fasta and accession to population assignment files) was prepared in advance. 
	- see [/SLiM](https://github.com/SantosJGND/SLiM)


Characteristics of the data prepared are to be stored in an `INFO_dict` file (python dictionary) within the package main directory. 

The package contains the following required items:

three scripts, namely :
- `mcount_stdlone_deploy.py`; 
- `VarInc.py`;
- `power_analysis.py`;

two dictionaries:
- `INFO_db.py`;
- `RATE_dict.py`

the directories:
- `tools`;
- `factor_tools`;


The remaining directories and files currently present within the package, `Figures` and `metadata` are produced in the course of analyses. 

## Pipeline

Two step analysis pipeline.

### I. Mutation count

> launch_mcount_stdlone_*.sh

Extract and store mutation-type counts per population and sample sizes across directories. 


`
		python -u mcount_sdtlone_dploy.py \
		--species chimp \
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
		--deployment local \
		--logdir logs
`

### II. Analyses

> launch_varInc.sh

Mutation spectrum variance and sample size.

`
		python -u VarInc_db.py \
		--species chimp \
		--data local_counts \
		--samp 100 \
		--steps 50 \
`

> launch_power_analysis.sh

Power to detect mutation-type proportion shifts.

`
		python -u power_analysis.py \
		--species chimp \
		--data local_counts \
		--collapsed \
		--norm
`


