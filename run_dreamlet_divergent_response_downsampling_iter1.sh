#! /usr/bin/env bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -N drmlt_div_dwnsmpl_i1
#$ -j y
#$ -o ./log
#$ -m ea
#$ -l mem_free=100G
#$ -l h_rt=30:00:00

conda deactivate
module load CBI
module load r
Rscript dreamlet_divergent_response_downsampling_iter1.R
