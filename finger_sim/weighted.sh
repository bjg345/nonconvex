#!/bin/bash

#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l h_rt=720:00:00
#$ -cwd
#$ -j y
#$ -R y
#$ -t 1501:3000
module load conda_R
Rscript weighted.R $SGE_TASK_ID

