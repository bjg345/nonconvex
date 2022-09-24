#!/bin/bash

#$ -l mem_free=32G
#$ -l h_vmem=32G
#$ -l h_rt=720:00:00
#$ -cwd
#$ -j y
#$ -R y
#$ -t 1:1
export _JAVA_OPTIONS=-Xmx3g

module load conda_R
Rscript test.R

