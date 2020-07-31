#!/bin/bash
#$ -S /bin/bash
#$ -cwd               # Run job from current directory

module load apps/R/3.2.0 # Load version 3.2.0 of R

R CMD BATCH readings.R test.R.o$JOB_ID
