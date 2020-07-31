#!/bin/bash
#
#  This script will run an array of jobs ( called the task arrays ) .
#   -t parameters specifies how many ( ie. 5 tasks )
#
#$ -j y
#$ -t 1-25
#$ -l h_rt=00:20:00
module load apps/matlab/02015a

matlab -nodisplay  < sgetasks.m 

