#!/bin/sh
#$ -cwd
#$ -N beatres
#$ -j y
#$ -t 1-5:1
#$ -l h_rt=00:01:00
#$ -q short.q

/usr/local/bin/matlab -nojvm -nosplash -nodisplay < beats.m
