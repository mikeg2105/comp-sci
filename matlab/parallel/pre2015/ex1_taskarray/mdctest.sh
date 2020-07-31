#!/bin/sh
#$ -cwd
#$ -N taskarraytest
#$ -j y
##$ -t 1-5:1
#$ -l h_rt=00:01:00
#$ -q short.q

/usr/local/packages5/matlab_r2009a/bin/matlab -nosplash -nodisplay < mdctest_ex1.m
