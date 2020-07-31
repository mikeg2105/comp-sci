The lastest version of R can be loaded with
module load apps/R/3.2.0

R can then be run with

$ R

dat <- read.csv(file = "../data/inflammation-01.csv", header = FALSE)
avg_day_inflammation <- apply(dat, 2, mean)
plot(avg_day_inflammation)

to exit R type 
q()




Serial (one CPU) Batch usage

Here, we assume that you wish to run the program my_code.R on the system.

First, you need to write a submission file. We assume you will call this rsge.sh

#!/bin/bash
#$ -S /bin/bash
#$ -cwd               # Run job from current directory

module load apps/R/3.2.0 # Load version 3.2.0 of R

R CMD BATCH readings.R test.R.o$JOB_ID

Note that R must be called with both the CMD and BATCH options which tell it to run an R program, in this case readings.R. If you do not do this, R will attempt to open an interactive prompt.

The final argument, test.R.o$JOBID, tells R to send output to a file with this name. Since $JOBID will always be unique, this ensures that all of your output files are unique. 
Without this argument R sends all output to a file readings.Rout.
Ensuring that readings.R and readings.sge are both in your current working directory, submit your job to the batch system

qsub rsge.sh

Further details on R at
http://rcg.group.shef.ac.uk/iceberg/software/apps/R.html
