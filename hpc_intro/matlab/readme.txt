Example: beats
Demonstrates the following features

1. Submission of matlab jobs to Sun Grid Engine
2. Submission of job arrays 
3. Running and using matlab


INSTRUCTIONS

1. Submit the job to grid engine by typing qsub beats.sh
2. Check the status of your job by typing qstat -u yourusername
3. Start an interactive matlab session (on iceberg )
4. If you are not running matlab interactively on iceberg make sure 
   then use sftp (or scp) to transfer the results from iceberg to where
   you are working.
5. From the matlab command line type the command plotbeats to display the output.
   The script plotbeats.m loads the data from each of the jobs 
   and displays data as a series of 2x2 subplots

