Example: beats
Demonstrates the following features

1. Submission of matlab jobs to Sun Grid Engine
2. Submission of job arrays 
3. Getting environment variables from UNIX system
4. Matlab cell arrays
5. Plotting array of matrices into  surface sub plots


INSTRUCTIONS
1. Copy the matlab course examples to your account on iceberg
2. From your matlab course directory set current directory to iceberg_sgearray.
3. Submit the job to grid engine by typing qsub beats.sh
4. Check the status of your job by typing qstat -u yourusername
5. Start an interactive matlab session (on iceberg or managed windows)
6. If you are not running matlab interactively on iceberg make sure 
   then use sftp (or scp) to transfer the results from iceberg to where
   you are working.
7. Edit the script plotbeats.m and load the data from each of the jobs 
   display the plots as a series of 2x2 subplots

