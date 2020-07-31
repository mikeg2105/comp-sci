Parallel MATLAB Examples:

Almost all the examples are constructed around a script called beats ( in file beats.m ) 
                       which calls the main function beat_wave      ( in file beat_wave.m ) 

This is a good demonstration program, involving two nested for-loops which we will attempt 
to parallelize in 6 different ways thus covering all possible forms of parallelization.
These methods are :

  (1) Using parfor loops  ( see beats_parf1.m  and beats_parf2.m )
  (2) Using parfeval syntax  ( see beats_pfeval.m )
  (3) Using SPMD construct ( see beats_spmd.m )
  (4) Offloading work from MATLAB to a single worker-node via matlab batch command ( batchbeats.m ) 
  (5) Offloading work from MATLAB to a pool of workers via matlab batch command. (batchparbeats.m ) 
  (6) Running multiple SGE batch jobs, by making use of the SGE task arrays (sgetasks.sh sgetasks.m ) 
 


INSTRUCTIONS
1. Copy the matlab course examples to your account on iceberg
2. From your matlab course directory set current directory to matlab_examples/parallel.
3. Start MATLAB up.
4. Select 'local' parallel configuration.
5. Take a copy of the file beats.m and parallelize it using ;
     - parfor loops
     - using parfeval function
     - using an SPMD construct
6. Select 'sge' as your default configuration.
7. Offload the beats.m script from MATLAB as a batch job to another worker ( using batch )     
8. Offload this job to a pool of workers :  using batch ( ... 'pool'  ) 
9. Study sgetasks.sh script
10. Modify sgetask.m so that it will pickup the taskid from MATLAB script and use it in the script. 
11. Submit it as a task-array to SGE  ( i.e. qsub ... ) 
 
