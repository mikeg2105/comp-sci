
 Exercise scripts for trying out the sun grid engine
 ---------------------------------------------------
 
 This directory contains a series of semi-complete script 
 files for exercise purposes.

Compile the program using the fortran compiler

First load the compiler module
module load compilers/intel/14.0

Compile the fish program
ifort fish.f90 -o fish
 
 Start with <test1> which is a simple job script. 
 Submit test1 to the batch queue and use qstat to follow 
 its progress.
 
 After running the test1 script have a look at all the 
 other test* scripts in the numeric order. 
 It is important that you follow the numeric order as
 the problems build up using the previous solutions. 
 The comments in each script file tells you what to do next.

 This directory contains an artificial program that 
 calculates the effect of quotas on the fish stocks.
 We will practice using Sun Grid Engine by running this
 technically simple problem, the politicians prefer  
 not to know about! The program is named <fish> and resides
 in this directory.
 
 When you come to test 3, 

  Your task is to run the fish program under SGE using qsub. 
  The following steps are advisable: 
  (1) Run the program interactively to see what it does 
      and what input it expects from the user.
  (2) Create an input file which will contain the 
      expected input data ' that you would normally 
      enter at the keyboard'.
  (3) Prepare a job submission script which will 
        -read the data from your input file
        -write the results into a job output file 
         or another file of your choice.
  (4) Submit the script to SGE.
  (5) Add commands to join the error and normal output into one single file
      and also to email you about the progress of your job. 


  The file fish.txt gives background information for the program
  that is run in all three exercises 3 4 and 5. 

For test 6 you are required to modify the example test 5 
so that the job can be run as a sun grid engine task array.

Test 6 provides an opportunity to practice submitting parallel jobs
to the scheduler.

To run testmpi7 you need to compile the mpi example

Load the openmpi compiler module
module load  mpi/intel/openmpi/1.8.3

comile the diffuse program
mpicc diffuse.c -o diffuse
 
  
	 
	 
 
