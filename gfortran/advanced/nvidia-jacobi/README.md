# nvidia-jacobi
NVIDIA Multi GPU example using jacobi iteration

https://developer.nvidia.com/blog/introduction-cuda-aware-mpi/


https://on-demand.gputechconf.com/gtc/2014/presentations/S4236-multi-gpu-programming-mpi.pdf

https://www.olcf.ornl.gov/tutorials/gpudirect-mpich-enabled-cuda/



So for Bede the best approach is to use MVAPICH2/2.3.5-2. ? I presume that compiler was built using nvcc and can build any cuda code.

The idea of CUDA-aware MPI seems to be necessary, but not sufficient. My code was also a halo-code (a shallow water code). I used sendrecv with the array elements directly in one direction (associated with the row / column ordering of the language) and then an MPI data construct to handle the sendrecv in the other direction.

 

Performance was dire, but it worked correctly. The same code works fine on standard processors and Xeon Phi.

 

Someone at Nvidia wrote a 2D “Jacobi” solver with MPI and that works well, but the way they deal with the data set-up for the MPI communications seems overly complicated.

 

I have attached the Jacobi code.

 

Make is straightforward but it wants to put the code in a bin directory above the directory where the source code exists.

 

I think this should work with any of the CUDA aware MPIs on bede like  OpenMPI/3.1.4-gcccuda-2019b -  haven’t tried it there – I just used the MPI that comes with the CUDA HPC SDK on our local system and that wraps nvcc. It builds both standard and CUDA-aware MPI versions.

 

Clearly nvcc is also needed for the kernel components.

 

Nvprof can give some insights into where time is being taken. With OpenMPI you use

nvprof --annotate-mpi openmpi ./jacobi_cuda_aware_mpi -t T G in your mpirun command.

 

This is really frustrating because a single GPU runs my shallow water code at a speed akin to 100 Skylake cores and the parallelisation is standard Fortran – not even OpenMP. Try to run on more than 1 GPU and performance is several orders of magnitude lower because all of the time is spent doing copy operations for the MPI.

