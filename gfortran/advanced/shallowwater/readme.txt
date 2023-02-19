From Cliff Addison University of Liverpool

Mike: Attached are MPI and mixed OpenMP-MPI variants of my shallow water
benchmark code.

The latter works fine with gfortran or ifort based implementations, but
it falls over with PGI when a true mixed set-up is required (so with MPI
processes = 1, things are OK or OMP_NUM_THREADS=1, things are OK).

This code is nice because automatic parallelisation works well. I
suspect it would run nicely on a GPU with proper language support.

Hope the codes are of some help.
