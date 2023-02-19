# Copyright (c) 1993-2015, NVIDIA CORPORATION. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of NVIDIA CORPORATION nor the names of its
#    contributors may be used to endorse or promote products derived
#    from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
# OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# #####
# ## Modified to run on Bede
# ## July 2021
# ## C. Addison
# ####

# 
# It is assumed that the following module file has been loaded before
# the Makefile is invoked.
# module load OpenMPI/3.1.4-gcccuda-2019b
# 
 
# Compilers
MPICC=$(PREP) mpicc
MPILD=$(PREP) mpic++
NVCC=$(PREP) nvcc
#
#
#These two environment variables have been defined consistent with:
# a) Bede environment active 08/07/2021
# b) module files used to define the CUDA libraries and MPI to use.
#MPI_HOME=/opt/software/eb/software/OpenMPI/3.1.4-gcccuda-2019b
#CUDA_INSTALL_PATH=/opt/software/eb/software/CUDA/10.1.243-GCC-8.3.0

MPI_HOME=/opt/software/builder/developers/libraries/openmpi/4.0.5/1/gcc-8.4.0-cuda-10.2.89
CUDA_INSTALL_PATH=/opt/software/builder/developers/compilers/cuda/10.2.89/1/default

# Flags - note original had -march=native, which gcc 8.3.0 complains about
CFLAGS=-std=c99 -O3 -mcpu=native -Wall
MPICFLAGS=-I${MPI_HOME}/include
CUDACFLAGS=-I${CUDA_INSTALL_PATH}/include

# 
GENCODE_SM70    := -gencode arch=compute_70,code=\"sm_70,compute_70\"
GENCODE_FLAGS   := $(GENCODE_SM70)

NVCCFLAGS=-O3 $(GENCODE_FLAGS) -Xcompiler -mcpu=native

CUDALDFLAGS=-L${CUDA_INSTALL_PATH}/lib64 -lcudart

# Description of binaries
BINDIR=$(PWD)
JACOBI_CUDA_NORMAL_MPI=$(BINDIR)/jacobi_cuda_normal_mpi
JACOBI_CUDA_AWARE_MPI=$(BINDIR)/jacobi_cuda_aware_mpi
BINARIES=$(JACOBI_CUDA_NORMAL_MPI) $(JACOBI_CUDA_AWARE_MPI)

# Commands
all: $(BINARIES)

device.o: Jacobi.h Device.cu Makefile
	$(NVCC) $(MPICFLAGS) $(NVCCFLAGS) -c Device.cu -o device.o

input.o: Jacobi.h Input.c Makefile
	$(MPICC) $(MPICFLAGS) $(CUDACFLAGS) $(CFLAGS) -c Input.c -o input.o

jacobi.o: Jacobi.h Jacobi.c Makefile
	$(MPICC) $(MPICFLAGS) $(CUDACFLAGS) $(CFLAGS) -c Jacobi.c -o jacobi.o

host.o: Jacobi.h Host.c Makefile
	$(MPICC) $(MPICFLAGS) $(CFLAGS) $(CUDACFLAGS) -c Host.c -o host.o

cuda_normal_mpi.o: Jacobi.h CUDA_Normal_MPI.c Makefile
	$(MPICC) $(MPICFLAGS) $(CFLAGS) $(CUDACFLAGS) -c CUDA_Normal_MPI.c -o cuda_normal_mpi.o

cuda_aware_mpi.o: Jacobi.h CUDA_Aware_MPI.c Makefile
	$(MPICC) $(MPICFLAGS) $(CFLAGS) $(CUDACFLAGS) -c CUDA_Aware_MPI.c -o cuda_aware_mpi.o

$(JACOBI_CUDA_NORMAL_MPI): jacobi.o input.o host.o device.o cuda_normal_mpi.o Makefile
	mkdir -p $(BINDIR)
	$(MPILD) -o $(JACOBI_CUDA_NORMAL_MPI) jacobi.o input.o host.o device.o  cuda_normal_mpi.o $(CUDALDFLAGS)
	
$(JACOBI_CUDA_AWARE_MPI): jacobi.o input.o host.o device.o cuda_aware_mpi.o Makefile
	mkdir -p $(BINDIR)
	$(MPILD) -o $(JACOBI_CUDA_AWARE_MPI) jacobi.o input.o host.o device.o cuda_aware_mpi.o $(CUDALDFLAGS)

doc: CUDA_Normal_MPI.c CUDA_Aware_MPI.c Device.cu Host.c Input.c Makefile Jacobi.c Jacobi.h Jacobi.doxygen
	doxygen Jacobi.doxygen
	
clean:
	rm -rf *.o *~ $(BINARIES)
