FC=gfortran
FFLAGS=-O3 -Wall -Wextra -fdefault-real-8 -fcheck=all -fbacktrace
SRC=io.f90 precision.f90 area.f90 main.f90
OBJ=${SRC:.f90=.o}

%.o: %.f90
	${FC} ${FFLAGS} -o $@ -c $<
	
main: ${OBJ}
	${FC} $(FFLAGS) -o $@ $(OBJ)
	
clean:
	rm *.o *.mod main

