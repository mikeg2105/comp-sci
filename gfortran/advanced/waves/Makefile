FC=gfortran
FFLAGS=-O3 -Wall -Wextra -fdefault-real-8 -fcheck=all -fbacktrace
SRC=lparamdata.f90 laxsolver.f90 io.f90 laxwend.f90
OBJ=${SRC:.f90=.o}

%.o: %.f90
	${FC} ${FFLAGS} -o $@ -c $<

laxwend: ${OBJ}
	${FC} $(FFLAGS) -o $@ $(OBJ)

clean:
	rm *.o *.mod laxwend

