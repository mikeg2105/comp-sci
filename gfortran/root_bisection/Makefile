#Makefile for findrtmod2

##### User configurable options #####
include make_mac
### End User configurable options ###




OBJ = findrtmod2.o\
      rootmod.o\
      testfunctionmod.o

#Compile


testfunctionmod.o : testfunctionmod.f90
	$(FF) -c testfunctionmod.f90 $(FLAGS)
	
rootmod.o : rootmod.f90
	$(FF) -c rootmod.f90 $(FLAGS)
	
findrtmod2.o : findrtmod2.f90
	$(FF) -c findrtmod2.f90 $(FLAGS)
		

	
#all : findrt
#	$(FF) -o findrt findrtmod2.o rootmod.o testfunctionmod.o

findrt: $(OBJ)
	$(FF) -o $(@) $(OBJ)  $(FLAGS) $(LINKFLAGS)
	
.f90.o:
	$(FF) $(FLAGS) -c $<

all:findrt
	
clean:
	-rm $(OBJ) findrt *.mod


