#Makefile for ising model

##### User configurable options #####
include make_iceberg
### End User configurable options ###




OBJ = mc.o \
      mc_system.o\
      mc_particle.o\
      mc_statistics.o

#Compile
mc: $(OBJ)
	$(CC) -o $(@) $(OBJ) $(LIBS) $(FLAGS) $(LINKFLAGS)
	
.cpp.o:
	$(CC) $(FLAGS) -c $<

all:mc
	
clean:
	-rm $(OBJ) mc


