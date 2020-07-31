#Makefile for c course

##### User configurable options #####
include make_include
### End User configurable options ###




OBJ = numerov_partialwave.o 

#Compile
pw: $(OBJ)
	$(CC) -o $(@) $(OBJ) $(INCLUDE) $(LIBS) $(FLAGS) $(LINKFLAGS)
	

#apply rule to any file appended with .c rule to generate target .o file
.c.o:
	$(CC) $(INCLUDE) $(LINKFLAGS) $(LIBS) -c $<

#another way of writing the compile rule more explicitly
#numerov_partialwave.o : numerov_partialwave.c
#	gcc -c numerov_partialwave.c

#all : game.o
#	gcc -o numerov_partialwave numerov_partialwave.o




all:pw
	
clean:
	-rm $(OBJ) pw


