#Makefile for c course

##### User configurable options #####
include make_include
### End User configurable options ###




OBJ = arith.o 

#Compile
prog: $(OBJ)
	$(CC) -o $(@) $(OBJ) $(INCLUDE) $(LIBS) $(FLAGS) $(LINKFLAGS)
	

#apply rule to any file appended with .c rule to generate target .o file
.c.o:
	$(CC) $(INCLUDE) $(LINKFLAGS) $(LIBS) -c $<

#another way of writing the compile rule more explicitly
#game.o : game.c
#	gcc -c game.c

#all : game.o
#	gcc -o game game.o




all:prog
	
clean:
	-rm $(OBJ) prog


