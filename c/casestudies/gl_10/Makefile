
LIB = -L/usr/X11R6/lib64 -lGL -lX11
INCLUDE = -I/usr/include 


OBJ = cube.o


cube: $(OBJ)
	gcc -o $(@) $(OBJ) $(LIB)

.c.o:
	gcc $(INCLUDE) -c $<


