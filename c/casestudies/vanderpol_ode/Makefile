##  Sample makefile
##  M.K.Griffiths
include make_include_windows
##Simple makefile

##  5/3/2004


##  Note 
##  use of ## for comments
##  use of \  for line continuation

##  Target is name of file generated by a program 
##  in this case the make utility


## To use the makefile
## To build the executable type
## make
##
## To use makefile to delete executable and 
## object files type
## make clean




## Here want to create an executable called
## target

## The rule for building the target folows immediately
## Can use CC, cc, f77, f95 to compile and link

## Note that the line follwoing the definition of the 
## target (the line defining the rule) must use tab chaacters and
## not white spaces


all : vdpodesolvers.o
	$(CC) -o vdpode vdpodesolvers.o $(LINKFLAGS)
 
## Rules for building each object file
## Note dependencies such as C/C++ header files
## contained on the lines defining the rule
## for a particular object file
## Note use of -c flag source is compiled only

vdpodesolvers.o : vdpodesolvers.c
	$(CC) -c vdpodesolvers.c $(FLAGS)



## Rule for removing executable and object files

clean :
	rm vdpode.exe vdpodesolvers.o
