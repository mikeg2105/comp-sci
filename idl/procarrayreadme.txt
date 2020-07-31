Attached is an example which shows how you can run idl under batch q

to submit to batch queue

type  
qsub getspec.sh

Here is an explanation
The idl script i'm running is
     getspecarray.pro

The getspec.sh script sets the IDL_STARTUP variable
export IDL_STARTUP=getspec_batch

getspec_batch is run as idl starts up

and note that this will execute the idl commands

.r getspecarray   (i.e. execute my script)
exit                   get out of idl

The attached idl script fetches some spectra from the SDSS database
its the first example I came across 
hope it's not too confusing

by the way the specsedscript
is just used by the idl job
