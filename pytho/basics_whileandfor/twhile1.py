# example while loop with a break
#
try1 = 1;
while try1 ==1:
	n = input( ' Enter a number: ' ) ;
	if  n < 0.5:
		print 'Too small !'
	elif n > 100:
		print ' Too big !'
	else :
		print ' It will do.'
#    break
# the below line will also work for terminating the loop
# the while is terminated by the final tab
	try1 = 0;


