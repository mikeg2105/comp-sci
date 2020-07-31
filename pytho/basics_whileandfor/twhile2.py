#simple python script
#example while statement and example prompt for input via disp command.

try1 = 1;
prompt2 = 'Try Again ! ';
while try1==1:
	n = input( ' Enter a number: ' ) ;
	if  n < 0.5:
		print 'Too small !';
		print ( prompt2);
	elif n > 100:
		print ' Too big !'
		print ( prompt2);
	else :
		print ' It will do.'
		try1=0
	
     