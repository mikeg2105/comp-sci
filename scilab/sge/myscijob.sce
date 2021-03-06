//function myscijob()
//env=getenv('SGE_TASK_ID');
//env='1';
//This example is used to demonstrate submission of a task array
//of scilab jobs using the EASA portal to submit an array of jobs
// to the ////White Rose grid portal 
//sgetid is define in the driver file by the EASA  application
//sgetid=sscanf(env,'%d');

	a=rand(2+sgetid,2+sgetid);
	b=rand(2+sgetid,1);
	A=sparse(a);
	
	[h,rk]=lufact(A);
	x=lusolve(h,b);
	
	res=a*x-b
	
	filename=sprintf('myscitestjob%d.out',sgetid);
	fprintfMat(filename,x);

        //define mtvars variables that will be saved to output
	myvars=['x', 'res'];

	//exit();

//endfunction
