poolSize = 4;
matlabpool(poolSize);


nsims=50;
in=rand(2,nsims);
tmpresults=zeros(1,nsims,distributor());

parfor lp=1:nsims
   tmpresults(lp)=iRosenbrock(in(:,lp)); 
end
results=gather(tmpresults,1);

matlabpool close;