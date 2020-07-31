


%poolSize = 4;
%matlabpool(poolSize);
tic

nsims=siz1*siz2;
xv=-1:0.1:1;
[XX,YY]=meshgrid(xv,xv);
[siz1,siz2]=size(XX);
newxx=reshape(XX,siz1*siz2,1);
newyy=reshape(YY,siz1*siz2,1);
tmpresults=zeros(1,nsims);

for lp=1:nsims
   tmpresults(lp)=mynewfun1(newxx(lp),newyy(lp)); 
end
results=reshape(tmpresults,siz1,siz2);
toc


t=toc
%matlabpool close;
surf(results);