tic
xv=-1:0.1:1;
[XX,YY]=meshgrid(xv,xv);
[siz1,siz2]=size(XX);

gxx   = gpuArray(reshape(XX,siz1*siz2,1)); 
gyy   = gpuArray(reshape(YY,siz1*siz2,1));


funres=arrayfun(@mynewfun1,gxx,gyy);
results = gather(funres);
res=reshape(results,siz1,siz2);
toc
t=toc
surf(res);

