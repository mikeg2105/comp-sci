function [xm1 ,xm2 , xm3 ] = mean(x1,x2,x3,x4,x5)
% MEAN : Calculate the mean value.
% For vectors: returns a value.
% For matrices: returns the mean value of each column.
for ii = 1:nargin 
   if ii == 1 
     x=x1;
   elseif ii == 2 
     x= x2;
   elseif ii == 3 
     x= x3;
   end
   [ m , n ] = size ( x );
   if m == 1
     m = n ;
   end
   if ii == 1 
     xm1=sum(x)/m;
   elseif ii == 2 
     xm2=sum(x)/m;
   elseif ii == 3 
     xm3=sum(x)/m;
   end
end
end
