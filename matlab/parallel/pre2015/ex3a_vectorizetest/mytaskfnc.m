function y = mytaskfnc(x)
 [nr,nc] = size(x);
 y = zeros(1, nc);
 for i = 1:nc
  
  y(i) = iRosenbrock(x(:,i));
 end
 
end
 
 