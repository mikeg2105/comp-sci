function y = iRosenbrock(x)
% The well-known optimization test function, the Rosenbrock function.
    y = 100*(x(2)-x(1)^2)^2+(1-x(1))^2;
 end