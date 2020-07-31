%xval is the value of x at which we wish to interpolate
%f is the series of values
function y=lagrange_interp(xval,f,x,i)  
    t1=(xval-x(i))*(xval-x(i+1))/((x(i-1)-x(i))*(x(i-1)-x(i+1)));
    t2=(xval-x(i-1))*(xval-x(i+1))/((x(i)-x(i-1))*(x(i)-x(i+1)));
    t3=(xval-x(i-1))*(xval-x(i))/((x(i+1)-x(i-1))*(x(i+1)-x(i)));
    y=t1*f(i-1)+t2*f(i)+t3*f(i+1);


%Here is the correct working version


          % t1=(xval-x(i+1))/(x(i)-x(i-1));
          % t2=(xval-x(i))/(x(i+1)-x(i));
          % y=t1*f(i)+t2*f(i+1);   



%endfunction
