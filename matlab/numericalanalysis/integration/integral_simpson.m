%Integration using Simpsons rule

%quadrature to higher order using lagrange interpolation

function y=testfunction (x)
    y=sin(x);
endfunction

%lagrange interpolation

%xval is the value of x at which we wish to interpolate
%f is the series of values
function y=lagrange_interp(xval,f,x,i)  
    t1=(xval-x(i))*(xval-x(i+1))/((x(i-1)-x(i))*(x(i-1)-x(i+1)));
    t2=(xval-x(i-1))*(xval-x(i+1))/((x(i)-x(i-1))*(x(i)-x(i+1)));
    t3=(xval-x(i-1))*(xval-x(i))/((x(i+1)-x(i-1))*(x(i+1)-x(i)));
    y=t1*f(i-1)+t2*f(i)+t3*f(i+1);
endfunction



n=10;  %size of interval

%use an odd number of elements
n=(n*2)+1;

%integration interval
a=0;
b=pi/2;

h=(b-a)/n;

x=a:h:b;

for i=1:n+1
    y(i)=testfunction(x(i));
end

sum=0;

%since we are using interpolation
%the loop performs the summation 
%over every other element
for i=2:2:n-1
    sum=sum+y(i-1)+4*y(i)+y(i+1);
end

result=sum*h/3;

disp('results is',result);

