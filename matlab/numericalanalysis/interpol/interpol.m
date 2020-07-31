%Integration using Simpsons rule

%quadrature to higher order using lagrange interpolation

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

disp('results is: ');
disp(result);

xval=1.43;
f=y;
ii=20;

val=lagrange_interp(xval,f,x,ii)




%lagrange interpolation




