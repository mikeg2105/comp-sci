%Integration using trapezium rule

function y=testfunction (x)
    y=sin(x);
    return y;
%endfunction



n=20;  %size of interval


a=0;
b=pi/2;

h=(pi/2)/(n+1);

x=a:h:b;

for i=1:n+1
    y(i)=testfunction(x(i));
end

sum=0;

%quadrature using the trapezoid rule
%accuracy to order h^2
for i=1:n
    sum=sum+y(i)+y(i+1);
end

result=sum*h/2;

disp('results is',result);
