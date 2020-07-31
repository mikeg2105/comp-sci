;Integration using Simpsons rule
;quadrature to higher order using lagrange interpolation

function testfunction,x
    return,sin(x);
end





n=10;  ;size of interval

;use an odd number of elements
n=(n*2)+1;

;integration interval
a=0;
b=!pi/2;

h=double((b-a)/n);
x=dblarr(n)
y=dblarr(n)

for i=0,n-1 do begin
       x[i]=a+h*i
       y[i]=testfunction(x[i]);
endfor

sum=double(0);

;since we are using interpolation
;the loop performs the summation 
;over every other element
for i=1:2:n-2 do begin
    sum=sum+y[i-1]+4*y[i]+y[i+1];
endfor

result=sum*h/3;

print,'results is',result;




end
