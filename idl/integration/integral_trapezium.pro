;Integration using Trapezium rule
;quadrature to higher order using lagrange interpolation

function testfunction,x
    return,sin(x);
end





n=20;  ;size of interval

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

//quadrature using the trapezoid rule
//accuracy to order h^2
for i=0,n-2 do begin
    sum=sum+y[i]+y[i+1];
endfor


result=sum*h/2;

print,'results is',result;

end


