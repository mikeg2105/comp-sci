%3 and 5 point functions to compute derivatives
%need to use extrapolation to determine derivative at boundary


function diff5p=diff5p(y,i,h)
    diff5p=(y(i-2)-8*y(i-1)+8*y(i+1)-y(i+2))/(12*h);
%endfunction

function diff3p=diff3p(y,i,h)
    diff3p=(y(i+1)-y(i-1))/(2*h);
%endfunction

function diffdd5p=diffdd5p(y,i,h)
    diffdd5p=(-y(i-2)+16*y(i-1)+16*y(i+1)-y(i+2)-30*y(i))/(12*h^2);
%endfunction

function diffdd3p=diffdd3p(y,i,h)
    diffdd3p=(y(i+1)+y(i-1)-2*y(i))/(h^2);
%endfunction

%find the derivative of sin(x)

n=100;
x=0.0:2*pi/(n-1):2*pi;
y=sin(x);

h=2*pi/(n);

ydash15p(1)=0;
ydash15p(2)=0;
for i=3:n-2
    ydash15p(i)=diff5p(y,i,h);
end
ydash15p(n-1)=0;
ydash15p(n)=0;

ydash13p(1)=0;
for i=2:n-1
    ydash13p(i)=diff3p(y,i,h);
end
ydash13p(n)=0;


yddash15p(1)=0;
yddash15p(2)=0;
for i=3:n-2
    yddash15p(i)=diffdd5p(y,i,h);
end
yddash15p(n-1)=0;
yddash15p(n)=0;






