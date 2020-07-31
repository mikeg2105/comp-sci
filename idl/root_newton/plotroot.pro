; here is an example use of the while statement
; which is used for finding the root of a polynomial 
; which is known to lie within a certain interval.
; a is the lower value of the range
; b is the upper value of the range
; This is an example use of the plot function 

DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

window, 0,xsize=700,ysize=512,XPOS = 200, YPOS = 350 

!p.multi = [0,1,1,0,1]


a = double(0) ;fa = -%inf;
b = double(3)   ;fb = %inf;
x=b;
fa = -double(9.9d99)
fb=double(9.9d99)
n=100

delta=(b-a)/n
xx=dblarr(n)
fxx=dblarr(n)
it=0
for i=0,n-1 do begin
	xx[i]=a+double(i)*delta
       fxx[i]=0.5*sin(2*(xx[i]-(!pi/4)))+0.5*sin(xx[i]);
endfor


plot, xx, fxx, charsize=1.5, title='Root Determination Using Newtons Method'


eps=double(0.001)

while (sqrt((x-a)^2) gt eps) do begin

  a=x

  fx = 0.5*sin(2*(x-(!pi/4)))+0.5*sin(x);
  dfx= cos(2*(x-(!pi/4)))+0.5*cos(x);
  x=a-(fx/dfx);

  ;oplot, x,fx, color=200
  label=strTrim(string(it),1)
  xyouts,x,fx, label, charsize=1.5,color=200
  it=it+1
endwhile;

print, ' The root is :',x

end
