; here is an example use of the while statement
; which is used for finding the root of a polynomial 
; which is known to lie within a certain interval.
; a is the lower value of the range
; b is the upper value of the range

function sign, x
    return, (long(x gt 0) - (x lt 0))
end


a = double(0) ;fa = -%inf;
b = double(3)   ;fb = %inf;
fa = -double(9.9d99)
fb=double(9.9d99)

eps=double(0.001)

while (abs(b-a) gt (eps*b)) do begin

  x = (a+b)/2;
  fx = x^3-2*x-5;
  if (sign(fx)eq sign(fa)) then begin
    a = x;
    fa = fx;
  endif else begin
    b = x;
    fb = fx;
  endelse;

endwhile;

print, ' The root is :',x

end
