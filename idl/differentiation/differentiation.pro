;3 and 5 point functions to compute derivatives
;need to use extrapolation to determine derivative at boundary


function diff5p,y,i,h
    nel=n_elements(y)
    if ((i ge (nel-2))  or   (i lt 2)) then begin
      print,'Inconsistant input, stop.',i,nel
      stop
    endif
    res=(y(i-2)-8*y(i-1)+8*y(i+1)-y(i+2))/(12*h);
    return,res
end

function diff3p,y,i,h
    nel=n_elements(y)
    if ((i ge (nel-1))  or   (i lt 1)) then begin
      print,'Inconsistant input, stop.'
      stop
    endif

    res=(y(i+1)-y(i-1))/(2*h);
    return,res
end

function diffdd5p,y,i,h
    nel=n_elements(y)
    if ((i ge (nel-2))  or   (i lt 2)) then begin
      print,'Inconsistant input, stop.',i,nel
      stop
    endif

    res=(-y(i-2)+16*y(i-1)+16*y(i+1)-y(i+2)-30*y(i))/(12*h^2);
    return,res
end

function diffdd3p,y,i,h
    nel=n_elements(y)
    if ((i ge (nel-1))  or   (i lt 1)) then begin
      print,'Inconsistant input, stop.'
      stop
    endif

    res=(y(i+1)+y(i-1)-2*y(i))/(h^2);
    return,res
end

;find the derivative of sin(x)



DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
PRINT, 'Date:      ', systime(0)
PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)

window, 0,xsize=700,ysize=256,XPOS = 200, YPOS = 350 

!p.multi = [0,3,1,0,1]






n=100;

x=dblarr(n)
y=dblarr(n)

for i=0,n-1 do begin
       x[i]=2*i*!pi/(n-1)
       y[i]=sin(x[i]);
endfor

h=2*!pi/(n);


ydash15p=dblarr(n)

ydash15p[1]=0;
ydash15p[2]=0;
for i=3,n-3 do begin
    ydash15p[i]=diff5p(y,i,h);
endfor

ydash15p[n-2]=0;
ydash15p[n-1]=0;

plot, x, y, charsize=1.5, title='Derivative of y - 5 point'
oplot, x, ydash15p,color=50



ydash13p=dblarr(n)
ydash13p[1]=0;
for i=2,n-2 do begin
    ydash13p[i]=diff3p(y,i,h);
endfor
ydash13p[n-1]=0;

plot, x, y, charsize=1.5, title='Derivative of y - 2 point'
oplot, x, ydash13p,color=50




yddash15p=dblarr(n)
yddash15p[1]=0;
yddash15p[2]=0;
for i=3,n-3  do begin
    yddash15p[i]=diffdd5p(y,i,h);
endfor
yddash15p[n-2]=0;
yddash15p[n-1]=0;

plot, x, y, charsize=1.5, title='Second Derivative of y - 5 point'
oplot, x, yddash15p,color=50




end
