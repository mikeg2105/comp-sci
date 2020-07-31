pro c20r2
;DEVICE, PSEUDO=8, DECOMPOSED=0, RETAIN=2
;WINDOW, /FREE, /PIXMAP, COLORS=256 & WDELETE, !D.WINDOW
;PRINT, 'Date:      ', systime(0)
;PRINT, 'n_colors   ', STRCOMPRESS(!D.N_COLORS,/REM)
;PRINT, 'table_size ', STRCOMPRESS(!D.TABLE_SIZE,/REM)



common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr
;loadct,5
;loadct,4
loadct,0
tvlct, r_orig, g_orig, b_orig
set_plot,'ps', /int
device, /landscape, bits=8,$
/inches, /color, filename='/data/gg1eh/c20thr/v2/PRMSL/dpabs24/c20thr_dpabs24_4xdaymean_trend_DJFM19582001.ps'

;device, /landscape, bits=8, xoffset=1, yoffset=10, xsize=9, ysize=7,$
;/inches, /color, filename='c20thr_dpabs24_4xdaymean_trend_DJFM19582001.ps'

;device, /landscape, bits=8, $
;/inches, filename='newc20thr_dpabs24_4xdaymean_trend_DJFM19582001.ps'

openu,1,'/data/gg1eh/c20thr/v2/PRMSL/dpabs24/c20thr_dpabs24_4xdaymean_trend_DJFM19582001.bin'
;openu,1,'c20thr_dpabs24_4xdaymean_trend_DJFM19582001.bin'


h1=fltarr(180,91)
readu,1,h1
close,1

print,h1

;px=!x.window*!d.x_vsize
;py=!y.window*!d.y_vsize
;sx=px(1)-px(0)
;sy=py(1)-py(0)
;POSITION=[px[0],py[0],px[0]+sx,py[0]+sy]

data = 255-bytscl(h1,min=-3,max=3,top=255)
data = rotate(data,7)




;mapStruct=MAP_PROJ_INIT(2,LIMIT=[-90,0,90,360],CENTER_LATITUDE=0)


;map_set,limit=[-90,0,90,360],$
;title='20CR dp(abs)24_4xdaymean trend: DJFM 1958-2001',$
;position=[0.094,0.089,0.89,0.9334], /grid,/noerase,color=128,/continents
map_set,limit=[-90,0,90,360],$
title='20CR dp(abs)24_4xdaymean trend: DJFM 1958-2001',$
position=[0.094,0.088,0.924,0.9334], /grid,/noerase,color=128,/continents

;position=[0.1,0.1,0.75,0.85]
;position=[0,0,1,1]

data = map_image(data,xs,ys,xsize,ysize, latmin=-90, latmax=90,$
lonmin=0, lonmax=358, compress=1)


;use rebin to smooth data

;*********************************************************;
;If the following error is reported
;REBIN: Result dimensions must be integer factor of original dimensions
;then adjust fac1 or fac2 below image dimensions are displayed by this program
fac1=1;
fac2=7;

szarr=intarr(5);
szarr=size(data);
print,'Image dimensions ',szarr(1),szarr(2),'  must be factors of ',fac1,fac2;


newdata=rebin(data,szarr(1)/fac1,szarr(2)/fac2);
newdata=congrid(newdata,szarr(1),szarr(2),/interp)
;newdata=congrid(data,szarr(1),szarr(2),/interp)



;tv, data, xs, ys, xsize=xsize, ysize=ysize





;PLOT, mapStruct.uv_box[[0,2]], mapStruct.uv_box[[1,3]],/nodata, /isotropic,xstyle=1,ystyle=1


;contour,position=POSITION,newdata , /follow , NLevels=3, C_COLORS = [120,180,240],/noerase
;contour, newdata,position=[0.05,0.05,.95,.95], /follow , NLevels=6,/noerase
contour, newdata,position=[0.094,0.088,0.924,0.9334], /follow , NLevels=6,/noerase
;contour, newdata,position=POSITION, /follow , NLevels=6,/noerase
;map_continents ,color=128
;map_continents ,/device
;use rebin to smooth data
print,size(data)



;contour, position=[0,0,24150,17800],rebin(data,483/1,356/4) , NLevels=6, /overplot , /follow 
;contour,position=POSITION, rebin(data,600/4,600/4) , NLevels=6 , /overplot, /follow,/DEVICE 



;map_grid,latdel=5,londel=10,/label


border = intarr(2,2)
;contour, border, position=[0,0,1,1],xticks=1,yticks=1,xstyle=1,ystyle=1, $
;XCharsize=0.001,YCharsize=0.001, /nodata, /noerase


levs=(findgen(101)*0.06-3)
levs2=fltarr(101,2)
levs2(*,0)=levs(*)
levs2(*,1)=levs(*)
cl=255-findgen(101)*2.5599
;contour,levs2,levs,[0,1],/cell_fill,levels=levs,position=[0.2,0.0,0.7,0.05],$
;/noerase,xticks=1,yticks=1,xstyle=1,ystyle=1, XCharsize=1.25, YCharsize=0.001,$
;c_colors=cl, xtitle = 'dp(abs)24 trend (hPa)'
;contour,levs2,levs,[0,1],levels=levs,position=[0.2,0.0,0.7,0.05],$
;/noerase,xticks=1,yticks=1,xstyle=1,ystyle=1, XCharsize=1.25, YCharsize=0.001,$
;c_colors=cl, xtitle = 'dp(abs)24 trend (hPa)'



device,/close







end
