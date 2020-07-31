;Radioactive Decay Model 
;Euler Method 
;Reference: Giordano pp 8-9 
;Samuel Jones 18/may/98

pro initialize,nuclei,t,tau,dt 
    print, "Initial number of nuclei:" 
    read, x 
    nuclei(0)=x 
    t(0)=0 
    print, "Time constant (tau):" 
    read, tau 
    print, "Timestep:" 
    read, dt 
end

pro calc,nuclei,t,tau,dt 
    for i=0, (n_elements(t)-2) do begin 
        nuclei(i+1)=nuclei(i)-(nuclei(i)*(dt/tau)) 
        t(i+1)=t(i)+dt 
    endfor 
end

pro displey,nuclei,t,tau,dt 
    exactsol=fltarr(100) 
    exactsol(0)=nuclei(0) 
    for i=0, (n_elements(exactsol)-2) do begin 
        exactsol(i+1)=(exactsol(0)*exp(-(i+1)*dt/tau)) 
    endfor 
    plot, t, nuclei, yrange = [0, max(nuclei)], $ 
    title='RadioActive Decay - Euler Method', xtitle = 'Time', $ 
    ytitle='Number of Nuclei', xrange=[0, max(t)],/xstyle 
    oplot, t, exactsol, line=1 
    end 
 

pro decay 
    nuclei=fltarr(100) 
    t=fltarr(100) 
    initialize,nuclei,t,tau,dt 
    calc,nuclei,t,tau,dt 
    displey,nuclei,t,tau,dt 
end