
function fx=myfunc(x,q)
        G=6.67e-11;
        M=2e30;
        cs=sqrt(2*1.38e-23*4e4/(1.67e-27)); %sqrt(2kT/m)*/
        au=6.955e8; %1AU in meteres
        
        u=q;
        r=x*au;
        cs=9e4;
        
        rc=G*M/cs^2
        %fx=(u/r)*(2-(G*M/(cs*cs*r)))/((u/cs)*(u/cs)-1);
        
        fx=2*((cs/r)^2)*(r-rc)/(u*(1-(cs/u)^2))
 %endfunction;


%function fx=myfunc(x,q)

   %FUN(X,Q)=-15.915494*Q/(2-X)**2
 % fx = -15.915494*q/((2-x)^2);

%endfunction;
