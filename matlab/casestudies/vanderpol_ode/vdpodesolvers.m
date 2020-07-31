%Van der pol equation using the Runge Kutta method
%Given:
%x'=f(t,x)
%x(0)=x0 initial condition
%and x(tn)=xn

%Use Runge Kutta to calculate xn+1 at time tn+h

%xn+1=xn+(h/6)*(a+2b+2c+d)

%a=f(tn,xn)
%b=f(tn+(h/2),xn+(h/2)a)
%c=f(tn+(h/2),xn+(h/2)b)
%d=f(tn+h,xn+h*c)

%take damped harmonic oscillator
%d2a/dt^2+ g(dq/dt)+w^2q    -1

%choose
%g=-g0(1-(q^2/q0^2))        -2

%replace 
%A=(q/q0)sqrt(g0/w)         -3
%tau=tw
%R=g0/w

%d2A/dt^2-(R-A^2)(dA/dt)+A=0  -4

%v=dA/dt                     -5

%de becomes
%dv/dt-(R-A^2)v+A=0          -6

%note:
%R<2   spiral repellor
%R>2   simple repellor
%R=0.3

%objective of the follwoing script to illustrate the convergence of runge kutta and effect
%of changing the time step


%exec('f1.sce');
%exec('f2.sce');
%exec('f3.sce');
%exec('solvers.sce');

function vanderpol()
ntests=1;
ha=zeros(1,ntests);
st=zeros(1,ntests);  %st is solver type rk4=0,  sciodes=1, adams=2, euler=3, eulerpredictorcorrector=4
st(1,1)=0;
ha(1)=0.001;
%ha(2)=0.001;
%ha(3)=0.0005;
ha(2)=0.2;
ha(3)=0.25;
ha(6)=0.0005;

h=0.1;
R=2.0;
ps0=0.5;
v0=0.3;
npoints=50000;
npoints=npoints+4; %for the adams moulton solver
gam0=1;
om=1;
ta=zeros(ntests,npoints);

ps=zeros(ntests,npoints);
ps(:,1)=ps0;

v=zeros(ntests,npoints);
v(:,1)=v0;

%t=2;
%xset('pixmap',1); % Activate pixmap
%xbasc();
%plot2d(ps(1),v(1),-1,rect=[-4,-6,4,6]);

for i=1:ntests
  h=ha(1,i);
  ta(i,1)=0;
  
  sol=st(i,1);
  if sol == 0
        [ps,v,ta,time]=rk4(ps,v,ta,h,R,gam0,om);
  elseif sol == 1 
        [ps,v,ta,time]=sciodes(ps,v,ta,h,R,gam0,om);
  elseif sol == 2
        [ps,v,ta,time]=adams(ps,v,ta,h,R,gam0,om);
  elseif sol == 3
        [ps,v,ta,time]=euler(ps,v,ta,h,R,gam0,om);
  elseif sol == 4
        [ps,v,ta,time]=eulerpredictorcorrector(ps,v,ta,h,R,gam0,om);        
  end,


  plot(ta(i,:),ps(i,:))
end

end%function


function [ps,v,ta,time]=adams(ps,v,ta,ha,R,gam0,om)
  %Runge-kutta for first three steps
  sz=size(ps);
  npoints=sz(2);
  h=ha;
  i=1;
  for j=1:4
    t=j*h;
    %ta(i,j)=t;
    rk1=f1(t,ps(i,j),v(i,j),R,gam0,om);
    rk2=f1(t+(h/2),ps(i,j),v(i,j)+(rk1*h/2),R,gam0,om);
    rk3=f1(t+(h/2),ps(i,j),v(i,j)+(rk1*h/2),R,gam0,om);
    rk4=f1(t+(h/2),ps(i,j),v(i,j)+(rk1*h/2),R,gam0,om);
    v(i,j+1)=v(i,j)+(h/6)*(rk1+2*rk2+2*rk3+rk4);

    rk1=v(i,j);
    rk2=v(i,j)+(rk1*h/2);
    rk3=v(i,j)+(rk2*h/2);
    rk4=v(i,j)+(h*rk3);
    ta(i,j+1)=t;
    ps(i,j+1)=ps(i,j)+(h/6)*(rk1+2*rk2+2*rk3+rk4);
   end;

  vam0=f1(t-(4*h),ps(i,1),v(i,1),R,gam0,om);
  vam1=f1(t-(3*h),ps(i,2),v(i,2),R,gam0,om);
  vam2=f1(t-(2*h),ps(i,3),v(i,3),R,gam0,om);

  %run adams moulton predictor-corrector steps
  tic
  for j=5:npoints-1
    t=j*h;
    ta(i,j)=t;
    vam3=f1(t-h,ps(i,j-1),v(i,j-1),R,gam0,om);
    vt=v(i,j)-(h/24)*(9*vam0-37*vam1+59*vam2-55*vam3);
    fvt=f1(t,ps(i,j),vt,R,gam0,om);
    %corrector step
    v(i,j+1)=v(i,j)+(h/24)*(vam1-5*vam2+19*vam3+9*fvt);

    fvt=f1(t,ps(i,j),v(i,j),R,gam0,om);

    vam0=vam1;
    vam1=vam2;
    vam2=vam3;



    pse=ps(i,j)-(h/24)*(9*v(i,j-4)-37*v(i,j-3)+59*v(i,j-2)-55*v(i,j-1));
    fpse=(pse-ps(i,j))/h; %numerical derivative using the estimated value
    %fpse=0;
    ps(i,j+1)=ps(i,j)+(h/24)*(v(i,j-3)-5*v(i,j-2)+19*v(i,j-1)+9*fpse);


    ta(i,j+1)=t;
   end;
  time=toc


end%function

function [ps,v,ta,time]=euler(ps,v,ta,ha,R,gam0,om)
  tic
  sz=size(ps);
  npoints=sz(2);
  h=ha;
  i=1;
  for j=1:npoints-1
    t=j*h;
    ta(i,j)=t;
    fvt=f1(t,ps(i,j),v(i,j),R,gam0,om);
    %euler step
    v(i,j+1)=v(i,j)+h*fvt;


    %euler pc step
    %v1=v(i,j)+h*fvt;
    %ps1=ps(i,j)+v1*h;
    %v(i,j+1)=v(i,j)+h*(fvt+f1(t,ps1,v1,R,gam0,om))/2;

 
    %Euler step
    ps(i,j+1)=ps(i,j)+h*v(i,j+1);

    %Euler pc step
    %p1=ps(i,j)+h*v(i,j);
    %ps(i,j+1)=ps(i,j)+h*(p1+v(i,j+1))/2;
    %ps(i,j+1)=ps(i,j)+h*v(i,j+1);

    ta(i,j+1)=t;

  end;
  toc
  time=toc




end%function

function [ps,v,ta,time]=eulerpredictorcorrector(ps,v,ta,ha,R,gam0,om)
  tic
  sz=size(ps);
  npoints=sz(2);
  h=ha;
  i=1;
  for j=1:npoints-1
    t=j*h;
    ta(i,j)=t;
    fvt=f1(t,ps(i,j),v(i,j),R,gam0,om);
 
    %euler pc step
    v1=v(i,j)+h*fvt;
    ps1=ps(i,j)+v1*h;
    v(i,j+1)=v(i,j)+h*(fvt+f1(t,ps1,v1,R,gam0,om))/2;
 
    %Euler pc step
    p1=ps(i,j)+h*v(i,j);
    ps(i,j+1)=ps(i,j)+h*(p1+v(i,j+1))/2;
    ps(i,j+1)=ps(i,j)+h*v(i,j+1);

    ta(i,j+1)=t;
  end;
  toc
  time=toc

end%function





function [ps,v,ta,time]=rk4(ps,v,ta,ha,R,gam0,om)
  tic
  sz=size(ps);
  npoints=sz(2);
  h=ha;
  i=1;
  for j=1:npoints-1
    t=j*h;
    %  fqr=(ps(i,j)/R);

    ta(i,j)=t;
    a=f1(t,ps(i,j),v(i,j),R,gam0,om);
    b=f1(t+(h/2),ps(i,j),v(i,j)+(a*h/2),R,gam0,om);
    c=f1(t+(h/2),ps(i,j),v(i,j)+(a*h/2),R,gam0,om);
    d=f1(t+(h/2),ps(i,j),v(i,j)+(a*h/2),R,gam0,om);
    

    % a=v(i,j)*gam0*(1-(fqr*fqr))-om*om*ps(i,j);
    %b=v(i,j)*gam0*(1-(fqr*fqr))-om*om*ps(i,j);
    %c=v(i,j)*gam0*(1-(fqr*fqr))-om*om*ps(i,j);
    %d=v(i,j)*gam0*(1-(fqr*fqr))-om*om*ps(i,j);   
    v(i,j+1)=v(i,j)+(h/6)*(a+2*b+2*c+d);
    a=v(i,j);
    b=v(i,j)+(a*h/2);
    c=v(i,j)+(b*h/2);
    d=v(i,j)+(h*c);
    ta(i,j+1)=t;
    ps(i,j+1)=ps(i,j)+(h/6)*(a+2*b+2*c+d);
    %plot2d(ps(t+1),v(t+1),-1,rect=[-4,-6,4,6]);
    %xset('wshow');
    %xset('wwpc');
  end;
  toc
  time=toc
%     plot(ta(i,:),ps(i,:))
  %   time=(tfinal-tinitial);


     %plot(ps(1,:),ta(2,:))
   
end%function



function [ps,v,ta,time]=sciodes(ps,v,ta,ha,R,gam0,om)
  tinitial=timer();
  for j=1:npoints-1
    t=j*h;
    %v(i,j+1)=ode("adams",v(i,j),t,t+h,list(f3,ps(i,j),R,gam0,om));
    ps(i,j+1)=ps(i,j)+h*v(i,j+1);
    %ps(i,j+1)=ode(ps(i,j),t,t+h,f2);
    ta(i,j+1)=t;
  end;
  tfinal=timer();
  time=(tfinal-tinitial);
   plot(ta(i,:),ps(i,:))
end%function

function [f1]=f1(t,ps,v,R,gam0,om)
  %f1=v*(R-ps^2)-ps;
  
  fqr=(ps/R);
  f1= v*gam0*(1-(fqr*fqr))-om*om*ps;
  %f1=v*gam0-om*om*ps;

end%function



  
  

