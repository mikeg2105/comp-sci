
// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);



clear;
//clf;
clc;

// Constants
g = 9.81;
u0 = 0;
v0 = 0;
b = 0.1;
h0 = 5030;
damp=0.0001;
force=0.0;
forcefreq=0.1;
k=2.5;

// Define the x domain
ni = 1001;
xmax = 10.0;
dx = xmax/(ni-1);
x = 0:dx:xmax;
dx=0.02;

// Define the y domain
nj = 51;
ymax = 1.0;
dy = ymax/(nj-1);
y = 0:dy:ymax;

// Define the wavespeed
wavespeed = u0+sqrt(k);

// Define time-domain
dt = 10*(0.68*dx)/wavespeed;
tmax = 100;
//t = [0:dt:tdomain];
t = 1:dt:tmax;
nt=max(size(t));
courant = (wavespeed*dt)/dx;
c=0.04;
fpu=0.0;
nt=10000;
// Build empty u, v, b matrices
u = zeros(ni,3);
v = zeros(ni,1);
tv = zeros(ni,1);
uh = zeros(ni,1);
vh = zeros(ni,1);

//u((ni+1)/2)=b;
ns=2+(ni/4);
nf=ni-(ni/4)-1;
  for i = 2:ni-1
  
  // for i = ns:nf
//Y[m][prevT] = Y[m][curT] = 2.5*sech((double)(m - GRIDLENGTH/2)/5.)*
    //  Math.sin(m*10.*PI/(GRIDLENGTH-1));



        //u(i,2) = 5.0*2.5*sech((i-(ni/2))/5)*sin(i*10*%pi/(ni-1));
        u(i,2) = 2.5*sin((i)*2*%pi/((ni)-1));
       //u(i,2) = 2.5*sin((i-ns)*2*%pi/((nf-ns)-1));
        u(i,1)=u(i,2);
  

end;

u(ni,1)=0.0;
u(ni,2)=0.0;
u(ni,3)=0.0;
u(1,1)=0.0;
u(1,2)=0.0;
u(1,3)=0.0;
 
curFig             = scf(100001);
clf(curFig,"reset");
//demo_viewCode("membrane.sce");

drawlater();

xselect(); //raise the graphic window


// set a new colormap
//-------------------
cmap= curFig.color_map; //preserve old setting
curFig.color_map = jetcolormap(64);

//plot3d1(x,25,u(:,25),35,45,' ');
plot2d(x,u(:,1));
s=gce(); //the handle on the surface
//s.color_flag=1 ; //assign facet color according to Z value
title("evolution of a 3d surface","fontsize",3)



drawnow();
// Employ Lax
//for n = 2:nt-1
    for n = 2:nt
   dt=0.1; 
   t=n*dt;

    ntt=2;
    //ntt=n;
 
  
  //Fermi, Pasta, Ulam Wave Equation
   for i = 2:ni-1
      
  t1 = u(i+1,ntt) - u(i,ntt);
  t2 = u(i,ntt) - u(i-1,ntt);
    
    u(i,ntt+1) = c*c*(dt/dx)*(dt/dx)*(u(i-1,ntt) - 2*u(i,ntt) + ...
      u(i+1,ntt) + fpu*(t1*t1 - t2*t2)) - u(i,ntt-1) + 2*u(i,ntt) - ...
      damp*dt*(u(i,ntt) - u(i,ntt-1))+((i>(ni/2)-2) & (i<((ni/2)+2)) )*force*(dt^2)*sin(forcefreq*n*%pi);
  end;
          u(:,1)=u(:,2);
          u(:,2)=u(:,3);
                  
 
  // Define Boundary Conditions
  //u(1,:,n+1) = 2.5*u(2,:,n+1)-2*u(3,:,n+1)+0.5*u(4,:,n+1);
  //u(max(size(x)),:,n+1) = 2.5*u(ni-1,:,n+1)-2*u(ni-2,:,n+1)+0.5*u(ni-3,:,n+1);
  //u(:,1,n+1) = 2.5*u(:,2,n+1)-2*u(:,3,n+1)+0.5*u(:,4,n+1);
  //u(:,max(size(y)),n+1) = 2.5*u(:,nj-1,n+1)-2*u(:,nj-2,n+1)+0.5*u(:,nj-3,n+1);
// outfile=sprintf('out/outfile_%d.out',n);
//save(outfile,u,v);
 //realtime(i); //wait till date 0.1*i seconds
  //s.data.z = (sin((I(i)/10)*x)'*cos((I(i)/10)*y))';
//  s.data.y = u(:,25);
  //s.data = u(:,25);
  clf;
  xset('pixmap',1);
//plot2d(x,u(:,n),rect=[0 -4 1.0 4.0]);
plot2d(x,u(:,ntt),rect=[0 -2.5 10.0 2.5]);
//plot2d(x,u(:,n));

xset('wshow');


end;
