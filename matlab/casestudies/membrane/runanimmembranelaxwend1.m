



clear;
%clf;
clc;

% Constants
g = 9.81;
u0 = 0;
v0 = 0;
b = 0.1;
h0 = 5030;
damp=0.0;
force=5;
forcefreq=0.05;
k=2.5;
jloc=10;
% Define the x domain
ni = 20;
xmax = 10.0;
dx = xmax/(ni-1);
x = 0:dx:xmax;
dx=0.02;

% Define the y domain
nj = 20;

ymax = 10.0;
dy = ymax/(nj-1);
y = 0:dy:ymax;
dy=0.02

% Define the wavespeed
wavespeed = u0+sqrt(k);

% Define time-domain
dt = 10*(0.68*dx)/wavespeed;
tmax = 100;
%t = [0:dt:tdomain];
t = 1:dt:tmax;
nt=max(size(t));
courant = (wavespeed*dt)/dx;
c=0.04;
fpu=0.0;
nt=100000;
% Build empty u, v, b matrices
u = zeros(ni,nj,3);
v = zeros(ni,1);
tv = zeros(ni,1);
uh = zeros(ni,1);
vh = zeros(ni,1);

%u((ni+1)/2)=b;
ns=2+(ni/4);
nf=ni-(ni/4)-1;
%u(ni/2,nj/2,2)=5.0;
%u(1+ni/2,1+nj/2,2)=-5.0;
  for i = 2:ni-1
  for j = 2:nj-1
  % for i = ns:nf
%Y[m][prevT] = Y[m][curT] = 2.5*sech((double)(m - GRIDLENGTH/2)/5.)*
    %  Math.sin(m*10.*PI/(GRIDLENGTH-1));



        %u(i,2) = 5.0*2.5*sech((i-(ni/2))/5)*sin(i*10*pi/(ni-1));
       % u(i,j,2) = 2.5*sin((i)*1*pi/((ni)-1))*sin((j)*1*pi/((nj)-1));
       %u(i,2) = 2.5*sin((i-ns)*2*pi/((nf-ns)-1));
        u(i,j,1)=u(i,j,2);
  
end;
end;

u(ni,:,1)=0.0;
u(ni,:,2)=0.0;
u(ni,:,3)=0.0;
u(1,:,1)=0.0;
u(1,:,2)=0.0;
u(1,:,3)=0.0;
u(:,nj,1)=0.0;
u(:,nj,2)=0.0;
u(:,nj,3)=0.0;
u(:,1,1)=0.0;
u(:,1,2)=0.0;
u(:,1,3)=0.0; 



u(1,1,2)=-6;
u(1,2,2)=6;

h=surf(u(:,:,2));                  %Create initial surface return a handle for this graphic
axis([0 50 0 50 -10 10]);   %define axis and keep fixed 
hold on;                    %hold the graphics

%plot3d1(x,y,u(:,:,2),35,45,' ',ebox=[0 10 0 10 -6.0 6.0]);
u(1,1,2)=0;
u(1,2,2)=0;

%s=gce(); %the handle on the surface
%s.color_flag=1 ; %assign facet color according to Z value
%s.color_mode=-1;
%title("evolution of a 3d surface","fontsize",3)



drawnow;




%plot3d1(x,25,u(:,25),35,45,' ');
%plot2d(x,u(:,nj/2,1));




%drawnow();
% Employ Lax
%for n = 2:nt-1
    for n = 2:nt
   dt=0.1; 
   t=n*dt;

    ntt=2;
    %ntt=n;
 
 
  %Fermi, Pasta, Ulam Wave Equation
   for i = 2:ni-1
   for j = 2:nj-1 
       
   if  j>(nj/2) then
     c=0.005;
 else
     c=0.01;
 end   
            
  t1 = u(i+1,j,ntt) - u(i,j,ntt);
  t2 = u(i,j,ntt) - u(i-1,j,ntt);
  t3 = u(i,j+1,ntt) - u(i,j,ntt);
  t4 = u(i,j,ntt) - u(i,j-1,ntt);
   
    %u(i,j,ntt+1) = c*c*(dt/dx)*(dt/dx)*(u(i-1,j,ntt) - 4*u(i,j,ntt) + u(i,j-1,ntt)  + ...
    %  u(i+1,j,ntt) + u(i,j+1,ntt) +fpu*(t1*t1 - t2*t2+t3*t3-t4*t4)) - u(i,j,ntt-1) + 2*u(i,j,ntt) - ...
    %  damp*dt*(u(i,j,ntt) - u(i,j,ntt-1))+((i>(ni/2)-2) & (i<((ni/2)+2)) & (j>(nj/2)-2) & (j<((nj/2)+2)) )*force*(dt^2)*sin(forcefreq*n*pi);
  
  u(i,j,ntt+1) = c*c*(dt/dx)*(dt/dx)*(u(i-1,j,ntt) - 2*u(i,j,ntt)   + ...
      u(i+1,j,ntt)  +fpu*(t1*t1 - t2*t2)) - u(i,j,ntt-1) + 2*u(i,j,ntt) + ...
      c*c*(dt/dy)*(dt/dy)*(  u(i,j-1,ntt) - 2*u(i,j,ntt) + ...
       u(i,j+1,ntt) +fpu*(t3*t3-t4*t4))  - ...
      damp*dt*(u(i,j,ntt) - u(i,j,ntt-1))+((i>(ni/2)-2) & (i<((ni/2)+2)) & (j==jloc) )*force*(dt^2)*sin(forcefreq*n*pi);
  
 % if n>500 then
 %     force=0.0;
 %end;
  
      
      end;
  end;
          u(:,:,1)=u(:,:,2);
          u(:,:,2)=u(:,:,3);
                  
 
  % Define Boundary Conditions
  %u(1,:,n+1) = 2.5*u(2,:,n+1)-2*u(3,:,n+1)+0.5*u(4,:,n+1);
  %u(max(size(x)),:,n+1) = 2.5*u(ni-1,:,n+1)-2*u(ni-2,:,n+1)+0.5*u(ni-3,:,n+1);
  %u(:,1,n+1) = 2.5*u(:,2,n+1)-2*u(:,3,n+1)+0.5*u(:,4,n+1);
  %u(:,max(size(y)),n+1) = 2.5*u(:,nj-1,n+1)-2*u(:,nj-2,n+1)+0.5*u(:,nj-3,n+1);
% outfile=sprintf('out/outfile_%d.out',n);
%save(outfile,u,v);
 %realtime(i); %wait till date 0.1*i seconds
  %s.data.z = (sin((I(i)/10)*x)'*cos((I(i)/10)*y))';
%  s.data.y = u(:,25);
  %s.data = u(:,25);
  
     %  nname=sprintf('out/output%d.out',n);
     %fprintfMat(nname, u(:,:,ntt), "%lg");

  %s.data.z = u(:,:,ntt);
 
    pause(0.15);
    set(h,'ZData',u(:,:,ntt))                    %reset the z data for the grahics window with grahics handle h
    drawnow


end;
