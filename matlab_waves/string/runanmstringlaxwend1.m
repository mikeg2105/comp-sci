

% Example of solving hydrodynamics equation to advect fluid
% Illustrates the problem of numerical instability
%
% Higher order terms for the finite element derivatives are n0t sufficient
% to remove discontinuities
%
%
% Use two step Lax-Wendroff method to stabilise solution
% http://en.wikipedia.org/wiki/Lax-Wendroff_method


%    The  time-dependent advection equation is:
%
%      du/dt +  du/dx = 0
%

%
%      du/dt + dF/dx = 0
%
%    For the advection equation,  define 
%
%      F(x,t) =  u,
%      A(x,t) = dF/dx = u
%
%    and then the Lax-Wendroff method approximates the solution 
%    using the iteration:
%
%      u(x,t+dt) = u(t) - dt dF/dx + 1/2 dt^2 d/dx A dF/dx
%
%    which can be written:
%
%      u(x,t+dt) = u(x,t) - dt ( F(x+dx,t) - F(x-dx,t) ) / ( 2 * dx )
%        + 1/2 dt^2/dx^2 ( A(x+dx/2,t) * ( F(x+dx,t) - F(x,t) )
%                        - A(x-dx/2,t) * ( F(x,t) - F(x-dx,t) )
%
%    Use the approximation:
%
%      A(x+dx/2,t) = 1/2 ( u(x+dx,t) + u(x,t) )
%      A(x-dx/2,t) = 1/2 ( u(x,t) + u(x-dx,t) )
%
%    There is a stability condition that applies here, which requires that
%
%      dt * max ( abs ( u ) ) / dx <= 1




clear;
%clf;
clc;

% Constants

u0 = 0;
v0 = 0;
b = 2;
h0 = 5030;

k=2.5;

% Define the x domain
ni = 1001;
xmax = 1.0;
dx = xmax/(ni-1);
x = 0:dx:xmax;

% Define the y domain
nj = 51;
ymax = 1.0;
dy = ymax/(nj-1);
y = 0:dy:ymax;

% Define the wavespeed
wavespeed=0.01;

% Define time-domain
dt = 50*(0.68*dx)/wavespeed;
tmax = 100;
%t = [0:dt:tdomain];
t = 1:dt:tmax;
dt=0.5
dt=0.001;
dt=0.3;
dx=10.0
courant=0.15
wavespeed=courant*dx/dt
wavespeed=10
%courant = (wavespeed*dt)/dx;
nt=10000;
hmax=3;

%parameters to define width of propagating shape
n1=75;
n2=100;
n3=125;
n3=400;
n1=30;

% Build empty u, v, b matrices
u = rand(ni,1);
v = zeros(ni,1);
tv = zeros(ni,1);
uh = zeros(ni,1);
vh = zeros(ni,1);



%u((ni+1)/2)=b;
  for i = 1:ni

        u(i) = 0.0;
        
        if i>n1
            if i<n2
               %u(i)=hmax*(i-n1)/(n2-n1);
               u(i)=hmax;
            end
        end
        
        if i>=n2
             if i<n3
                %u(i)=hmax*(n3-i)/(n3-n2);  
                u(i)=hmax;             
             end
        end
        
        ymax=1.2*hmax;
        ymin=-0.2*hmax;
      %remove comments below to set up a sine wave  
        if (i>n1) && (i<n3)
          u(i)=hmax*sin(2*pi*(i-n1)/((n3-n1)));
           ymin=-1.2*hmax;
        end
      
 

  end;
  

oldu1=u(1);
u(1)=5;
oldu2=u(2);
u(2)=-2;
h=plot(x,u(:));

hold on;
u(1)=oldu1;
u(2)=oldu2;

vw=0.1;

uold=u;


for n = 1:nt
   %dt=0.001; 
   t=n*dt;
   c=-wavespeed*dt/dx;
    pause(0.001);
    
    
    set(h,'YData',u);
    set(gca,'YLim',[ymin ymax]);

     %first order 
     starti=2;
     finishi=ni-1;
     
      %second order %uncomment these to test second order differencing
     %starti=3;
     %finishi=ni-2;
     
     %calculate hyper viscosity
     for i = starti:finishi
         
         
     end
     
     
    for i = starti:finishi
        
    %comment the above     
    %for i = 3:ni-2
      %Lax-Wendroff   
     v(i) = u(i)+c*(u(i+1)-u(i-1))/2+c^2*(u(i+1)-2*u(i)+u(i-1))/2;
     
     
    
      
      
     
      %first order central differencing
      %v(i) = u(i)+c*(u(i+1)-u(i-1))/2;%+c*(uold(i+1)-uold(i-1))/2;
  
        %second order central differencing %includes Lax-Wedroff correction
       %term
       %v(i) = u(i)+c*(u(i-2)+8*u(i+1)-8*u(i-1)-u(i+2))/12;%+c^2*(u(i+1)-2*u(i)+u(i-1))/2;
     
      
    end;
     uold=u;
     
     for i = starti:finishi
         u(i)=v(i);
      end;
      
  %1st order periodic boundary condition
  u(ni)=u(2);
  u(1)=u(ni-1);
  
  %2nd order periodic boundary condition
  %u(ni-1)=u(3);
  %u(2)=u(ni-2);
  %u(ni)=u(4);
  %u(1)=u(ni-3);
  
  
      
 
drawnow;
end;
