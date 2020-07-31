% This provides a solution for a boundary-value problem for a
%first order ODE with one unknown parameter to be determined.
% The ODE ids of the form dy/dx=f(x,k) k is an unknown parameter.
% The boundary conditions are y=Ya for x=a and y=Yb for x=b


   
%  THE FUNCTION f(x,Q) IS PROVIDED AS A myfunc1(q,x) in SCILAB 
%  FILE myfunc.sci.


%Parkers fluid solar wind model*/
%We guess the initial condition and speed of solar wind when it leaves the sun*/
%The final boundary condition is known the measured speed of the solar wind 
%  at the earths orbit*/

%r/rsun=215 corresponds to earths orbit*/

%r/rsun     1.03    1.5    10    215*/
%n(m^-3)    2e14    1.6e13 4e9   7e3*/
%T(K)       1.7e6   1e6    4e5   4e4*/
%B(Teslas)  1e4     0.5e-4 1e-6  5e-9*/
%U(m/s)     0.6e3   3e3    3e5   4e5*/



%exec('myfunc.sci');
      y(1:1000)=zeros(1000,1);
      w(1:4)=zeros(4,1);
      dy(1:4)=zeros(4,1);
      dat.d1=0;
      dat.d2=0.5;
      dat.d3=0.5;
      dat.d4=1.0;

          
      
      %INPUT ESTIMATE OF THE UNKNOWN PARAMETER, Q.
      q=30.0;  %too high
      q=20;   %too low
      
      q=25.1325;
      
      
      
      
      
      
       
      q=0.8e6;  %too high
      %q=0.4e-3;  %too low      
      %q=0.6e3;
      
      
      
            %REAL Y(0:1000),W(4),DY(4)
      %CHARACTER ANS*1
      %DATA W/0.0,0.5,0.5,1.0/
      %FUN(X,Q)=-15.915494*Q/(2-X)**2
      %INPUT THE FIRST BOUNDARY CONDITION AS a,Ya
      xi=0;
      y(1)=500;
       xi=1.03;
      y(1)=0.6e3;
      %INPUT THE FINAL BOUNDARY CONDITION AS b,Yb
      %WHERE y = Yb WHEN x = b.'')')
      xf=1.0;
      yf=300.0;
      
      xf=215.0;
      yf=4e5;      
      
      
      y(1)=q;
      
      %INPUT THE INTEGRATION STEP LENGTH h IN THE 
      %FORM OF AN INTEGER N WHERE h = [b - a]/N.
      %IF YOU WANT OUTPUT OF Y AT INTERVALS OF [B-A]/M
      %THEN MAKE N A MULTIPLE OF M.
      n=1000;
      h=(xf-xi)/n;      

      
      
      
      %THE RUNGE-KUTTA INTEGRATION NOW BEGINS.F
      for i=2:n
       x=(i-1)*h;
          for j=1:4
              xx=x+w(j).*h;
              %y(i-1)=q;
              %dy(j)=h*myfunc(xx,q);
              dy(j)=h*myfunc(xx,y(i-1));
          end
       y(i)=y(i-1)+(dy(1)+dy(4)+2*(dy(2)+dy(3)))/6.0;
      end
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      