%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script to solve the heat equation in soil with
%%% seasonal temperature forcing
% See:  http://www.colorado.edu/geography/class_homepages/geog_4023_s07/labs/html/PDE_lab.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initialize variables
dz = .25; %each depth step is 1/4 meter
Nz = 400; % Choose the number of depth steps (should go to at least 100 m)
Nt = 5000; % Choose the number of time steps
dt = (365*24*60*60)/Nt; %Length of each time step in seconds (~ 6.3*10^3 seconds, or ~105 minutes)
% K =ones(1,Nz) * 2*10^-6; % "canonical" K is 10e-6 m^2/s ( range is 0.4
% to 4 ) times 10^6
K =ones(1,Nz) * 4.0*10^-6;
K= K' ;
%  K(1:200)=2*10^-6;

T = 15*ones(Nz+1,Nt+1);  %Create temperature matrix with Nz+1 rows, and Nt+1 columns
                         % Initial guess is that T is 15 everywhere.
                         % first dimension is distance , second dimension is time 
time = 0:12/Nt:12;
T(1,:)=15;
T(1,1:2000)= 300 ;  % termal shock ;
% T(1,:) = 15-20*sin(2*pi*time/12);  %Set surface temperature 
%  T(1,:) = 15-10*sin(2*pi*time/12);  %Set surface temperature
%          surface-layer temperature varying around 15C between 5C and 25C sinosoidally  
Tlast = T;
maxiter = 500 ;
err = ones(1,maxiter);
alpha=0.3 ;

for iter = 1:maxiter
    % lets use acceleration factors : T = Told + alpha( Tnew-Told) where
    % 0.0 < alpha <1 
    
    T = (1.0-alpha)*Tlast + alpha*T ;
    Tlast = T; %Save the last guess
%DS    T(:,1) = Tlast(:,end); %Initialize the temp at t=0 to the last temp
    for i=2:Nt+1,
        % d2(Temp)/dz2 is given below.
        depth_2D = (T(1:end-2,i-1)-2*T(2:end-1,i-1)+T(3:end,i-1))/dz^2;
        time_1D = K(1:Nz-1).*depth_2D;
        %  Temp(i) - Temp(i-1)    d2(Temp)
        %  ------------------  = K ------
        %      dt                   dz2
        % Therefore: 
        % Temp(i) = K*( d2(Temp)/dz2)  *dt + Temp(i-1) diffusion equation 
        T(2:end-1,i) = time_1D*dt + T(2:end-1,i-1);
%DS         T(end,i) = T(end-1,i); % Enforce bottom BC
        T(end,i) = 0.0;
    end
    err(iter) = max(abs(T(:)-Tlast(:)));  %Find difference between last two solutions
    if err(iter)<1E-4
        break; % Stop if solutions very similar, we have convergence
    end
end
if iter==maxiter;
    warning('Convergence not reached')
end
pause 
figure
plot(log(err(1:iter))), title('Convergence plot')
figure
imagesc([0 12],[0 100],T); title('Temperature plot (imagesc)')
colorbar
figure
depth = [0:dz:Nz*dz];
contourf(time,-depth,T); title('Temperature plot (contourf)')
colorbar
