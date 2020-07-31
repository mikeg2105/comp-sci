

def wave2d(time, wavetype, maxamplitude, wavenumber,waveshift, wavefreq, delta, n):
#2d wave
#return a 2d matrix of wave amplitudes

	nx=n(1);
	ny=n(2);

	wave2d=zeros(nx,ny);
	dx=delta(1);
	dy=delta(2);

  
	k1=wavenumber(1);
	k2=wavenumber(2);
  
	ws1=waveshift(1);
	ws2=waveshift(2);
  
#stationary/standing wave 
	if wavetype == 0: 
		for i=1:nx
		for j=1:ny
			wave2d(i,j)=(sin((k1*i*dx)+(k2*j*dy))+sin((k1*(1+ws1)*i*dx)+(k2*(1+ws2)*j*dy)))*sin(wavefreq*time);
#travelling/progressive wave
	elif wavetype == 1:
		for i=1:nx
		for j=1:ny
			wave2d(i,j)=sin(((k1*i*dx)+(k2*j*dy))-(wavefreq*time))+sin(((k1*(1+ws1)*i*dx)+((k2*(1+ws2))*j*dy))-(wavefreq*time));

  
return wave2d





jobname='wav2dtest';


%Read input

wavetype=1; %travelling=1 stationary=1
nsteps=6;
maxamplitude=20;
wavenumber(1)=20;
wavenumber(2)=10;
waveshift(1)=1.4;
waveshift(2)=0.5;
wavefreq=24;
delta(1)=0.01; %
delta(2)=0.01;
nmax(1)=100;
nmax(2)=100;
deltat=0.05;
tstep=1;

outfile=[jobname,'.out'];

x=1:1:nmax(1);
y=1:1:nmax(2);

fd=fopen(outfile,'w');

%set up surface plot of output
z=zeros(nmax(1),nmax(2));
h=surf(z, 'LineStyle', 'none');

for i=tstep:tstep+nsteps
    z=wave2d(i*deltat, wavetype, maxamplitude, wavenumber,waveshift, wavefreq, delta,nmax);

    %update surface plot
    set(h,'ZData',z);
    pause;

    %Write data to output
    %the following line was used to write header an descriptive lines for dx
    %fprintf(fd, '%d %d %d\n',i, nmax(1), nmax(2));
    for j1=1:nmax(1)
      for j2=1:nmax(2)
          fprintf(fd, '%f ',z(j1,j2));
      end
      fprintf(fd, '\n');
    end
end %end of cycling over steps
fclose(fd);
%exit;

