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

