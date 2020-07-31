%
%  This example performs fourier series studies for analysing a wave-form
% We apply fast fourier transforms to a signal x(t) with various degrees of
% fast fourier series , in the sequence of increasing powers of 2.
% Each application is performed by using the parfeval function to pass the
% work to one of the workers. 
% 
t= 0.0:0.1:6.0 ;  % time interval is 0.0 to 6.0 


Fs = 1000;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = 1000;                     % Length of signal
t = (0:L-1)*T;                % Time vector

%% Prepare a sample signal with an artificial noise introduced to it.
%
% Pure signal is the sum of a 50 Hz sinusoid and a 120 Hz sinusoid wave.
   X = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t); 
   
%   artifical noise is normally distributed around 0.0 with amplitude=0.5 
   noise = 0.5*randn( size(t) ) ; 
   xdirty = X + noise;     % Sinusoids plus noise

%%  plot the polluted signal.. 
h1= plot(Fs*t(1:50),xdirty(1:50) ); 
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('time (milliseconds)')


%% preliminaries 
powfft = nextpow2(L) ; % Next power of 2 from length of pure signal x.
%
Y = zeros(1,2^powfft ) ; % pre allocate Y vector for efficiency reasons...

% 
%% Serial Execution :
tic; 
for ii = 3:powfft 
    NFFT = 2^ii;
    Y= myfourier( xdirty , L , Fs, NFFT ) ;
    
end
toc; 
%
%%  Parallel execution : 
 pool= gcp();  %   start the parallel pool 
% parallel evaluate myfourier function 
% each parfeval invokation returns a different FFT series.
%  defined by  2^ii terms in the series. We start with 2^3 i.e. 8 terms in
%  the series.
% 
%  Each invocation returns a different parallel-evaluation handle RES
 tic; 
 for ii = 3:powfft 
     NFFT = 2^ii;
     RES(ii) = parfeval(  pool , @myfourier ,1 , xdirty , L , Fs , NFFT ) ;
 end
%  Now we recover the results  by using the handles RES returned at each call to
%  parfeval by passing it to fetchNext . 
%  Note that fetchNext waits for a result to return. The returned result is
%  then set to the value of aresult below.  idcomplete returns the index
%  number of the parfeval call whose result is being returned. 
%  There is no guarantee that the results will be returned in the same
%  sequence as the parfeval calls are made and that is why we need the
%  idcomplete index. 

%%  Fetch results  
%fetchNext 
 for ii = 3:powfft
     [ idcomplete , aresult ] = fetchNext(RES(3:powfft) ) ;
      yresult{idcomplete} = aresult ;
 end
 toc ;
 %%  Plot Results
 hf= figure ; % start a new figure ...
 for ii = 3:powfft 
     Y = yresult{ii-2} ; 
     NFFT = 2^ii;
   
% Plot single-sided amplitude spectrum.
% figure
    f = Fs/2*linspace(0,1,NFFT/2+1);
    plot(f,2*abs(Y(1:NFFT/2+1)))   
    nffc = num2str(NFFT );
    title([ 'Single-Sided Amplitude Spectrum of y(t) using ' nffc ' terms' ])
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
    pause
end
 close(hf) ; % close the figure 
