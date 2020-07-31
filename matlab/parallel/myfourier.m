function Y = myfourier( x , L , Fs , NFFT )
% Returns the fast fourier transforms of time-domain signal x for NFFT  
% terms ( frequencies ) . Note NFFT must always be  a power of 2. 
%  x =  the original signal ( time domain )
%  L = length of signal
%  Fs = sampling frequency
%  NFFT = Number of terms in the fourier series

Y = fft(x,NFFT)/L;

%%
%the section below that use to do plotting is commented out 
% as the workers do not have plotting capabilities. 
% Plot single-sided amplitude spectrum.
% figure
% f = Fs/2*linspace(0,1,NFFT/2+1);
% plot(f,2*abs(Y(1:NFFT/2+1))) 
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')