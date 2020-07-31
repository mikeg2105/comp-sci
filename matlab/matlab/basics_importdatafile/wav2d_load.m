%Example importing data from ascii file using the load command
%use with matlab script wav2dtest

%Load results and assign to cell array

load 'wav2dtest.out';
for ic=1:7
    res{ic}=wav2dtest(1+(ic-1)*100:ic*100,1:100);
end