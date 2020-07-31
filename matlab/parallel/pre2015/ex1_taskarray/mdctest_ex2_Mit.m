%simple demo script for matlab dce
%test by Mike Griffiths 22nd January 2007
%updated by Mitra 10th September 2013 

%the dfeval function allows you to evaluate a function in a cluster of
%workers having to individually define jobs and tasks yourself.


%The above script file can be executed as a single command using the dfeval
%function as follows

results=dfeval(%........%)




save('myresults2','results');



for ic=1:4
    subplot(2,2,ic);
    surf(results{ic});
end