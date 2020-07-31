%beats script
%  This is the parfeval implementation of the beats script.
% before running this start a parallel pool explicitly using parpool.
%   
%id='4';
t=1:100;
% create four wave patterns each with different phase shifts. 
pool=gcp ;
for k=1:4

    id=num2str( k );   
 
    wave=zeros(20,100);
% create 20 parallel function evaluations  to be executed by the pool of workers .... 
    for shift=1:20
         
 % %   ------        wave( shift ,:)= beat_wave(       t,1,1,0.05,k+(shift/20.0),0);
      f(shift) = parfeval( pool , @beat_wave , 1 , t,1.0, 1.0 ,0.05,(k+shift)/20.0, 0.0 );
    end
  
% fetch the outputs from the above 20 evaluations of beat_wave function....
% note that the order of output may not be in sequence ....
    for shift=1:20 
        [idfetched,value] = fetchNext(f );
        wave( idfetched ,: ) = value ; 
    end    

    outfile= [ 'beat', id , '.mat' ] ;
    save(outfile, 'wave');

end;
%  Now plot the waves by calling plotbeats script.
% 
   plotbeats ;
 
 