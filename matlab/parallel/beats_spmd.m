%beats script
%  This script parallelizes the scalar MATLAB code in beats.m using SPMD model. 
% 
% Before running this code open a parallel pool with 4 workers < parpool(4) >
%  or using  <spmd(4) > 
% create four wave patterns each with different phase shifts. 

 spmd (4 , 4) 
     t=1:100;
    id=num2str( labindex ) ; 
 
    wave=zeros(20,100);
    for shift=1:20
          wave( shift ,:)=beat_wave(t,1,1,0.05,labindex+(shift/20.0),0);
    end
    
 %   outfile= [ 'beat', id , '.mat' ] ;
 %   save(outfile, 'wave');

 end
%  Now plot the waves by calling plotbeats script.
% 
   plotbeatspmd ;
 
 