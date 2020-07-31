%beats script
%  This is a scalar piece of MATLAB code. 
% We shall attempt to parallalize it using; 
%  (1) parfor loops  ( see beats_parf1.m  and beats_parf2.m )
%  (2) by using parfeval syntax  ( see beats_pfeval.m )
%  (3) SPMD construct ( see beats_spmd.m )
%  (4) a batch job 
%  (5) a batch job containing multiple tasks
%  (6) An array of serial tasks , by using the SGE task_arrays. 
%id=getenv('SGE_TASK_ID');  %id is a string variable
%id='4';
t=1:100;
% create four wave patterns each with different phase shifts. 
for k=1:4

    id=num2str( k );
    disp( [ 'calculating layer ' , id ] ) 
 
    wave=zeros(20,100);
    for shift=1:20
          wave( shift ,:)=beat_wave(t,1,1,0.05,k+(shift/20.0),0);
    end
    
    outfile= [ 'beat', id , '.mat' ] ;
    save(outfile, 'wave');

end;
%  Now plot the waves by calling plotbeats script.
% 
   plotbeats ;
 
 