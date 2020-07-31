% this will purely submit the beats script for batch processing.
%  This implements just one loop of beats.m serial code as a batch task.
% 
% Submit a batch j
     j = batch( 'beats' , 'CaptureDiary' , true ) 
     wait (j)
     j.State
     plotbeats 
     shg
     disp(' The listing of the job diary follow: ' )
     diary(j)
     delete ( j)
     
 
 