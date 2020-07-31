% This script generates a wave pattern for a given set of parameters.
%  The wave pattern is saved into a file named beatsnnn.mat where
%  nnn is obtained from the sge task id.  
% Therefore each task write to a differently named file. 
% Once all the tasks are run , the output can be read back and plotted
% by using the plotall.m script. 

 id = getenv('SGE_TASK_ID') 
 ii =  str2num(id) ;
 outfile = [ 'wave' id '.mat' ]

    wave1=zeros(20,100);
    t=1:1:100;
    for shift=1:1:20;
         wave1( shift ,:)=beat_wave(t,1,1,0.05, ii+(shift/20), 0);
    end
         save(outfile, 'wave1');

% now use plotall.m to read and plot the saved waves .

