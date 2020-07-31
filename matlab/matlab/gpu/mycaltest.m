
meas = ones(1000)*3; % 1000-by-1000 matrix



%The function allows the gain and offset to be arrays of the same size as rawdata, 
%so that unique corrections can be applied to individual measurements. 
%In a typical situation, you might keep the correction data on the GPU so that
%you do not have to transfer it for each application:
gn   = gpuArray(rand(1000))/100 + 0.995; 
offs = gpuArray(rand(1000))/50  - 0.01;



%Run your calibration function on the GPU:
corrected = arrayfun(@myCal, meas, gn, offs);


%This runs on the GPU because the input arguments gn and offs are already in GPU memory.
%Retrieve the corrected results from the GPU to the MATLAB workspace:
results = gather(corrected);


 