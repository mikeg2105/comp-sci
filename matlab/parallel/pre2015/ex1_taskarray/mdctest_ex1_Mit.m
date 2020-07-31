%simple demo script for matlab dce
%test by Mike Griffiths 22nd January 2007
%updated by Mitra 10th September 2013 


%In old examples function findreasource has been used instead of parcluster
%function.
%resource = findResource('scheduler', 'configuration', config);

%In matlab2013a parcluster is used to create cluster object
%config='sge';
%resource = parcluster('scheduler', 'configuration', config);

defaultProfile=parallel.defaultClusterProfile
resource=parcluster(defaultProfile);
%The set function can be used to set different propoerties of the resource
%set(resource, 'configuration', config);

job=createJob(resource);


for i=1:4
    createTask(%.......%);
end

%submit all the tasks
submit(job);

%wait for job to complete before continuing
%new release waitForState changed to wait
%waitForState(job);
wait(job);

%The above script file can be executed as a single command using the dfeval
%function as follows
%results=dfeval(@rand,{1 2 3 4},{20 20 20 20},{100 100 100 100}, 'Configuration', 'sge')



results=getAllOutputArguments(job)
save('myresults1','results');



for ic=1:4
    subplot(2,2,ic);
    surf(results{ic});
end