%simple demo script for matlab dce
%test by Mike Griffiths 22nd January 2007
%updated by Mitra 17th October 2013 

%In old examples function findreasource has been used instead of parcluster
%function.
%resource = findResource('scheduler', 'configuration', 'sge');


%In matlab2013a parcluster is used to create cluster object
%resource = parcluster('scheduler', 'configuration', 'sge');
%set(resource, 'configuration', 'sge');

defaultProfile=parallel.defaultClusterProfile
resource=parcluster(defaultProfile);

%defaultProfile=parallel.defaultClusterProfile
%resource=parcluster(defaultProfile);
%In matlab2013a createCommunicatingJob has been replaced with
%createParallelJob.
%parjob=createParallelJob(resource);

parjob=createCommunicatingJob(resource);

% In the new release, use NumWorkersRange instead of NumWorkersRange and MaximumNumberOfWorkers
%set(parjob,'MinimumNumberOfWorkers',4);
%set(parjob,'MaximumNumberOfWorkers',4);

%Set the number of workers required for parallel execution:
parjob.NumWorkersRange = [4 4];


%createTask(parjob, 'rand', 1, {3});
%createTask(parjob, 'myrandpar', 2, {5});
createTask(parjob, 'quadpi', 1, {});

submit(parjob);

%wait for job to complete before continuing
%new release waitForState changed to wait
%waitForState(job);
wait(parjob);

% In the new release fetchOutputs is replaced with getAllOutputArguments
%parout = getAllOutputArguments(parjob);

parout = parjob.fetchOutputs;


save('results1.mat', 'parout', '-v6');
%a1=parout{1,1};
%a2=parout{1,2};
%a3=parout{2,1};
%a4=parout{2,2};
%save('jobdat1.out', 'a1','a2','a3','a4', '-ascii')

