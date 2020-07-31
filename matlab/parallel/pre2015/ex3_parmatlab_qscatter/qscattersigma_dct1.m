%qscatter
%partial wave analysis of scattering
%simple demo script for matlab dce
%test by Mike Griffiths 22nd January 2007


deltah=0.01;
nsteps=200;
parconfig='local';
numprocs=4;

resource = findResource('scheduler', 'configuration', parconfig);
set(resource, 'configuration', parconfig);

%set(resource, 'SubmitFcn', @sgeSubmitFcn);
%set(resource, 'ParallelSubmitFcn', @sgeParallelSubmitFcn);
%set(resource, 'ClusterMatlabRoot', '/usr/local/packages5/matlab_r2007b');
parjob=createParallelJob(resource);

set(parjob,'MinimumNumberOfWorkers',numprocs);
set(parjob,'MaximumNumberOfWorkers',numprocs);
createTask(parjob, 'qscatterpar1', 2, {1});


submit(parjob);
tic
waitForState(parjob);
myruntime=toc
parout = getAllOutputArguments(parjob);

%save('results1.mat', 'parout', '-v6');
sumouter=parout{1,1};
destroy(parjob);
plot(sumouter);

