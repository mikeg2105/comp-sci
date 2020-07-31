%simple demo script for matlab dce
%test by Mike Griffiths 22nd January 2007
resource = findResource('scheduler', 'configuration', 'local');
set(resource, 'configuration', 'local');

nprocs=[1 2 3 4 5 6 8 10 12 20 24];
ntests=2;
timea=zeros(ntests,3);
%for j=1:11
for k=1:ntests
for j=4:4
numprocs=nprocs(j);
parjob=createParallelJob(resource);
set(parjob,'MaximumNumberOfWorkers',numprocs);
set(parjob,'MinimumNumberOfWorkers',numprocs);

createTask(parjob, 'transporttestpar2', 1, {1});

submit(parjob);
tic
waitForState(parjob);
myruntime=toc

parout = getAllOutputArguments(parjob);



time=0;
for i=1:numprocs
    time=time+1000*parout{i,1};
end
timea(k,:)=time'/numprocs

destroy(parjob);

end

end %looping over tests

save('jobdat1.mat', 'timea', '-v6')
%save('jobout.mat','-v6');
  

