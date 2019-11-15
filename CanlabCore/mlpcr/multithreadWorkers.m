% function [workers, threads] = multithreadWorkers()
%
% If numWorkers < numCores, distributes surplus cores as threads across
% workers. Matlab has built in multithreading capabilities for some
% functions (e.g. some matrix math), so this can speed up computation when
% using a less than maximal parallel pool (which can occur when jobs don't
% easily divide into the number of CPUs).

function [workers, threads] = multithreadWorkers()
    numCores = feature('numcores');
    pool = gcp('nocreate');
    if isempty(pool)
        maxNumCompThreads(numCores);
        
        workers = 1;
        threads = numCores;
    else
        pctRunOnAll(['maxNumCompThreads(' int2str(floor(numCores/pool.NumWorkers)) ');']);
        
        workers = pool.NumWorkers;
        threads = floor(numCores/pool.NumWorkers);
    end
    
    
end
