function [taskclass,maxlike,likeratio, taskprob] = classify_choose_most_likely(ptask,testvec)
% :Usage:
% ::
%
%     [taskclass,maxlike,likeratio, taskprob] = classify_choose_most_likely(ptask,testvec)
%
% :Inputs:
%
%   **ptask:**
%        likelihood of each task given activation, p(task | activation)
%        classses x variables (features, brain voxels)
%
%   **testvec:**
%        activation values across features variables x 1
%
% :Output:
%
%   **taskclass:**
%        integer for which is max likelihood class
%
%   **maxlike :**
%        log likelihood of chosen class given data (if testvec is 1/0
%               indicator)
%
%   **likeratio :**
%        likelihood ratio for most likely vs. least likely class
%
% :Examples:
% ::
%
%    % voxels in original image space that were in dataset
%    whsave = sum(indx > 0, 2) > 0;
%
%    [tc, ml, lr] = classify_choose_most_likely(ptask, MC_Setup.unweighted_study_data(whsave,1))

    ptask = log(ptask);    % likelihoods; summing these = taking product of probabilities

    taskprob = ptask * testvec;

    % choose most likely task
    [maxlike,taskclass] = max(taskprob);
    if all(taskprob == 0)
        likeratio = NaN;
    else
        likeratio = max(taskprob) ./ min(taskprob);
    end

    % get relative probabilities of being in each class; not sure about
    % scale!!
    taskprob = exp(taskprob ./ size(ptask, 2));
    taskprob = taskprob ./ sum(taskprob);
    
    return
    
