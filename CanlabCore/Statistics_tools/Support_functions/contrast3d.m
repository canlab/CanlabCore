function [xccon,levels] = contrast3d(xc,contrast)
% [xccon,levels] = contrast3d(xc,contrast)
%
% apply a contrast to the 3rd dimension of xc
% xc is a k x k x n matrix, where n has values to calculate contrasts over
% nested within participant 
% 
% e.g., if S = subject, t = task state
% the 3rd dim of xc might be [s1t1 s1t2 s1t3 s2t1 s2t2 s2t3] and so on.
% The contrast should have 3 elements for the 3 task states, and might be
% [1 -1 0]
% This function would return t1 - t2 matrices for each subject
%
% xccon is a 3-D matrix of contrast values
% levels is a vector of which entries in xc go with which elements of
% contrast
%
% Tor Wager
%
    ind = 1;
    levels = [];
    len = length(contrast);       % number of state matrices per subject
    
    if size(contrast,1) ~= len, contrast = contrast';, end     % make column vector
    
    
        
    for i = 1:len:size(xc,3)      % contrast must equal # of states; for each subject...
        
        tmp = xc(:,:,i:i+len-1);  % matrices for this subject
        
        tmp = permute(tmp,[1 3 2]); % make so cols are states
        
        for j = 1:size(tmp,3)
            tmp2(:,:,j) = tmp(:,:,j) * contrast;       % apply contrast to cols
        end
        
        xccon(:,:,ind) = permute(tmp2,[1 3 2]);
        ind = ind + 1;
        
        levels = [levels; (1:len)'];
        
    end
    
    
    return
    
    