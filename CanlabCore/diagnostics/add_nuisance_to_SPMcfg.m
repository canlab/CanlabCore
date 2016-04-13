function add_nuisance_to_SPMcfg(Xn)
% Adds a matrix Xn to the end of covariates of interest in
% xX structure in SPMcfg.mat
%
% :Usage:
% ::
%
%     function add_nuisance_to_SPMcfg(Xn)
% 
% :Inputs:
%
%   **oXn:**
%        should contain ALL nuisance covariates and intercepts
%        as in output of tor_get_physio.m
%
% This function is automatically run by tor_get_physio.m
%
% ..
%    Tor Wager
%..

    if exist([pwd filesep 'SPMcfg.mat']) == 2
        load SPMcfg
    else
        warning(['No SPMcfg.mat file found in current directory:' pwd])
    end
    
    xX.X = [xX.X(:,xX.iC) Xn];

    XNamesB = xX.Xnames(xX.iB);
    
    for i = max(xX.iC)+1:size(xX.X,2)-length(XNamesB), xX.Xnames{i} = ['Nuisance ' num2str(i)];, end
    
    xX.Xnames(size(xX.X,2)-length(XNamesB)+1:size(xX.X,2)) = XNamesB;
    xX.iB = max(xX.iC)+1:size(xX.X,2);
    F_iX0.iX0 = xX.iB;
    save SPMcfg F_iX0 SPMid Sess VY xGX xM xX xsDes

    fprintf(1,'\n SPMcfg xX modified and saved.\n')
    
return
    
