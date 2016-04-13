function reset_SPMcfg()
% resets columns in SPMcfg by removing all non-intercept nuisance covariates.
% runs on the SPMcfg.mat file in the current directory

    if exist([pwd filesep 'SPMcfg.mat']) == 2
        load SPMcfg
    else
        warning(['No SPMcfg.mat file found in current directory:' pwd])
    end
    
    nsess=length(Sess);
    xX.X(:,xX.iB(1:end - nsess)) = [];
    
    if length(xX.Xnames) > size(xX.X,2)
        xX.Xnames(xX.iB(1:end - nsess)) = [];
    end
    
    xX.iB = size(xX.X,2) - nsess + 1: size(xX.X,2); 
    
    F_iX0.iX0 = xX.iB;
    save SPMcfg F_iX0 SPMid Sess VY xGX xM xX xsDes

    fprintf(1,'\n SPMcfg xX modified to remove nuisance covariates and saved.\n')
    
return
    
