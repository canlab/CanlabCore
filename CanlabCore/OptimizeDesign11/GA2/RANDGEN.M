function [X,bestp,eff] = randgen(gensize,conditions,mspec)
%[X,bestp,meff] = randgen(gensize,conditions,mspec)
%
% everything in conditions should be in number of frames (TRs), not seconds
%
% recommended: test conditions setup before running with
% [X,pv] = construct_model(mspec,conditions,[]);
%
% Tor Wager
% April 19, 2002

numgen = 1;
t0 = clock;
for j = 1:gensize, newvec{j} = [];,end

% intrinsic autocorrelation
Vi = eye(mspec.numframes); 

f1 = figure; set(gcf,'Color','w'); hold on; title('Max efficiency'),xlabel('Generation')

for i = 1:numgen
    
    t1 = clock;
    
    for j = 1:gensize
    
        % --------------------------------------------------------------------
        % * make models
        % --------------------------------------------------------------------
        [X,paramvec{j}] = construct_model(mspec,conditions,newvec{j});
        
        if i == i & j == 1
            cweights = ones(1,size(X,2) - 1);
        end
        
        eff(j) = calcEfficiency(cweights,[],pinv(X),Vi);
        
    end
    
    meffi = max(eff);
    meffind = find(eff == meffi); meffind = meffind(1);
    meff(i) = meffi;
    
    t2 = clock;
    fprintf(1,'Gen %3.0f: best eff is %3.4f, model %3.0f, time = %3.2f, elapsed = %3.0f\n',i,meffi,meffind,etime(t2,t1),etime(t2,t0));
    figure(f1); plot(1:i,meff,'LineWidth',2); drawnow
    
    % --------------------------------------------------------------------
    % * save best 
    % -------------------------------------------------------------------- 
    bestp = paramvec{meffind};
    
    if i == numgen, break, end
    if ~(any(eff > median(eff))), disp(['System converged at generation ' num2str(i)]), break, end
    
    
end

X = construct_model(mspec,conditions,bestp);
figure; imagesc(X); colormap(gray)


return
