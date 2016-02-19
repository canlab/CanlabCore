function [paramvals,z80_ind,z80_grp] = onsets2efficiency2(ons,n,tr,hp,rho,paramvals)
% [paramvals,z80_ind] = function onsets2efficiency2(ons,n,tr,hp,rho,paramvals)


% build model
mhrf = spm_hrf(1./16,[6 16 1 1 6 0 32]); mhrf = mhrf ./ max(mhrf);
% [X] = onsets2delta(ons,1,n,mhrf);
[X] = onsets2delta(ons,n);

% get smoothing parameters
[S,V,svi,KH] = getSmoothing(hp,0,tr,n,[1 rho rho^2 rho^4]);

ind = 1;
for d = paramvals

    % get true-response model with other HRF
    hrf = spm_hrf(1./16,[6 16 1 1 6 d 32]); hrf = hrf ./ max(hrf);
%     [trueX] = onsets2delta(ons,1,n,hrf);
    trueX = onsets2delta(ons,n);

    % get power
    [z80_ind(ind),z80_grp(ind),OUT] = xpower(X,[.5 -.5],[1 0],30,.5,V,S,trueX);
    
    [z2] = xzpower(X,[.5 -.5],[1 0],30,.5,V,S,trueX);
    
    ind = ind+1;

    %d
    %sum(abs(X(:) - trueX(:)))
    %figure; plot(mhrf,'r'); hold on; plot(hrf,'k--'); title(num2str(d));
    %drawnow
end


%tor_fig; 
plot(paramvals,z80_ind,'k--','LineWidth',2); 
plot(paramvals,z2,'k-','LineWidth',2); 
xlabel('Delay'); ylabel('80% achievable Z-score');
