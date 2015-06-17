function [ipow,gpow,ihrf,ghrf] = ga2power(M)
% ga2power(M)
%
% tor wager
% takes output from GA, M structure, and converts to power
% uses first contrast only to develop 'true' response!
%
% returns average individual and group power estimates for contrast 1 and
% HRF, avgd across trial types that show a true response in the HRF case
%
% Ignores smoothing and autocorrelation for now!!
% Fixed at 8 time points for HRF for now!!
% plots gpow
% 
% batch example:
% figure('Color','w')
% for i = 1:length(MM)
%   [dum,pow(i,1),dum,pow(i,4)] = ga2power(MM{i})
% end

[X,d] = getPredictors(M.stimlist,spm_hrf(M.ga.TR)./max(spm_hrf(M.ga.TR)));

X(:,end+1) = 1;

[ipow,gpow] = xpower(X,M.contrasts);

P = hrf_power(M.ga.TR,8,d,[1 0 0 0 0]); %M.contrasts);

numpos = sum(P.contrasts(1,:) > 0);% how many trial types show "true" responses
ghrf = sort(P.gpow); ghrf = mean(ghrf(end-numpos+1:end));

ihrf = sort(P.ipow); ghrf = mean(ihrf(end-numpos+1:end));

hold on;
plot(ghrf,gpow,'bd','MarkerFaceColor','y','LineWidth',2),