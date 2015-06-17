function [out,out_se,out_dist] = testM(M,n)
% function [out,out_se,out_dist] = testM(M,n)
%
% tests an M structure against n random permutations
% returns structures with info on the design goodness
% of the random designs.
%
% tor wager, 11/18/01

M = modelDiagnostics2(M,'noplot');

    out.energy = zeros(1,size(M.energy,2));
    out.conEnergy = zeros(1,size(M.conEnergy,2));
    out.eff_fitness = zeros(1,size(M.eff_fitness,2));
    out.eff = zeros(1,size(M.eff,2));
    out.se_contrasts = zeros(1,size(M.se_contrasts,2));
    out.resample_eff_loss = zeros(1,size(M.resample_eff_loss,2));
    out.smoothing_eff_loss = zeros(1,size(M.smoothing_eff_loss,2));
    out.vif = zeros(1,size(M.vif,2));
    out.conVif = zeros(1,size(M.conVif,2));
    out.cond = zeros(1,size(M.cond,2));
    out.conCond = zeros(1,size(M.conCond,2));
    out.conColin = zeros(size(M.conColin));
    out.hrf_eff_fitness = zeros(1,size(M.hrf_eff_fitness,2));
    out.hrf_eff = zeros(1,size(M.hrf_eff,2));
    out.hrf_eff_avg = zeros(1,size(M.hrf_eff_avg,2));
    out.cBal = zeros(1,size(M.cBal,2));
    
    out_se = out;
    out_dist = out;
    
    N = fieldnames(out);
    
    for j = N'
        eval([j{1} ' = [];']);
    end

    
    
warning off
    
for i = 1:n
    
    Mr = M;
    Mr.stimlist = getRandom(Mr.stimlist);
    
    Mr = modelDiagnostics2(Mr,'noplot');
    
    for j = N'
        eval(['out.' j{1} ' = out.' j{1} ' + Mr.' j{1} ';'])
        eval(['out_dist.' j{1} ' = [out_dist.' j{1} '; (Mr.' j{1} ')];'])
    end
    
        
    if mod(i,100) == 0, fprintf(1,'.'), end
    
end

warning on


for j = N'
    eval(['out.' j{1} ' = out.' j{1} ' / n;'])
    eval(['out_se.' j{1} ' = std(out_dist.' j{1} ') ./ sqrt(n);'])
end
    
return
    
