function K = getK(numsamps,ISI,TR)
% function K = getK(numsamps,ISI,TR)   
    khrf = spm_hrf(TR)';
    %khrf = khrf / max(khrf); this does weird things
    K(1,:) = [khrf zeros(1,numsamps - size(khrf,2))];
    for i = 1:numsamps-1
        row = zeros(1,i);
        row = [row khrf];
        row = [row zeros(1,numsamps - size(row,2))];
        row = row(1:numsamps);
        K(i+1,:) = row;
    end
return