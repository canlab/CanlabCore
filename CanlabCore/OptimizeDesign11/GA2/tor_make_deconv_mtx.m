function [DX,timeRes] = tor_make_deconv_mtx(sf,tp,eres)
% function [DX] = tor_make_deconv_mtx(sf)
%   sf: cell array of stick functions, one per condition
%       all sf cells should be of the same length
%
%   tp: number of timepoints to estimate in hrf deconvolution matrix
%   eres: timebins in sf array to shift each estimate - e.g., 4 of 16
%         at TR = 1 gives res of .25 s.  estimates at 1st,5th,etc. 
%
%   DX: deconvolution matrix
%       estimates O.tp time points for each condition
%       each time point is timeRes TRs from last.
%
% Tor Wager, 10/20/01

shiftElements = eres;
timeRes = 1 ./ eres;
% each time point is timeRes TRs.

myzeros = zeros(size(sf{1}));
xlen = size(sf{1},1);

% -------------------------------------------------------------------
% * make deconvolution matrix DX
% -------------------------------------------------------------------

index = 1;
for i = 1:size(sf,2)
    DX(:,index) = sf{i};
    index = index + 1;
    inums = find(sf{i} == 1);
    
    for j = 2:tp
        inums = inums + shiftElements;
        reg = myzeros;
        reg(inums) = 1;
        reg = reg(1:xlen);
        DX(:,index) = reg;
        index  = index + 1;
    end
    
end



return