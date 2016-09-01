function [DX,sf] = tor_make_deconv_mtx2(sf,tp,eres)
% function [DX,sf] = tor_make_deconv_mtx(sf,tp,eres)
%   sf: cell array of stick functions, one per condition
%       all sf cells should be of the same length
%	Or matrix of stick functions, 1 column per condition 
%
%
%   tp: number of timepoints to estimate in hrf deconvolution matrix
%   eres: timebins in sf array for each TR
%
%   DX: deconvolution matrix
%       estimates O.tp time points for each condition
%       Time resolution is in TRs
%   
%   sf: stick function resampled at TR
%
%   No parametric modulation of sf's allowed.
%
% Tor Wager, 10/20/01
%warning off


if ~iscell(sf)
	%sf = mat2cell(sf,size(sf,1),ones(1,size(sf,2)));
    for i = 1:size(sf,2), sf2{i} = sf(:,i);, end
    sf = sf2;
end

shiftElements = eres;
% each time point is timeRes TRs.

% -------------------------------------------------------------------
% * downsample sf to number of TRs
% -------------------------------------------------------------------
numtrs = ceil(length(sf{1}) ./ eres);
myzeros = zeros(numtrs,1);
origsf = sf;

for i = 1:length(sf)
    Snumtrs = length(sf{i}) ./ eres;
%     if Snumtrs ~= round(Snumtrs), warning(['sf{ ' num2str(i) '}: length not evenly divisible by eres.']),end
%     if numtrs ~= ceil(Snumtrs), warning(['sf{ ' num2str(i) '}: different length than sf{1}.']),end

    inums = find(sf{i} > 0);
    if ~isempty(inums)
        inums = inums ./ eres;      % convert to TRs
        inums = ceil(inums);        % nearest TR
                                % i'm getting weird effects with round
                                % and stimuli onsets in between TRs -
                                % if the 1 in the matrix occurs before the
                                % onset of the actual event, the est. HRF is inverted -
                                % thus i'm now using ceiling.
        inums(inums == 0) = 1;      % never use 0th element
        inums(inums > numtrs) = []; % get rid of long indices
    end
    sf{i} = myzeros;
    sf{i}(inums) = 1;           % always use 1 for sf
      
end

% plot to check sampling of delta function
%figure; 
%for i = 1:length(sf)
%    subplot(length(sf),1,i);hold on
%    plot(1:1/length(origsf{i}):2-1/length(origsf{i}),origsf{i})
%    plot(1:1/length(sf{i}):2-1/length(sf{i}),sf{i},'r')
%end

% -------------------------------------------------------------------
% * make deconvolution matrix DX
% -------------------------------------------------------------------

index = 1;
for i = 1:size(sf,2)
    try
        DX(:,index) = sf{i};
    catch
        lasterr
        keyboard
    end
    index = index + 1;
    inums = find(sf{i} == 1);
    
    for j = 2:tp
        inums = inums + 1;      % + 1 because we've downsampled already.  + shiftElements;
        reg = myzeros;
        reg(inums) = 1;
        reg = reg(1:numtrs);
        while length(reg) < size(DX,1), reg = [reg;0];,end % add 0's if too short
        try
            DX(:,index) = reg;
        catch
            whos DX
            whos reg
            error('Different column lengths!')
        end
        index  = index + 1;
    end
    
end

% -------------------------------------------------------------------
% * add intercept
% -------------------------------------------------------------------
DX(:,end+1) = 1;
%warning on
return