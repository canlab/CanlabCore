function [DX,sf] = tor_make_deconv_mtx3(sf,tp,eres,varargin)
% :Usage:
% ::
%
%     function [DX,sf] = tor_make_deconv_mtx(sf,tp,eres,[opt] TRs before stim onset,[num. sessions],[docenter],[scanspersess])
%
% :Inputs:
%
%   **sf:**
%        cell array of stick functions, one per condition
%        all sf cells should be of the same length
%
%        Or matrix of stick functions, 1 column per condition 
%
%   **tp:**
%        number of timepoints to estimate in hrf deconvolution matrix
%
%   **eres:**
%        timebins in sf array for each TR
%
%   **DX:**
%        deconvolution matrix
%        estimates O.tp time points for each condition
%        Time resolution is in TRs
%
% :Optional:
%
%   1. TRs before: 0 or number of time-points to shift LEFT
%   2. number of sessions; if > 1, adds session-specific intercepts
%   3. docenter, 1/0 for do/do not center columns, default 0
%   4. scanspersess: how many scans per session?  prevents regressors from
%       running over into the next session (recursive).
%
%   No parametric modulation of sf's allowed.
%
% :Outputs:
%
%   **sf:**
%        stick function resampled at TR
%
% ..
%    Tor Wager, 10/20/01   modified 9/20/02 for variable tp's for diff evt types
%    modified 4/22/04  to center columns, 2/9/05 for multi-session boundary
%    respect
% ..

docenter = 0;
if nargin > 5, docenter = varargin{3};,end
    
if ~iscell(sf)
	%sf = mat2cell(sf,size(sf,1),ones(size(sf,2)));
    for i = 1:size(sf,2), sf2{i} = sf(:,i);, end
    sf = sf2;
end

if length(tp) == 1, tp = repmat(tp,1,length(sf));, end
if length(tp) ~= length(sf), error('timepoints vectors (tp) and stick function (sf) lengths do not match!'), end

tbefore = 0;
nsess = size(sf,1);

if nargin > 4, nsess = varargin{2};, end
if nargin > 3, tbefore = varargin{1};, end

shiftElements = eres;
% each time point is timeRes TRs.




% multiple sessions
% prevent regressors from running across session lines
if nargin > 6, numframes = varargin{4};, 

    st = cumsum([1 numframes]);   
    en = st(2:end) - 1;         % ending values
    st = st(1:end-1);           % starting values

    for sess = 1:length(numframes)
        
        % get ons for this session only
        for i = 1:length(sf)
            sfsess{i} = sf{i}(st(sess):en(sess));
        end
        
        % get DX for this session only, no centering
        [DXs{sess,1}] = tor_make_deconv_mtx3(sfsess,tp,eres,varargin{1},varargin{2},0);
    end
    
    % concatenate across sessions
    DX = cat(1,DXs{:});
   
    % intercepts
    %warning('intercept term removed');
    DX = DX(:,1:end-1);                     % remove overall intercept
    DX = [DX intercept_model(repmat(numframes, 1, length(numframes)))];   % get session-specific
        
else    
    % run the single-session model


% -------------------------------------------------------------------
% * downsample sf to number of TRs
% -------------------------------------------------------------------
numtrs = round(length(sf{1}) ./ eres);
myzeros = zeros(numtrs,1);
origsf = sf;

for i = 1:length(sf)
    Snumtrs = length(sf{i}) ./ eres;
    if Snumtrs ~= round(Snumtrs), warning(['sf{ ' num2str(i) '}: length not evenly divisible by eres.']),end
    if numtrs ~= Snumtrs, warning(['sf{ ' num2str(i) '}: different length than sf{1}.']),end

    inums = find(sf{i} > 0);
    inums = inums ./ eres;      % convert to TRs
    inums = ceil(inums);        % nearest TR
                                % i'm getting weird effects with round
                                % and stimuli onsets in between TRs -
                                % if the 1 in the matrix occurs before the
                                % onset of the actual event, the est. HRF is inverted -
                                % thus i'm now using ceiling.
    inums(inums == 0) = 1;      % never use 0th element
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

    if tbefore ~= 0
	for j = tbefore:-1:1
		mysf = [sf{i}(j+1:end); zeros(j,1)];
		DX(:,index) = mysf;
		index = index + 1;
	end
    end

    DX(:,index) = sf{i};
    index = index + 1;
    inums = find(sf{i} == 1);
    
    for j = 2:tp(i)
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
if nsess < 2
    DX(:,end+1) = 1;
else
    index = 1;
    scanlen = size(DX,1) ./ nsess;
    if round(scanlen) ~= scanlen, warning('Model length is not an even multiple of scan length.'),end
  
	for startimg = 1:scanlen:size(DX,1)
		X(startimg:startimg+scanlen-1,index) = 1;
		index = index + 1;
	end
    
    DX = [DX X];
end


end     % multi-session vs. single

if docenter
    % center columns (not intercepts)
    wh = 1:size(DX,2)-nsess;
    DX(:,wh) = DX(:,wh) - repmat(mean(DX(:,wh)),size(DX,1),1);
end


return
