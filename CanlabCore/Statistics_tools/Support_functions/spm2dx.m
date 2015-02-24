function c = spm2dx(SPM,varargin)
% c = spm2dx(SPM,varargin)
% make DX FIR deconv model from SPM (spm2 format) structure
% tor wager
%
% should work for spm99/old spm2 OR new spm2

if isstr(SPM), load(SPM),end

if length(varargin) > 0, c = varargin{1};,end

if length(varargin) > 1, wh = varargin{2};,end

% --------------------------------------------------
if exist('Sess') == 1 % SPM 99 way
    [c.ons,c.delta,c.names] = spm2delta(SPM,Sess);
    for i = 1:length(Sess), spersess(i) = length(Sess{i}.row);,end
else
    [c.ons,c.delta,c.names] = spm2delta(SPM);
    
    % this will give an error if SPM.Sess is not a cell array
    try
        for i = 1:length(SPM.Sess), spersess(i) = length(SPM.Sess{i}.row);,end  % correct?
    catch
        for i = 1:length(SPM.Sess), spersess(i) = length(SPM.Sess(i).row);,end 
    end
end

try
    c.TR = SPM.xY.RT;
catch
    % old way!?
    c.TR = xX.RT;
end

c.bf = round(32./c.TR);      % number of basis functions in FIR

c.nsess = length(c.ons) ./ size(c.delta,2);
% --------------------------------------------------


DX = tor_make_deconv_mtx3(c.delta,c.bf,1,0,1,0,spersess); % account for session cutoffs

c.numframes = repmat(c.bf,1,length(c.delta));
c.model_desc = 'DX deconv FIR model';
c.model = DX;
c.px_desc = 'pinv of c.model';
c.px = pinv(c.model); 
c.spersess = spersess;


if length(varargin) > 1, 
    % Take only certain event types of interest
    
    wh = varargin{2};,
    len = zeros(length(c.delta),1);
    len(wh) = 1;    % make indicator for events of interest
    len = repmat(len,1,c.nsess);
    c.ons = c.ons(find(len));
    
    c.delta = c.delta(wh);
    c.names = c.names(wh);
    c.numframes = c.numframes(wh);

    DX = tor_make_deconv_mtx3(c.delta,c.bf,1,0,1,0,spersess); % account for session cutoffs
    c.model = DX;
    c.px = pinv(c.model); 
end



return
