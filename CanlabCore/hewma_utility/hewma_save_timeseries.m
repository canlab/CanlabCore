function cl = hewma_save_timeseries(varargin)
% :Usage:
% ::
%
%     cl = hewma_save_timeseries([mask to extract from or cl],[k],[[set of:grpmean,grpste,xdim,ydim]])
%
% load hewma_timeseries
% run this to save data in clusters format
%
% Function: extract timeseries data from a mask (and extent thr k)
% Gets group avg data for clusters, not individual.  For indiv, use
% hewma_extract_voxel.m
%
% Last args contain the data for the group for each voxel.
% If last args are not entered, attempts to load from hewma_timeseries.mat

% defaults
p = 'hewma_sig.img';    % map of significant voxels to extract
k = 10;                 % at least this many voxels

if length(varargin) > 0, p = varargin{1};,end
if length(varargin) > 1, k = varargin{2};,end

if length(varargin) > 2
    grpmean = varargin{3};
    grpste = varargin{4};
    xdim = varargin{5};
    ydim = varargin{6};
else
    fprintf(1,'Loading hewma_timeseries. ');
    load hewma_timeseries xdim ydim grpmean grpste
end

T = size(grpmean{1},2);
zdim = length(grpmean);

fprintf(1,'Reshaping. ');
grpmean = cat(1,grpmean{:});
grpste = cat(1,grpste{:});

%for i = 1:length(grpmean)
%    tmp = reshape(grpmean{i},xdim,ydim,T);
%    grpdat(:,:,i,:) = tmp;
%end

%for i = 1:length(grpste)
%    tmp = reshape(grpste{i},xdim,ydim,T);
%    stedat(:,:,i,:) = tmp;
%end

if isstruct(p)
    cl = p; % already a structure
else
    fprintf(1,'Getting clusters. ');
    cl = mask2clusters(p);
end

voxs = cat(1,cl.numVox);
whcut = find(voxs < k);
cl(whcut) = [];

% which are significant
% wh = reshape(v,xdim*ydim*zdim,1);


fprintf(1,'Saving cluster data. ');
for i = 1:length(cl)
    fprintf(1,'%3.0f.',i);
    
    vm = voxel2mask(cl(i).XYZ,[xdim ydim zdim]);
    wh = reshape(vm,xdim*ydim*zdim,1); wh = find(wh);

    cl(i).all_data = grpmean(wh,:)';
    cl(i).timeseries = nanmean(grpmean(wh,:))';
    
    cl(i).all_ste = grpste(wh,:)';
    cl(i).timeseries_ste = nanmean(grpste(wh,:))';
    
    
    %ts = timeseries4(cl(i).XYZ,grpdat);
    %cl(i).all_data = ts.indiv;
    %cl(i).timeseries = ts.avg;
end
fprintf(1,'\n');



return

