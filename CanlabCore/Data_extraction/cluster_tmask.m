function [cl, varargout] = cluster_tmask(cl, tm, si, varargin)
% Given clusters and a string name of a t-image,  finds voxels that exceed a
% specified t-threshold
%
% :Usage:
% ::
%
%    [cl, varargout] = cluster_tmask(cl, tm, si, varargin)
%
% :Inputs:
%
%   **cl:**
%        clusters
%   **tm:**
%        t-image
%   **si:**
%        subject index integer
%   **[dat]:**
%        cluster_barplot data structure
%
% creates new XYZ in the space of t-image using cl.XYZmm coordinates in mm.
% Required fields of cl:  XYZmm
%
% Calculates and saves single-subject data avgd over voxels if:
%    A) cl.all_data field is present
%           THIS WORKS if all_data has individual subject contrast estimates
%           in it,  with rows as subjects and columns as voxels
%           indiv data saved in cl(region).timeseries(subject)
%
%    B) cl.raw_data is present
%           raw_data should be time x voxels x subjects,  a 3D matrix
%           see output of extract_raw_data.
%           indiv data saved in cl(region).indiv_timeseries(:, subject)
%
% :Note: Retains upper 50% of voxels; highest t-values


tt = [50 100];      % % of voxeis to save,  use with read_t2; 1st number is typically used,  2nd for very small regions.

% get thresholds

%tt = [3 2];        % t-thresholds in absolute values,  use with read_t

reverse_vals = 0;
if length(varargin) > 1,  reverse_vals = varargin{2}; end

V = spm_vol(tm);

for i = 1:length(cl)
    % save mask name
    % make sure it appends it in the right place
    if ~isfield(cl(i), 'INDIV'),  cl(i).INDIV = []; end
    if ~isfield(cl(i).INDIV, 'tname'),  
        cl(i).INDIV.tname = tm;
    else
        cl(i).INDIV.tname = str2mat(cl(i).INDIV.tname, tm);
    end

    %cl(i).INDIV.tname = cl(i).INDIV.tname(1:end-1, :);
    %cl(i).INDIV.tname(si, :) = tm;
end
    
    
for i = 1:length(cl)

    % get voxel coordinates for clusters

    XYZ = mm2voxel(cl(i).XYZmm, struct('M', V.mat), 1)';

    % get t-values
    % save values of t in tm,  so that we pass in values next time to save
    % time
        [t, cl(i).INDIV.XYZ{si}, cl(i).INDIV.XYZmm{si}, cl(i).INDIV.sigt(si, :), cl(i).INDIV.maxt(si), tm,  ...
        cl(i).INDIV.center(si, :), cl(i).INDIV.mm_center(si, :)] = ...
     read_t2(tm, V, XYZ, tt, reverse_vals);   

    % Calculates and saves single-subject data avgd over voxels if:
    % A) cl.all_data field is present
    % THIS WORKS if all_data has individual subject contrast estimates in
    % it,  with rows as subjects and columns as voxels
    % indiv data saved in cl(region).timeseries(subject)
    % 
    % B) cl.raw_data is present
    % raw_data should be time x voxels x subjects,  a 3D matrix
    % see output of extract_raw_data.
    % indiv data saved in cl(region).indiv_timeseries(:, subject)
    
    if isfield(cl, 'all_data')
        try
            sigt = cl(i).INDIV.sigt(si, :);
            if any(sigt)
                cl(i).timeseries(si) = mean(cl(i).all_data(si, find(sigt)));
            else
                cl(i).timeseries(si) = NaN;
            end
        catch
            % if the data in all_data is in the wrong format (e.g., 
            % timeseries data time x voxels) we may get an error; just move
            % on.
        end
    end
  
    if isnan(cl(i).timeseries(si)),  disp('NaNs when they shouldn''t be there!'),  keyboard, end
    
    if isfield(cl, 'raw_data')
        sigt = cl(i).INDIV.sigt(si, :);
        if any(sigt)
            tmp = squeeze(cl(i).raw_data(:, :, si));
            cl(i).indiv_timeseries(:, si) = mean(tmp(:, find(sigt))')';    
        else
            cl(i).indiv_timeseries(:, si) = repmat(NaN, size(cl(i).raw_data, 1), 1);
        end
    end
    
    
    % average selected columns of data matrix entered as variable input
    % dat is cell array,  each cell is a matrix of subjects x voxels
    % argument; compatible with cluster_barplot
    
    if length(varargin) > 0
        dat = varargin{1};
        sigt = cl(i).INDIV.sigt(si, :);
        
        if ~iscell(varargin),  error('optional dat argument must be cell array,  one cell per cluster within a cell per image list'), end
        
        if any(sigt)
              datout(i) = mean(dat{i}(si, find(sigt)));
        else
            datout = NaN .* zeros(size(dat));
        end
        
        varargout{1} = datout;
    end
    
end 


return




% based on t-value threshold
function [t, XYZ, XYZmm, ts, maxt, v, com, com_mm] = read_t(tm, V, XYZ, tt)

if isstr(tm)
    v = spm_read_vols(V);
    v = v(:);
else
    v = tm;
end

ind = sub2ind(V.dim(1:3), XYZ(1, :)', XYZ(2, :)', XYZ(3, :)');

t = v(ind)';

ts = t > tt(1);

maxt = max(t(t~=0));

% lower threshold,  if necessary
if length(tt) > 1 && sum(ts) < 5     % 5 voxel cutoff
    ts = t > tt(2);
end

% save stats
if any(ts)
    XYZ = XYZ(:, find(ts));
    XYZmm = voxel2mm(XYZ, V.mat);
    com = center_of_mass(XYZ, t(find(ts)));
    com_mm = center_of_mass(XYZmm, t(find(ts)));
    %com = mean(XYZ, 2)';
    %com_mm = mean(XYZmm, 2)';
else
    XYZ = [];
    XYZmm = [];
    com = [NaN NaN NaN];
    com_mm = com;
end

return


% get voxels based on percentage (tt) of region to save
function [t, XYZ, XYZmm, ts, maxt, v, com, com_mm] = read_t2(tm, V, XYZ, tt, reverse_vals)


if isstr(tm)
    v = spm_read_vols(V);
    v = v(:);
else
    v = tm;
end

ind = sub2ind(V.dim(1:3), XYZ(1, :)', XYZ(2, :)', XYZ(3, :)');

t = v(ind)';

if reverse_vals, t = -t; end

% NaN out exactly zero voxels
t(t == 0) = NaN;

% translate percentages into t-thresholds
tt = prctile(t, 100-tt);

ts = t > tt(1);

maxt = max(t(t~=0));

% lower threshold,  if necessary,  for smaller regions
if length(tt) > 1 && sum(ts) < 5     % 5 voxel cutoff
    ts = t > tt(2);
end

% save stats
if any(ts)
    XYZ = XYZ(:, find(ts));
    XYZmm = voxel2mm(XYZ, V.mat);
    com = center_of_mass(XYZ, t(find(ts)));
    com_mm = center_of_mass(XYZmm, t(find(ts)));
    %com = mean(XYZ, 2)';
    %com_mm = mean(XYZmm, 2)';
else
    XYZ = [];
    XYZmm = [];
    com = [NaN NaN NaN];
    com_mm = com;
end

return
