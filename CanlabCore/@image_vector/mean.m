function [m, varargout] = mean(obj, varargin)
% Mean across a set of images. Returns a new image_vector object.
% Creates an image_vector object with mean values for each voxel (cols)
% across images (rows) of an image_vector (e.g., fmri_data) object.
%
% :Usage:
% ::
%
%    function [m, imagemeans, voxelmeans]  = mean(obj, [optional args])
%
% m is an image_vector object whose data contains the mean values.
%
% - Average available valid data in each voxel. Some images may have
% missing data for some voxels.
% - Treats values of 0 as missing (after SPM) and excludes from mean.
%
% :Optional Inputs:
%   - 'group', 'group_by' -> followed by vector of integers to define
%       groups (e.g., images from the same participant). Averages across images
%       within each group.
%   - 'write', followed by file name
%   - 'path', followed by location for file (default = current directory)
%   - 'orthviews' -> show orthviews for this image, same as orthviews(m)
%   - 'histogram' -> show histogram for this image, same as histogram(m)
%   - 'plot' -> do both
%
% :Examples:
% ::
%
%    % If sdat is an fmri_data object with multiple images,
%    m = mean(sdat, 'plot', 'write', anatmeanname, 'path', maskdir);
%

% Programmers' notes:
% - Tor: Average available valid data in each voxel (June 2017)
% - Michael Sun, PhD: Zero may actually be data, so NaN-ing them can be
% destructive, so add a flag 'treat_zero_as_data' 04/21/2025

fname = [];
fpath = pwd;
dohist = false;
doorth = false;
doplot = false;
dogroup = false;
group_by = [];
u = [];
treat_zero_as_data = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case 'write', fname = varargin{i+1}; varargin{i+1} = [];
            case 'path', fpath = varargin{i+1}; varargin{i+1} = [];
            case 'plot', dohist = 1; doorth = 1;

            case {'hist', 'histogram'}, dohist = 1;
            case {'orth', 'orthviews'}, doorth = 1;

            case {'group', 'group_by'}
                dogroup = true;
                group_by = varargin{i+1};
                group_by = categorical(group_by);  % this should work for strings or numeric arrays
                u = unique(group_by);
            case {'treat_zero_as_data'}, treat_zero_as_data = 1;

            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


% get missing values and replace with NaNs, so we average only valid data in each
% voxel
% wh = obj.dat == 0; % NaNs are already OK | isnan(obj.dat);
% obj.dat(wh) = NaN;

% This is potentially dangerous as 0 may be data - Michael Sun 05/28/2025
if ~treat_zero_as_data
    wh = obj.dat == 0; % NaNs are already OK | isnan(obj.dat);
    obj.dat(wh) = NaN;
end


% return output in the same format as input object

if isa(obj, 'fmri_data')

    if dogroup
        for i = 1:length(u)
            mydat(:, i) = nanmean(obj.dat(:, group_by == u(i))', 1)';
        end

    else
        mydat = nanmean(obj.dat', 1)';
    end

    m = image_vector('dat', mydat, 'volInfo', obj.mask.volInfo, 'noverbose');

    m = fmri_data(m, [], 'noverbose');
    m.mask = obj.mask;


    % metadata_table, if entered
    if ~isempty(obj.metadata_table)
    
        m.metadata_table = obj.metadata_table;

        % m = average_metadata_table(m, group_by, u); % can handle single-group case (empty grouping var) or separate groups. Need to do this even if group_by = [] so we end up with 1 row

        try
            m = average_metadata_table(m, group_by, u); % can handle single-group case (empty grouping var) or separate groups. Need to do this even if group_by = [] so we end up with 1 row
        catch
            warning('Warning: Error Constructing Metadata Table')
        end

    end


    copyfield = {'images_per_session','Y','Y_names','Y_descrip',...
        'covariates','covariate_names','additional_info'};
    for i = 1:length(copyfield)
        m.(copyfield{i}) = obj.(copyfield{i});
    end

else
    if dogroup, error('Grouping only permitted for fmri_data objects...extend code if needed'), end

    m = image_vector('dat', mean(obj.dat', 1, 'omitnan')', 'volInfo', obj.volInfo, 'noverbose');
end

% Copy selected other fields
% -------------------------------
copyfield = {'source_notes','history', 'removed_voxels'};
for i = 1:length(copyfield)
    m.(copyfield{i}) = obj.(copyfield{i});
end

if ~iscell(m.history) || isempty(m.history), m.history = {m.history}; end
m.history{end + 1} = 'Averaged images';
if dogroup, m.history{end} = 'Averaged images by group'; end

% Selective average by group of other fields
% -------------------------------
u = [];
m.Y = average_var_by_group(m.Y, group_by, u);                       % can handle single-group case (empty grouping var) or separate groups
m.covariates = average_var_by_group(m.covariates, group_by, u);     % can handle single-group case (empty grouping var) or separate groups
m.X = average_var_by_group(m.X, group_by, u);     % can handle single-group case (empty grouping var) or separate groups

% Recast if needed
% -------------------------------
if isa(obj, 'fmri_data_st')
    m = fmri_data_st(m);
end

% Not completed for statistic_image
% if isa(obj, 'statistic_image')
%     m = statistic_image(m);
% end

% Other items
% -------------------------------
% convert NaNs in means back to 0s for compatibility
m.dat(isnan(m.dat)) = 0;

% Plots and other output
% -------------------------------

if doplot || doorth
    orthviews(m);
end

if doplot || dohist
    create_figure('image_histogram');
    histogram(m);
end

if ~isempty(fname)
    fullp = fullfile(fpath, fname);
    fprintf('Writing mean image to disk')

    m.filename = fname;
    m.fullpath = fullp;

    write(m);
end

% Optional outputs:
% Row means, Column means, double-centered object

if nargout > 1

    imagemeans = nanmean(obj.dat);
    varargout{1} = imagemeans;

end

if nargout > 2
    voxelmeans = nanmean(obj.dat, 2);
    varargout{2} = voxelmeans;
end


end % function





function obj = average_metadata_table(obj, group_by, u)

% average values in entire table if group_by is missing
if isempty(group_by)
    group_by = ones(size(obj.metadata_table, 1), 1);
    u = 1;
end

newtable = obj.metadata_table(1, :); % placeholder

for i = 1:length(u)
    % for each group

    t = obj.metadata_table(group_by == u(i), :); % this group/participant
    newt = t(1, :); % placeholder

    for j = 1:size(t, 2) % for each column

        tcolumn = table2array(t(:, j));

        if isnumeric(tcolumn)
            newt(1, j) = table(nanmean(tcolumn));
        elseif iscell(tcolumn) | ischar(tcolumn)
            newt(1, j) = table({char(mode(categorical(cellstr(tcolumn))))});
        else
             warning(['No policy for averaging datatype ''' class(tcolumn(1)) ''' in metadata_table. Dropping column ''' t.Properties.VariableNames{j} '''']);
             
             % newt{1, j} = [];
             % newt(1, j) = NaN;
             newt(1, j) = table(NaN);
        end
    end

    newtable(i, :) = newt;

end % groups/participants

obj.metadata_table = newtable;

end % function

% NOTES: from fmri_data_st, in case we want to use later for char and to
% print warnings
% if ~isempty(t)
%         tnames = t.Properties.VariableNames;
%         vars = cell(1,length(tnames));
%         for i = 1:length(tnames)
%             if isnumeric(t.(tnames{i}))
%                 vars{i} = nanmean(t.(tnames{i}));
%             elseif iscell(t.(tnames{i}))
%                 vars{i} = t.(tnames{i})(1);
%             elseif ischar(t.(tnames{i}))
%                 vars{i} = t.(tnames{i})(1,:);
%             else
%                 warning(['No policy for averaging datatype ''' class(t.(tnames{i})) ''' in metadata_table. Dropping column ''' tnames{i} '''']);
%                 vars{i} = [];
%             end
%         end


function newv = average_var_by_group(v, group_by, u)  % can handle single-group case (empty grouping var) or separate groups
% average values in table
if isempty(group_by)
    group_by = ones(size(v, 1), 1);
    u = 1;
end

newv = NaN * zeros(length(u), size(v, 2));

if isempty(v), return, end
assert(isnumeric(v))
if size(v, 1) ~= length(group_by), error('Variable and group_by vector must be the same length. Check input object fields and grouping variable'); end


for i = 1:length(u)
    % for each group

    wh = group_by == u(i); % this group/participant

    newv(i, :) = nanmean(v(wh, :), 1); 

end

end % function
