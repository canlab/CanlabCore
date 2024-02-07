function [obj_separated_matched, levels_table] = unstack_by_condition(obj, condition_names, sid_var_name)
% unstack image into separate objects for each condition based on 1 or more metadata_table variables
%
% [obj_separated_matched, levels_table] = unstack_by_condition(obj, condition_names, sid_var_name)
%
% This object method uses the metadata_table field of an image_vector object (e.g., fmri_data object) 
% to separate images belonging to different groups or conditions. It can be used along with the mean() 
% method to create average images for subgroups (e.g., participants) separated by condition. This turns 
% a "long format" object, for example one with multiple trials per condition per participant, into a 
% "wide format" dataset with one object per condition, with the order of images matched on participant 
% across conditions. This "wide format" is suitable for estimating within-person contrasts and doing 
% related statistical tests.
%
% e.g., 
% condition_names = {'stimLvl' 'heat' 'reg'};  % condition names
% sid_var_name = 'subject_id';               % must be numeric 
%
% This will create a cell array with one cell per unique combination of the
% levels of variables in condition_names
% e.g., if there are 3 levels of reg, 2 levels of stimLvl, and 1 level of heat
% will create 3 x 2 = 6 separate image objects, each in a cell
%
% To handle missing data, assuming conditions are within-subject:
% Match images in each object on a subject (or other grouping) variable
% * All variables must be numeric *
% This is an alpha-version code stub; please finish me using
% documentation_template.m

% Programmers' notes:
% Created by Tor Wager, Jan 2024

% matched rows on subject

dealtable = obj.metadata_table(:, condition_names);
levels_table = unique(dealtable, 'rows');

obj_separated = cell(1, size(levels_table, 1));
for i = 1:length(obj_separated)

    wh = ismember(dealtable, levels_table(i, :), 'rows');
    obj_separated{i} = get_wh_image(obj, find(wh));


end

% make sure we have matched rows for all subjects
% sid must be numeric!
for i = 1:length(obj_separated)
    sid{i} = obj_separated{i}.metadata_table.(sid_var_name); 
end
sid = unique(cat(1, sid{:}));

obj_separated_matched = cell(1, size(levels_table, 1));
all_sids = [];

for i = 1:length(obj_separated)
    [~, ~, indx] = intersect(sid, obj_separated{i}.metadata_table.(sid_var_name), 'stable');
    obj_separated_matched{i} = get_wh_image(obj_separated{i}, indx);

    % check
    all_sids(:, i) = obj_separated_matched{i}.metadata_table.(sid_var_name);
end

if any(any(diff(all_sids'))), error('uh oh, subjects do not match'); end

end % function

