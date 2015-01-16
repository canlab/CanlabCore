function saveplots(fmri_dat, varargin)

% Output dir
if isempty(varargin)
    savedir = pwd;
else
savedir = varargin{1};
end

if ~exist(savedir, 'dir'), mkdir(savedir); end


% Register of possible plot names associated with this object type

plotnames = {'fmri data matrix' ...
    'Orthviews_means_by_unique_Y' ...
    'means by condition (unique Y values)' ...
    'Montage_coeff_of_var_across_conditions' ...
    'Montage_mean_across_conditions' ...
    'Orthviews_fmri_data_mean_and_std'
    };

for j = 1:length(plotnames)
    
    han = findobj('Name', plotnames{j});
    
    if isempty(han), continue, end
    
    han = han(end); % only take last instance
    if ishandle(han) && strcmp(get(han, 'Type'), 'figure')
        
        fprintf('Saving figure: %s\n', plotnames{j});
        figure(han);
        scn_export_papersetup(500);
        
        name = get(han, 'Name');
        name(name == ' ') = '_';
        
        name = fullfile(savedir, name);
        
        saveas(han, name, 'png');
        
    end
    
end

end