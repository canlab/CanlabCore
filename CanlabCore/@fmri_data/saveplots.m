function saveplots(fmri_dat, varargin)
% saveplots Save the standard fmri_data diagnostic figures to disk as PNGs.
%
% :Usage:
% ::
%
%     saveplots(fmri_dat)
%     saveplots(fmri_dat, savedir)
%
% Iterates over a registered list of figure-window names produced by the
% standard fmri_data plotting / QC routines (e.g., plot, orthviews,
% montage) and writes each one that currently exists to a PNG file in
% the requested directory. Figures whose registered names don't currently
% exist on screen are silently skipped.
%
% :Inputs:
%
%   **fmri_dat:**
%        fmri_data object. Currently used only to dispatch the method;
%        the figures saved are the ones currently open on screen, not
%        regenerated from the object.
%
%   **savedir:** *(optional)*
%        Output directory. Default: pwd. Created with mkdir if
%        it does not already exist.
%
% :Outputs:
%
%   None. PNG files are written into savedir with names derived from
%   the figure window names (spaces replaced by underscores).
%
% :Registered figure names:
%
%   ::
%
%       'fmri data matrix'
%       'Orthviews_means_by_unique_Y'
%       'means by condition (unique Y values)'
%       'Montage_coeff_of_var_across_conditions'
%       'Montage_mean_across_conditions'
%       'Orthviews_fmri_data_mean_and_std'
%
% :See also:
%   - plot (creates several of these figures), orthviews, montage
%   - scn_export_papersetup (page-setup helper used internally)

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