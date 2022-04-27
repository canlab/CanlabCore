function saveplots(fmri_model, varargin)

if isempty(varargin) % Output dir
    savedir = pwd;
else
savedir = varargin{1};
end

if ~exist(savedir, 'dir'), mkdir(savedir); end

% Register of possible plot names associated with this object type

plotnames = {'basis sets' ...
    'design matrix image' ...
    'design matrix' ...
    };

for j = 1:length(plotnames)
    
    han = findobj('Name', plotnames{j});
    if ~isempty(han) && ishandle(han) && strcmp(get(han, 'Type'), 'figure')
        
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
