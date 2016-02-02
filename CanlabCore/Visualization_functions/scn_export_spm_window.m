% :Usage:
% ::
%
%    scn_export_spm_window(meth, [savefilename] or [overlay])
%
% Modes:
%   - scn_export_spm_window('setup')
%   - scn_export_spm_window('setup', overlayimgname)
%
% Set up SPM figure window paper size, etc. for saving
%
% The optional argument sets the spm_orthviews display range based on your
% image
%
% scn_export_spm_window('save', 'myfile')
% saves SPM window to myfile.png
%
% Run setup first, then save.
%
% ..
%    tor wager
%    aug 3, 06
% ..

function scn_export_spm_window(meth, savefilename)
    fh = findobj('Tag', 'Graphics');
    if isempty(fh) || ~ishandle(fh)
        disp('No SPM figure window.  Using current figure.');
        fh = gcf;
    end

    switch meth
        case 'setup'
            if nargin > 1 && ~isempty(savefilename)
                % set display window on orthviews based on input image
                v = spm_read_vols(spm_vol(savefilename));
                v = v(:);
                ub = mean(v)+4*std(v);
                %spm_orthviews('Window', 1, [0 ub])
            end

            % set colors
            set(fh, 'InvertHardCopy', 'off');
            set(fh, 'Color', 'black');

            scn_export_papersetup;

            % turn off crosshairs
            spm_orthviews('Xhairs', 'off');

            
            % set axis color and ticks of colorbars
            global st
            for i = 1:length(st.vols)

                % set color bar colors, if color bar exists
                if ~isempty(st.vols{i}) && ~isempty(st.vols{i}.blobs) && isfield(st.vols{i}.blobs{1}, 'cbar') && ~isempty(st.vols{i}.blobs{1}.cbar)
                    % colors
                    set(st.vols{i}.blobs{1}.cbar, 'YColor', 'g', 'XColor', 'g', 'FontSize', 16)
                    
                    % set axis ticks of colorbars
                    minmax = get(st.vols{i}.blobs{1}.cbar, 'YLim');
                    set(st.vols{i}.blobs{1}.cbar, 'YTick', linspace(minmax(1), minmax(2), 5));

                end
            end

        case 'save'
            % now we're ready to save
            saveas(fh, savefilename, 'png');

        otherwise
            error('Meth should be ''setup'' or ''save''.')
    end
end




function scn_export_papersetup
    % make sure that min size of fig. is 400 pixels
    sz = get(gcf, 'Position');
    sz = sz(3:4);   % width and height
    szratio = max(sz) ./ min(sz);   % max to min size
    wh = find(sz == min(sz));
    sz(wh) = 400;
    wh = find(sz == max(sz));
    sz(wh) = 400 * szratio;

    % set paper size using current screen
    set(gcf, 'PaperUnits', 'points', 'PaperPosition', [0 0 sz], 'PaperType', 'usletter');
end

