function spm_orthviews_white_background()
% :Usage:
% ::
%
%    spm_orthviews_white_background
%
% Set up orthviews to draw with a white background and softened grayscale
% brain
%
% No input arguments.  Updates spm's global variable st
% 
% Uses SPM_OV_BLACK2WHITE.M

spm('Defaults','fmri')

global st

for i = 1:24 % for all vols
    if ~isstruct(st.vols{i}), continue, end
    st.vols{i}.black2white = 1;
end

bwexist = strfind(st.plugins, 'black2white');
bwexist = any(cat(2, bwexist{:}));
if ~bwexist
    st.plugins{end+1} = 'black2white';
end

f1 = findobj('Type', 'figure', 'Tag', 'Graphics');
figure(f1); set(f1, 'Color', 'w');

spm_orthviews('redraw')

end


% for i = 1:24 % for all vols
%
%     if ~isstruct(st.vols{i}), continue, end
%
%     for a = 1:3
%
%         cd = get(st.vols{i}.ax{a}.d, 'CData');
%
%         if size(cd, 3) == 3 % is true color matrix
%
%             % turn black to white
%
%             isblack = double(all(cd == 0, 3));
%             cd = cd + repmat(isblack, [1 1 3]);
%             set(st.vols{i}.ax{a}.d, 'CData', cd);
%
%             % lighten the dark edges
%             for k = [.01:.01:.35] % for each of these darkness values
%                 w = 1 - k; % weights, for wtd average of orig and inverted, 1 is all inverted
%                 % (bilinear weighting)
%
%              isdark = all(cd < k, 3);
%              for j = 1:3
%                  %invert
%                  colorslice = cd(:, :, j);
%                  colorslice(isdark) = w .* (1 - colorslice(isdark)) + (1 - w) .* colorslice(isdark);
%
%                  cd(:, :, j) = colorslice;
%              end
%
%             end % the weight value loop
%
%             set(st.vols{i}.ax{a}.d, 'CData', cd);
%
%         end
%
%     end
%
% end
%
