function cmap = canlab_perceptual_colormap(name, n)
% Return a perceptual (perceptually-uniform) colormap as an n x 3 RGB matrix.
%
% Provides the popular sequential perceptual colormaps used by NiiVue / matplotlib
% (viridis, inferno, magma, plasma) plus MATLAB's built-in turbo and parula, as a
% single keyword-resolved helper. Intended for use as a CONTINUOUS colormap in the
% CANlab display pipeline (e.g. set_colormap(o2, 'colormap', canlab_perceptual_colormap('viridis'))),
% where data values are mapped continuously through the LUT.
%
% :Usage:
% ::
%
%     cmap = canlab_perceptual_colormap(name, [n])
%
% :Inputs:
%
%   **name:**
%        Colormap name (case-insensitive): 'viridis', 'inferno', 'magma',
%        'plasma', 'turbo', or 'parula'.
%
% :Optional Inputs:
%
%   **n:**
%        Number of rows (colours). Default: 256.
%
% :Outputs:
%
%   **cmap:**
%        n x 3 RGB matrix in [0, 1].
%
% :Examples:
% ::
%
%     cmap = canlab_perceptual_colormap('inferno');
%     figure; colormap(cmap); colorbar;
%
%     % In a managed display:
%     o2 = montage(load_image_set('emotionreg'));            %#ok<*NASGU>
%     o2 = set_colormap(o2, 'colormap', canlab_perceptual_colormap('viridis'));
%
% :See also:
%   - canlab_colormap (the value->colour module that consumes these), set_colormap, turbo, parula
%
% ..
%    2026 visualization overhaul — central colour pipeline
% ..

if nargin < 2 || isempty(n), n = 256; end
name = lower(strtrim(name));

% MATLAB built-ins
switch name
    case 'turbo',  cmap = turbo(n);  return
    case 'parula', cmap = parula(n); return
end

% matplotlib sequential maps as anchor points (RGB 0-255), interpolated below.
switch name
    case 'viridis'
        anchor = [ 68   1  84; 72  36 117; 65  68 135; 53  95 141; 42 120 142; ...
                   33 145 140; 34 168 132; 68 191 112; 122 209  81; 189 223  38; 253 231  37];
    case 'inferno'
        anchor = [  0   0   4; 22  11  57; 58   9  99; 96  19 110; 133  33 107; ...
                  169  46  94; 203 65  73; 230 93  47; 247 131  17; 252 181  24; 252 255 164];
    case 'magma'
        anchor = [  0   0   4; 22  11  57; 59   9 103; 99  21 110; 139  35 105; ...
                  179  48  94; 219 67  77; 246 105 75; 253 155 108; 254 205 144; 252 253 191];
    case 'plasma'
        anchor = [ 13   8 135; 75   3 161; 125  3 168; 168 34 150; 203  70 121; ...
                  229 107  93; 248 148 65; 253 195 40; 240 249  33; 240 249  33; 240 249  33];
    otherwise
        error('canlab_perceptual_colormap:unknownName', ...
            ['Unknown colormap ''%s''. Options: viridis, inferno, magma, plasma, ' ...
             'turbo, parula.'], name);
end

anchor = anchor / 255;
x  = linspace(0, 1, size(anchor, 1));
xi = linspace(0, 1, n);
cmap = [interp1(x, anchor(:, 1), xi, 'pchip')', ...
        interp1(x, anchor(:, 2), xi, 'pchip')', ...
        interp1(x, anchor(:, 3), xi, 'pchip')'];
cmap = min(max(cmap, 0), 1);

end
