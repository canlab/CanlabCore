function [cm, kpos, kneg] = hotcool_split_colormap(nvals, clim, axis_handle, varargin)
%
% cm = hotcool_split_colormap(nvals, clim, [lowhot hihot lowcool hicool]), each is [r g b] triplet
%
% Get hotcool split colormap (around 0)
% note: colormaps with these colors match spm_orthviews_hotcool_colormap used in orthviews()
% [.4 .4 .3], [1 1 0] (hot/pos) and [0 0 1], [.4 .3 .4] (cool/neg)

% Split blue/yellow
% lowhot = [.4 .4 .3];
% hihot = [1 1 0];
% lowcool = [0 0 1];
% hicool = [.4 .3 .4];

% Default addblobs - orange/pink
hihot = [1 1 0]; % max pos, most extreme values
lowhot = [1 .4 .5]; % [.8 .3 0]; % min pos
hicool = [0 .8 .8]; % [.3 .6 .9]; % max neg
lowcool = [0 0 1]; % min neg, most extreme values

graybuffer = 20;  % for border - fade to gray in lowest k values

if ~isempty(varargin)
    if length(varargin) < 4
        error('Enter no optional arguments or 4 [r g b] color triplets');
    end
    lowhot = varargin{1};
    hihot = varargin{2};
    lowcool = varargin{3};
    hicool = varargin{4};
end

cmgray = gray(nvals);

thr = 0;    % Threshold to split colormap into two colors

% vec = linspace(clim(1), clim(2), nvals);
% pos = find(vec > thr);
% np = length(pos);

[hotcm, coolcm] = deal([]);

if graybuffer
    
    hotcm = [colormap_tor([.5 .5 .5], lowhot, 'n', graybuffer); ...
        colormap_tor(lowhot, hihot, 'n', nvals-graybuffer)];
    
else
    
    hotcm = colormap_tor(lowhot, hihot, 'n', nvals); %np);
    
end

%hotcm(1, :) = [.5 .5 .5]; % first element is gray
% end

% neg = find(vec < -thr);
% nn = length(neg);

% if nn > 0
%coolcm(end, :) = [.5 .5 .5]; % last element is gray
% end

if graybuffer
    
    coolcm = [colormap_tor(lowcool, hicool, 'n', nvals-graybuffer); ...
        colormap_tor(hicool, [.5 .5 .5], 'n', graybuffer)];
    
else
    
    coolcm = colormap_tor(lowcool, hicool, 'n', nvals); % nn);
    
end


% for interpolation issues
% 1:256 is gray-scale anat; 257:512 is buffer; 513:768 is negcm; 769:1024
% is buffer; 1025:1280 is poscm
buffercm = repmat([.5 .5 .5], nvals, 1);

% which block of nvals (usually 256) has colored values
kpos = 61;
kneg = 55;


%cm = [coolcm; buffercm; buffercm; cmgray; hotcm];
cm = [cmgray; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    buffercm; buffercm; buffercm; buffercm; coolcm; ...
    buffercm; buffercm; buffercm; buffercm; buffercm; ...
    hotcm];

% Kludgy fix for interpolation issues
% cm(128+3*256:4*256, :) = .5;
% cm(1+3*256:50+3*256, :) = .5;

% cm(128+0*256:1*256, :) = .5;
% cm(1+0*256:50+0*256, :) = .5;

% cm = [cmgray; cm];

% Change colormap(s)
if isempty(axis_handle)
    colormap(cm)
    
elseif iscell(axis_handle)
    
    for i = 1:length(axis_handle)
        colormap(axis_handle{i}, cm);
    end
    
else
    colormap(axis_handle, cm);
end

end