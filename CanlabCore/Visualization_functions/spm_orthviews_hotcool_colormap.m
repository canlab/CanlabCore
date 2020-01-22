function cm = spm_orthviews_hotcool_colormap(t, thr)
% :Usage:
% ::
%
%    cm = spm_orthviews_hotcool_colormap(t, thr)
%
% Create split-level colormap with hot colors for positive values and
% cool colors for negative ones
%
% Apply this colormap to the spm Graphics window (or any figure with the
% tag 'Graphics')
%
% :Inputs:
%
%   **t:**
%        range of input statistic values (doesn't have to be t-values)
%        *NOTE: in the 2020 update, this is no longer used, so you can
%        enter []
%
%   **tthr:**
%        threshold more extreme than which (+/-) colors should be used
%       enter a positive value
%
% Designed to work on blobs added using spm_orthviews 'addblobs' feature
%
% :Examples: for using 'clusters' in torlab format
% ::
%
%    % Threshold an image to get clusters:
%    [dat, volInfo, cl] = iimg_threshold('test_statistic.img', 'thr', 3.1440, 'k', 20);
%
%    % Display the clusters using spm_orthviews and apply the new color map:
%    cluster_orthviews(cl); cm = spm_orthviews_hotcool_colormap(cat(2,cl.Z), 3);
%
% ..
%    Tor Wager, April 2007
%    Updated color map Jan 2020 - set zero point based on range shown in
%    SPM st
% 
% ..

    global st
    if isempty(st), reset_st; end
    
    %myfig = findobj('Tag','Graphics');
    axishandles = findobj(st.fig, 'Type', 'Axes');

    % range of input statistic values (called 't' here, but doesn't have to be)
    mx = max([eps max(t)]);
    mn = min([0 min(t)]);

    % We need to find the zero-point of the data range displayed
    % ----------------------------------------------------------
    % indices in color map for range of data (specified in spm_orthviews)
    cvals = [1:64]' + 64;

    % data range shown (specified in spm_orthviews)
    % If not defined, exit - we may have a different kind of image, e.g.,
    % colored blobs with no colormap
    if ~isfield(st.vols{1}.blobs{1}, 'cbar') || isempty(st.vols{1}.blobs{1}.cbar)
        return
    end
    
    interval = st.vols{1}.blobs{1}.cbar.YLim;
    vec = linspace(interval(1), interval(2), 64);
    
    % t-vals corresponding to c-vals
    %tvals = linspace(mn, mx, length(cvals));

    % make new split colormap: hot colors for pos, cool colors for neg

    newcm = NaN .* zeros(64, 3);  %initialize

    %pos = find(tvals > thr);       % Changed Jan 2020
    pos = find(vec > thr);
    
    n = length(pos);

    if n > 0
        hotcm = hotmap(n);
        newcm(end - n + 1 : end,:) = hotcm;
    end

    % neg = find(tvals < -thr);
    neg = find(vec < -thr);
    n = length(neg);

    %coolcm = cool(n);

    if n > 0
        
        %coolcm = [linspace(0, 0.5, n)' linspace(0, 0.7, n)' linspace(1, .5, n)' ];
        % replaced Jan 2020
        coolcm = colormap_tor([0 0 1], [.4 .3 .4], 'n', n);
        
        newcm(1:n,:) = coolcm;
    end

    % specify split color-map of 128 values; first 64 are gray, last 64 have
    % colors for activations; this is specified by spm_orthviews

    cm = [gray(64) ; newcm];

    spm_figure_canlab('Colormap',cm)
    
%     for i = 1:length(axishandles)
%         colormap(axishandles(i), cm);
%     end

% Set Transparency map
% This will not work because SPM creates single-layer image for anatomy and
% colors. So this will be left for a later date.
%
% hh = findobj(st.fig, 'Type', 'Image');
% for i = 1:length(hh)
% set(hh(i), 'AlphaDataMapping', 'scaled', 'AlphaData', abs(hh(i).CData));
% end

end % function



function hotcm = hotmap(n)

%     offset = 3;
%     hotcm = hot(n + offset);
%     hotcm = hotcm(offset+1:end,:);
% 
%     hotcm = hotcm + greyfade(n, offset+2);  % blend hot map with fade-out grayscale map
%     hotcm(hotcm > 1) = 1;

% Replaced jan 2020:
    hotcm = colormap_tor([.4 .4 .3], [1 1 0], 'n', n);
    
end


function grayfadeout = greyfade(n, grayels)

    %     div = 5;
    %     grayels = floor(n ./ div);
    %
    %     grayels = 3;

    grayfadeout = linspace(.5, 0, grayels)';
    grayfadeout = repmat(grayfadeout, 1, 3);
    grayfadeout = [grayfadeout; zeros(n - grayels, 3)];

    grayfadeout = grayfadeout(1:n, :);
end




%     coolcm = [linspace(0, 0, n)' linspace(0, 0, n)' linspace(1, 0, n)' ];
%     coolcm = coolcm + flipud(greyfade(n));  % blend hot map with fade-out grayscale map
%     coolcm(coolcm > 1) = 1;
