function cm = spm_orthviews_hotcool_colormap(t, thr)
    % cm = spm_orthviews_hotcool_colormap(t, thr)
    %
    % Create split-level colormap with hot colors for positive values and
    % cool colors for negative ones
    % Apply this colormap to the spm Graphics window (or any figure with the
    % tag 'Graphics')
    %
    % Inputs:
    %   t:  range of input statistic values (doesn't have to be t-values)
    %   thr: threshold more extreme than which (+/-) colors should be used
    %       enter a positive value
    %
    % Designed to work on blobs added using spm_orthviews 'addblobs' feature
    %
    % Tor Wager, April 2007
    %
    % Examples for using 'clusters' in torlab format:
    % ----------------------------------------------------------------------
    % Threshold an image to get clusters:
    % [dat, volInfo, cl] = iimg_threshold('test_statistic.img', 'thr', 3.1440, 'k', 20);
    %
    % Display the clusters using spm_orthviews and apply the new color map:
    % cluster_orthviews(cl); cm = spm_orthviews_hotcool_colormap(cat(2,cl.Z), 3);


    myfig = findobj('Tag','Graphics');
    axishandles = findobj(myfig, 'Type', 'Axes');
    %cm = get(myfig,'Colormap');


    % range of input statistic values (called 't' here, but doesn't have to be)
    mx = max([eps max(t)]);
    mn = min([0 min(t)]);

    % indices in color map for range of data (specified in spm_orthviews)
    cvals = [1:64]' + 64;

    % t-vals corresponding to c-vals
    tvals = linspace(mn, mx, length(cvals));


    % make new split colormap: hot colors for pos, cool colors for neg

    newcm = NaN .* zeros(64, 3);  %initialize


    pos = find(tvals > thr);
    n = length(pos);

    if n > 0
        hotcm = hotmap(n);

        % gray-yellow
        %hotcm = [linspace(0.5, 1, n)' linspace(0.5, 1, n)' linspace(.5, 0, n)' ];


        newcm(end - n + 1 : end,:) = hotcm;
    end

    neg = find(tvals < -thr);
    n = length(neg);

    %coolcm = cool(n);

    if n > 0
        coolcm = [linspace(0, 0.5, n)' linspace(0, 0.7, n)' linspace(1, .5, n)' ];

        newcm(1:n,:) = coolcm;
    end


    % specify split color-map of 128 values; first 64 are gray, last 64 have
    % colors for activations; this is specified by spm_orthviews

    cm = [gray(64) ; newcm];

    for i = 1:length(axishandles)
        colormap(axishandles(i), cm);
    end

end



function hotcm = hotmap(n)

    offset = 3;
    hotcm = hot(n + offset);
    hotcm = hotcm(offset+1:end,:);

    hotcm = hotcm + greyfade(n, offset+2);  % blend hot map with fade-out grayscale map
    hotcm(hotcm > 1) = 1;

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
