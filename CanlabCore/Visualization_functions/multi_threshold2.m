% ::
%
%    cl = multi_threshold(P, type, df, ['overlay', overlay file name], ['thresholds', thrs], ['large montage prompt', 0|1], ['title', titlestring], ['save images', 0|1])
%

% F contrast:  df = xSPM.df;
%
% type = 'F' or 'T'or 'none'
%


function cl = multi_threshold2(P, type, df, varargin)

    % ---------------------------------------------------
    % defaults
    % ---------------------------------------------------
    % Colors

    red = [1 0 0];
    yellow = [1 1 0];
    orange = [.9 .5 0];
    orange2 = [1 .6 .1];
    orange3 = [1 .7 .3];
    buff = [.7 .6 .4];
    dkred = [.7 .3 .1];

    blue = [0 0 1];
    ltblue = [0 .5 .9];
    ltblue2 = [.2 .7 1];
    aqua = [0 .5 .5];
    dkblue = [0 .1 .7];
    dkblue2 = [0 .2 .5];

    % low thresh colors
    colors  = {yellow orange3 orange2};
    colors2 = {blue ltblue ltblue2};


    thr = [Inf .001 .005];   % Inf = FDR, thr should be increasing p-values
    sizes = [3 5 10];     % minimum sizes
    direction = 'both';  % 'pos', 'neg', 'both'
    ovl = [];
    prompt_for_large_montages = 0;
    saving_images = 0;
    maskimg = [];       % empty, or specify image name

    
    for i=1:length(varargin)
        if(ischar(varargin{i}))
            switch(varargin{i})
                case {'ovl', 'overlay'}
                    ovl = varargin{i+1};
                case 'large montage prompt'
                    prompt_for_large_montages = varargin{i+1};
                case 'title'
                    fig_title = varargin{i+1};
                case {'saveimages' 'save images'}
                    saving_images = varargin{i+1};
                case 'thresholds'
                    thr = varargin{i+1};
            end
        end
    end


    
    % read image
    V = spm_vol(P);
    v = spm_read_vols(V);

    % mask image, if necessary
    if ~isempty(maskimg)
        disp(['Masking with: ' maskimg])
        
        % sample mask data in space of input images
        vmask = scn_map_image(maskimg, P(1,:));

        % make sure mask is 1's or 0's
        vmask = double(vm > 0);
        
        v = v .* vmask;
        
    end

    % convert to p-values
    switch type
        case 'F'
            p = 1 - fcdf(v(:), df(1), df(2));
            direction = 'pos';
        case 'T'
            v(isnan(v)) = 0;
            p = 1 - tcdf(v(:), df);
        case {'P' 'p' 'none'}
            direction = 'pos';  % for p-image
            % do nothing
            p = v(:);
            %thr = [3 2 1.2];
            colors = {[1 1 0] [.9 .5 0] [.7 0 0]};   % pos
    end

    cumwh =  false(size(p)); % define cumulative significant voxels
    cumwh2 = false(size(p)); % define cumulative significant voxels

    for i = 1:length(thr)
        % find significant voxels - only those sig at this thresh, but not
        % higher ones tested previously

        if isinf(thr(i))
            tmp = p;
            tmp(v == 0 | isnan(v)) = [];   % eliminate out-of-analysis voxels
            pt = FDR(tmp, .05);

            switch type
                case 'F'
                    disp(['Height thresh FDR: F = ' num2str(finv(1-pt, df(1), df(2))) ', p = ' num2str(pt)])
                case 'T'
                    disp(['Height thresh FDR: T = ' num2str(tinv(1-pt, df(1))) ', p = ' num2str(pt)])
            end
        else
            pt = thr(i);
        end

        if isempty(pt), pt = -Inf; end

        switch direction    % applies only for t-values!!
            case 'pos', 
                wh = logical(p <= pt);
                doboth = 0;
            case 'neg'
                wh = logical(p >= 1-pt);
                doboth = 0;
            case 'both'
                wh = logical(p <= pt);
                wh2 = logical(p >= 1-pt);
                doboth = 1;
        end

        % special for "none" - just threshold -- this would be for t-image, not
        % p-image!
        %if strcmp(type, 'none'), wh = logical(p >= pt); doboth = 0;, end

        % positive response, or neg only response
        mask{i} = zeros(V.dim(1:3));
        mask{i}(wh) = 1;
        mask{i}(cumwh) = 0; % voxels previously sig get 0

        cumwh = cumwh | wh; % update cumulative which sig

        % write an image of it
        V.fname = 'tmp_mask.img';
        spm_write_vol(V, mask{i});

        % get clusters
        cl{i} = mask2clusters(V.fname);


        % negative response (do both)
        if doboth
            mask2{i} = zeros(V.dim(1:3));
            mask2{i}(wh2) = 1;
            mask2{i}(cumwh2) = 0; % voxels previously sig get 0

            cumwh2 = cumwh2 | wh2; % update cumulative which sig


            % write an image of it
            V.fname = 'tmp_mask.img';
            spm_write_vol(V, mask2{i});

            % get clusters
            cl2{i} = mask2clusters(V.fname);
        end

        % delete is slower, but cross-platform
        delete('tmp_mask.img');
        delete('tmp_mask.hdr');
    end


    % size threshold
    for i = 1:length(cl)
        if ~isempty(cl{i})
            wh = find(cat(1, cl{i}.numVox) < sizes(i));
        else
            wh = [];
        end

        if ~isempty(wh)
            cl{i}(wh) = [];
        end
    end

    if doboth
        for i = 1:length(cl2)
            if ~isempty(cl2{i})
                wh = find(cat(1, cl2{i}.numVox) < sizes(i));
            else
                wh = [];
            end
            
            if ~isempty(wh)
                cl2{i}(wh) = [];
            end
        end
    end


    % add together, if both
    if doboth
        num_thr = length(thr);
        for i = 1:num_thr   % for legend
            if isempty(cl{i}), saveme(i) = 0; else saveme(i) = 1; end
        end
        len = sum(saveme);   % for legend

        thr = [thr thr];
        colors = [colors(1:min(num_thr, length(colors))) colors2(1:min(num_thr, length(colors)))];
        cl = [cl cl2];
    end

    % remove empties
    for i = 1:length(thr)
        if isempty(cl{i}), saveme(i) = 0; else saveme(i) = 1; end
    end
    saveme = logical(saveme);
    cl = cl(saveme);
    colors = colors(saveme);
    thr = thr(saveme);


    if isempty(cl), return; end

    % save clusters
    disp(['Saving clusters at multiple thresholds in multi_clusters.mat']);
    fprintf(1,'To re-create orthviews, load the file and execute the lines below:\n');
    fprintf(1,'e.g.,\nif ~isempty(cl{1}), cluster_orthviews(cl{1}, colors(1), ''overlay'', ovl); end\n');
    fprintf(1,'for i = 2:length(cl), if ~isempty(cl{i}), cluster_orthviews(cl{i}, colors(i), ''add''); end, end\n');
    
    save multi_clusters cl colors ovl
    
    % now display image clusters
    if ~isempty(cl{1}), cluster_orthviews(cl{1}, colors(1), 'overlay', ovl); end
    for i = 2:length(cl)
        if ~isempty(cl{i}), cluster_orthviews(cl{i}, colors(i), 'add'); end
    end


    % make legend string and legend
    for i = 1:length(cl)
        lstr = [];
        if doboth
            if i <= len
                lstr = '+ ';
            else
                lstr = '- ';
            end
        end
        if isinf(thr(i)), 
            lstr = [lstr 'p < .05 FDR'];
        else 
            lstr = [lstr 'p < ' num2str(thr(i))];
        end
        legstr{i} = lstr;
    end

    h = axes('Position', [.55 .2 .25 .2]); hold on;
    set(gca, 'FontSize', 24)
    for i = 1:length(colors)
        hh(i)=plot(0, 0, 'Color', colors{i}, 'LineWidth', 10);
    end
    axis off
    legend(legstr);


    % check for cases in which we DO NOT want montage

    if strcmp(type, 'none'), return; end

    tmp = [];
    for i = 1:length(cl)
        tmp = [tmp cat(2, cl{i}.XYZ)];
    end
    tmp = unique(tmp(3, :));

    if(prompt_for_large_montages && length(tmp) > 30)
        cont = input('More than 30 slices in montage.  Make montage figure? (1/0) ');
        if(~cont)
            return
        end

    end

    % montage
    mch = montage_clusters(ovl, cl{:}, colors);

    if(~isempty(mch))
        axial_fig_title = sprintf('%s - axial', fig_title);
        set(mch, 'Name', axial_fig_title);
        uicontrol(mch, 'Style', 'text', 'String', axial_fig_title, 'Units', 'normalized', 'Position', [0 .97 1, .03], ...
            'BackgroundColor', 'white', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized', 'FontSize', .9);
        if(saving_images)
            saveas(mch, [strrep(axial_fig_title, '/', '') '.png']);
        end
    end
    cl{1}(1).thr = thr;
    cl{1}(1).colors = colors;



    % montage medial slices
    tmp = unique(tmp(1, :));
    if(prompt_for_large_montages && length(tmp) > 30)
        cont = input('More than 30 slices in medial montage.  Make montage figure? (1/0) ');
        if ~cont
            return
        end
    end

    
    mcmh = montage_clusters_medial(ovl, cl{:}, colors);

    if(~isempty(mcmh))
        medial_fig_title = sprintf('%s - medial', fig_title);
        set(mcmh, 'Name', medial_fig_title);
        uicontrol(mcmh, 'Style', 'text', 'String', medial_fig_title, 'Units', 'normalized', 'Position', [0 .97 1, .03], ...
            'BackgroundColor', 'white', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized', 'FontSize', .9);
        if(saving_images)
            saveas(mcmh, [strrep(medial_fig_title, '/', '') '.png']);
        end
    end
    cl{1}(1).thr = thr;
    cl{1}(1).colors = colors;
end


