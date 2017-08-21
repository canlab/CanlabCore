function mch = montage_clusters(ovl, clusters, varargin)
% :Usage:
% ::
%
%    fig_handle = montage_clusters(ovl, clusters, varargin)
%
% :varargin: (in any order) =
%   a) additional clusters structures
%   b) cell array of colors (text format), must be ROW vector {'r' 'g'} etc...
%       if length of color string is longer than number of clusters inputs,
%       additional colors are interpreted as '2-intersection' and 'all-intersection'
%       colors, in that order.  This overrides single string argument (c, below) for
%       color input
%   c) single string color argument for overlaps (intersections)
%       plots intersections of ANY TWO clusters right now.
%       also color for plotting points, if entered.
%       use 'nooverlap' as an input argument to suppress this.
%   d) [n x 3] matrix of points to plot on map
%   e) text labels for points, must be cell array in COLUMN vector
%   f) single number, 1/0 for whether to plot overlapping coordinates in overlap colors
%       default is 1.
%   g) color limit vector [min max] indicates color mapping for blobs
%   rather than
%       solid colors.  This will do hot/cool mapping
%
% Intersections of 2 colors are magenta, and ALL colors are yellow
% unless otherwise specified
% try:
% ::
%
%    CLU = clusters2clu(clusters);
%    spm_orthviews('AddColouredBlobs', 1, CLU.XYZ, CLU.Z, CLU.M, [1 0 0])
%
% ..
%    Tor Wager
%
%    Edited: Jan 2010, to add functionality to display underlay image only if
%    no clusters are entered (with 9 mm slice spacing)
% ..


    % ..
    %    defaults
    % ..
    
    if ischar(ovl)
        [path, name, ext] = fileparts(ovl);
        if(isempty(path))
            ovl = which(ovl);
        end
    end

    if isempty(ovl), ovl = which('spm2_single_subj_T1_scalped.img'); end
    mch = [];
    myc = {'r' 'b' 'g' 'c' 'm' 'w'};
    customcolors = 0;
    bcolor = 'm';
    acolor = 'y';
    XYZmm_both = [];
    clindx = 1;
    XYZpts = [];
    myplab = [];
    cl = [];
    plotovl = 1;
    doverb = 1;
    colorlimits = [];

    % recursive operation if clusters is a cell array
    % indicating multiple clusters structures
    %if iscell(clusters),
    %    for i = 1:length(clusters),
    %        montage_clusters(ovl, clusters{i}, varargin)
    %    end
    %    return
    %end
    
    if nargin < 2, clusters = []; end
    
    if isempty(clusters)
        cl{1} = [];
        disp('montage_clusters: no activation blobs; displaying underlay image.');
        XYZmm = zeros(3, 16);
        XYZmm(3, :) =  -60:9:80;
    else
        % we have clusters

        cl{1} = clusters;
        clindx = length(cl)+1;

        XYZmm = cat(2, clusters.XYZmm);
        XYZmm_all = XYZmm;
    end

    % ----------------------------------------
    % * process input arguments
    % ----------------------------------------

    %why not turn this into a switch statement? - SS 5/20/10
    if length(varargin) > 0
        for i = 1:length(varargin)
            if isempty(varargin{i})
                % ignore it
            elseif isstruct(varargin{i}) || isa(varargin{i},'region')
                % struct inputs interpreted as clusters variables
                cl{clindx} = varargin{i};
                clXYZmm = cat(2, cl{clindx}.XYZmm);
                clindx = clindx + 1;
                a = XYZmm'; b = clXYZmm';
                XYZmm_both = [XYZmm_both intersect(a, b, 'rows')'];
                if ~isempty(XYZmm_both), XYZmm_all = intersect(XYZmm_all', b, 'rows')'; else XYZmm_all = []; end
                XYZmm = [XYZmm clXYZmm];

            elseif ischar(varargin{i})
                % string arguments interpreted as overlap colors
                if strcmp(varargin{i}, 'nooverlap'), plotovl = 0;
                else
                    bcolor = varargin{i};
                end
            elseif iscell(varargin{i})
                % cell array vectors with one row interpreted as color inputs
                if size(varargin{i}, 1) == 1
                    customcolors = 1;
                    myc = varargin{i};
                elseif size(varargin{i}, 2) == 1
                    % cell array vectors with one column interpreted as point coordinates to plot
                    myplab = varargin{i};
                else
                    error('cell array must be row (for colors) or column (for text labels) vector.')
                end
            elseif prod(size(varargin{i})) == 1
                % single integers interpreted as 'do overlap plot' flag, can be 1 or 0
                % default is 1
                plotovl = varargin{i};
            elseif any(size(varargin{i}) == 3)
                % any 3-vector matrix interpreted as list of points to plot on figure
                XYZpts = varargin{i};
                %XYZmm = [XYZmm XYZpts];
                % then it's a list of points to plot
            elseif size(varargin{i}, 1) == 1 && size(varargin{i}, 2) == 2
                % 1 x 2 double vector, means we want height-mapped colors instead of solid
                colorlimits = varargin{i};
            else
                error('Unknown input argument type.')
            end
        end
    end

    % get rid of empty clusters
    allempty = 0;
    for i =1:length(cl), wh(i) = isempty(cl{i}); end; wh = find(wh);
    if length(wh) == length(cl), allempty = 1; plotovl = 0; end
    cl(wh) = []; myc(wh) = [];

    % do not plot overlaps if only one clusters struct entered
    if length(cl) == 1, plotovl = 0; end

    if customcolors
        if length(myc) > length(cl)+1, acolor = myc{length(cl) + 2}; end
        if length(myc) > length(cl), bcolor = myc{length(cl) + 1}; end
    end

    if doverb
        disp([num2str(length(cl)) ' Clusters found.'])
        if plotovl
            disp(['Two cluster overlap: ' bcolor ', found ' num2str(length(XYZmm_both)) ' voxels.'])
            disp(['All cluster overlap: ' acolor ', found ' num2str(length(XYZmm_all)) ' voxels.'])
        else
            disp('No overlap plotting.')
        end
    end

    % ----------------------------------------
    % * overlay image
    % ----------------------------------------
    try
        V = spm_vol(ovl);
        oimg = spm_read_vols(V);
    catch
        fprintf(1, 'Could not use spm_vol or spm_read_vols. Either image "%s" is missing or you don''t have SPM.', ovl);
        [oimg, hdr, h] = readim2(ovl);
        V.mat = diag([hdr.xsize hdr.ysize hdr.zsize 1]); V.mat(:, 4) = [hdr.origin(1:3); hdr.SPM_scale];
        orig = (hdr.origin(1:3) - [hdr.xdim hdr.ydim hdr.zdim]') .* [hdr.xsize hdr.ysize hdr.zsize]';
        V.mat(1:3, 4) = orig;
        V.mat = [
            2     0     0   -92
            0     2     0  -128
            0     0     2   -74
            0     0     0     1];
    end

    %V.mat = scn_mat_conform(V.mat);
    V.M = V.mat;
    
    textx = size(oimg, 1) - 50;
    texty = 6; %size(oimg, 2) - 6;

    %[array, hdr, h, whichslices, rows, cols, figh] = readim2(ovl, 'p');

    
    % how many slices, which ones

    XYZ = mm2voxel(XYZmm, V, 1)';    % 1 no re-ordering, allows repeats             %2 THIS RE-ORDERS VOXELS, BUT IS FAST, AND ORDER SHOULDN'T MATTER HERE.
    whsl = unique(XYZ(3, :));
        
    whsl(whsl <=0) = [];
    whsl(whsl >= size(oimg, 3)) = [];

    nsl = length(whsl) + 1;
    rc = ceil(sqrt(nsl));
    h = [];

    mch = create_figure('SCNlab_Montage');
    colormap gray;
    set(mch, 'Color', [1 1 1], 'MenuBar', 'none');
%     ssize = get(0, 'screensize');
%     set(mch, 'Position', ssize ./ [0.0303    0.0128    2.0447    1.2577]);
    
% colormap works to make white background for continuous-valued images, but
% does not work well for binary images
cc = colormap(gray); cc(1:10, :) = repmat([1 1 1], 10, 1);    % [0 0 .3] for dark blue
colormap(cc)

if length(unique(oimg(:))) < 10  % binary or other strange image type
    cm = colormap(gray); cm = cm(length(cm):-1:1, :); % invert colormap
    colormap(cm)
end

    display_underlay(mch, oimg, whsl, rc, V, textx, texty), drawnow

    if allempty, return, end
    
    % ----------------------------------------
    % * plot first cluster structure
    % ----------------------------------------

    index = plot_cluster(mch, cl{1}, oimg, rc, V, whsl, myc{1}, textx, texty, 1, size(oimg), colorlimits);
    if ~isempty(XYZpts), ph = plot_points(mch, XYZpts, rc, V, whsl, bcolor, myplab); end
    % plot points for single-voxel clusters
    if isempty(colorlimits)
        for i = 1:length(clusters)
            try
                if clusters(i).numVox == 1
                    plot_points(mch, clusters(i).XYZmm, rc, V, whsl, myc{1}, myplab);
                end
            catch
            end
        end
    end

    % ----------------------------------------
    % * plot additional cluster structures
    % ----------------------------------------

    if length(cl) > 1
        for i = 2:length(cl)
            index = plot_cluster(mch, cl{i}, [], rc, V, whsl, myc{i}, textx, texty, i, size(oimg), colorlimits);
            if isempty(colorlimits)
                for j = 1:length(cl{i}),
                    try 
                        if cl{i}(j).numVox == 1, plot_points(mch, cl{i}(j).XYZmm, rc, V, whsl, myc{i}, myplab);end 
                    catch
                    end
                end
            end
        end
    end


    % ----------------------------------------
    % * plot overlap areas
    % ----------------------------------------

    if plotovl

        if ~isempty(XYZmm_both)
            bcl.XYZmm = XYZmm_both; bcl.M = cl{1}(1).M;
            plot_cluster(mch, bcl, [], rc, V, whsl, bcolor, textx, texty, i+1, size(oimg), colorlimits, 1);
        end

        if ~isempty(XYZmm_all)
            bcl.XYZmm = XYZmm_all; bcl.M = cl{1}(1).M;
            plot_cluster(mch, bcl, [], rc, V, whsl, acolor, textx, texty, i+1, size(oimg), colorlimits, 1);
        end

    end

    % fix bug at end - replot XYZ slice text
    %for z = whsl
    %    subplot(rc, rc, index);
    %    zmm = voxel2mm([1 1 z]', V.mat);
    %    text(textx, texty, ['z = ' num2str(zmm(3))], 'Color', 'w')
    %    zi(z) = zmm(3);
    %end
    %zi(zi > 0)

    try
        enlarge_axes(mch)
    catch
        disp('Error enlarging axes. Skipping.')
    end

end




% ----------------------------------------
%
% * Sub-functions
%
% ----------------------------------------

function index = plot_cluster(mch, clusters, oimg, rc, V, whsl, myc, textx, texty, clind, odims, colorlimits, varargin)
    % varargin suppresses end text

    if ~isfield(clusters, 'Z'), for i = 1:length(clusters), clusters(i).Z = ones(1, size(clusters(i).XYZmm, 2));end, end

    % take only first row of Z scores
    for i = 1:length(clusters), clusters(i).Z = clusters(i).Z(1, :);end


    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    % Set up
    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    XYZmm = cat(2, clusters.XYZmm);
    XYZ = mm2voxel(XYZmm, V, 1)';     % no unique suppression

    % get relative voxel sizes for fill
    xs = diag(clusters(1).M(1:3, 1:3));
    xs = xs ./ diag(V.mat(1:3, 1:3));
    % add 25% to make sure no gaps in plot
    xs = xs + .25*xs;

    if ~isempty(colorlimits)

        if (colorlimits(1) < 1 && colorlimits(2) < 1) || all(colorlimits - mean(colorlimits)) < eps

            % -------------------------------------------------------------
            % define 4 split color maps - red/yellow, green/blue
            % -------------------------------------------------------------
            colorlimits = 0;  % define maps within function; otherwise, use limits entered
        end


        % -------------------------------------------------------------
        % define color maps - biscale hot/cool
        % -------------------------------------------------------------

        % color map - hot
        % --------------------------------------------
        h1 = linspace(.7, 1, 250)';
        h2 = linspace(0, .8, 175)';
        h2 = [h2; linspace(.8, 1, 75)'];

        h3 = zeros(size(h1));
        h = [h1 h2 h3];

        % color map - winter
        % --------------------------------------------
        h1 = linspace(0, .5, 250)';
        h2 = linspace(1, .7, 250)';
        h3 = zeros(size(h1));
        hc = [h3 h1 h2];

        % 	% color map - hot
        % 	% --------------------------------------------
        % 	h1 = (0:1/99:1)';
        % 	h2 = ones(size(h1));
        % 	h3 = zeros(size(h1));
        % 	h = [h1 h3 h3; h2 h1 h3; h2 h2 h1];
        % 	h(1:50, :) = []; % take only red values to start
        % 	% in new matlab: h = colormap(hot(300));
        %
        % 	% color map - winter
        % 	% --------------------------------------------
        % 	h1 = (0:1/249:1)';
        % 	h2 = (1:-1/(249*2):.5)';
        % 	h3 = zeros(size(h1));
        % 	hc = [h3 h1 h2];

        % -------------------------------------------------------------
        % determine overall z-score range
        % -------------------------------------------------------------

        allz = cat(2, clusters.Z);
        zrange = allz;
        tmp = zrange(zrange > 0);
        tmpc = zrange(zrange < 0);

        if ~isempty(tmp)
            if colorlimits
                zrange = colorlimits;
            else
                zrange = [min(tmp) max(tmp)];
            end
            zh = zrange(1):(zrange(2)-zrange(1))./249:zrange(2);
            zh = round(zh*100);
            if isempty(zh)  % we probably have all same values
                zh = ones(1, 250);
            end
        else
            zh = [];
        end

        if ~isempty(tmpc)
            if colorlimits
                zrangec = -colorlimits;
            else
                zrangec = [min(tmpc) max(tmpc)];
            end
            zhc = zrangec(1):(zrangec(2)-zrangec(1))./249:zrangec(2);
            if isempty(zhc), zhc = 1:length(hc);end
            zhc = round(zhc*100);
            if isempty(zhc)  % we probably have all same values
                zhc = ones(1, 250);
            end
        else
            zhc = [];
        end

        if isempty(zh) && isempty(zhc), error('No Z-scores in cluster'), end

    else

        % surface patch method
        % ----------------------------------------------------------------------------------------
        vol = voxel2mask(XYZ', odims);
        vol = smooth3(vol);

    end

    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    % Loop through slices
    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    if length(whsl) > 12, fsz = 14; else fsz = 18;end
    index = 1;
    for z = whsl
        set(0, 'CurrentFigure', mch);
        subplot(rc, rc, index);

% This is now done earlier, in display_underlay
%         if ~isempty(oimg)
%             imagesc(oimg(:, :, z)')
%             set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
%             hold on; axis image; axis off
%             zmm = voxel2mm([1 1 z]', V.mat);
%             text(textx, texty, ['z = ' num2str(zmm(3))], 'Color', 'k', 'FontSize', fsz)
%         else
%             hold on
%         end

        doimgpatch = 0;

        if doimgpatch   %isempty(colorlimits)

            % --------------------------------------------------------------------------------
            % solid color patches
            % --------------------------------------------------------------------------------
            if z>1
                mvol = vol(:, :, z-1:z); for i = 1:size(mvol, 3), myvol(:, :, i) = mvol(:, :, i)';end
                %FVC = isocaps(vol(:, :, z-1:z)', 0, 'zmax');
            else
                mvol = vol(:, :, z:z+1); for i = 1:size(mvol, 3), myvol(:, :, i) = mvol(:, :, i)';end
                %FVC = isocaps(vol(:, :, z:z+1)', 0, 'zmax');
            end
            FVC = isocaps(myvol, 0, 'zmax');

            try
                patch(FVC, 'EdgeColor', 'none', 'FaceColor', myc, 'FaceAlpha', 1)
            catch
                patch(FVC, 'EdgeColor', 'none', 'FaceColor', myc)
            end

        else
            % --------------------------------------------------------------------------------
            % color-mapped points
            % --------------------------------------------------------------------------------

            myxyz = XYZ(1:2, XYZ(3, :) == z);

            if ~isempty(colorlimits)
                myz = allz(XYZ(3, :) == z);
                clear h2, clear wh
            end

            %plot(myXYZ(2, :), myXYZ(1, :), [myc 's'], 'MarkerFaceColor', myc, 'MarkerSize', 3)
            hold on
            for j = 1:size(myxyz, 2)

                % added to replace image patch, solid color point fill for pos Z values only
                if isempty(colorlimits)
                    Xfill = [myxyz(1, j) myxyz(1, j) myxyz(1, j)+xs(1) myxyz(1, j)+xs(1)];
                    Yfill = [myxyz(2, j) myxyz(2, j)+xs(2) myxyz(2, j)+xs(2) myxyz(2, j)];
                    h2(j) = patch(Xfill, Yfill, myc, 'EdgeColor', 'none');

                else

                    if myz(j) < 0
                        tmp = find((zhc-round(myz(j)*100)).^2 == min((zhc-round(myz(j)*100)).^2));
                        wh(j) = tmp(1);
                        %h2(j) = plot(myxyz(1, j), myxyz(2, j), 'Color', hc(wh(j), :), 'MarkerSize', 3, 'MarkerFaceColor', hc(wh(j), :));


                        h2(j) = fill([myxyz(1, j) myxyz(1, j) myxyz(1, j)+xs(1) myxyz(1, j)+xs(1)], [myxyz(2, j) myxyz(2, j)+xs(2) myxyz(2, j)+xs(2) myxyz(2, j)], hc(wh(j), :), 'EdgeColor', 'none');
                    else
                        tmp = find((zh-round(myz(j)*100)).^2 == min((zh-round(myz(j)*100)).^2));
                        wh(j) = tmp(1);
                        %h2(j) = plot(myxyz(1, j), myxyz(2, j), 'Color', h(wh(j), :), 'MarkerSize', 3, 'MarkerFaceColor', h(wh(j), :));
                        h2(j) = patch([myxyz(1, j) myxyz(1, j) myxyz(1, j)+xs(1) myxyz(1, j)+xs(1)], [myxyz(2, j) myxyz(2, j)+xs(2) myxyz(2, j)+xs(2) myxyz(2, j)], h(wh(j), :), 'EdgeColor', 'none');
                    end
                    %if exist('h2') == 1, set(h2, 'Marker', 'square'), end

                end %if isempty colorlimits
            end

        end % if solid or color-mapped


        % text numbers, if plotting text at xyz coords of centers (varargin)
        %if length(varargin) > 0
        %   cencooz = cencoo(:, cencoo(3, :) == z);
        %  for i = 1:size(cencooz, 2)
        %     text(cencooz(2, i), cencooz(1, i), num2str(cenind), 'Color', 'k')
        %    cenind = cenind + 1;
        %end
        %end

        index = index + 1;
        %drawnow

    end % slice loop


    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    % Text and color bars
    % * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    if isempty(varargin)
        set(0, 'CurrentFigure', mch);
        subplot(rc, rc, index)
        a = pwd; a = a(end-6:end);
        
        if isfield(clusters(1), 'threshold') && isnumeric(clusters(1).threshold)
            b = num2str(clusters(1).threshold);
        else
            b = 'none';
        end

        axis off
        c = num2str(length(clusters));
        if ~isfield(clusters, 'title') || isempty(clusters(1).title)
            clusters(1).title = 'Clusters';
        end

        text(0, clind-1, [myc ': ' clusters(1).title ' ' a ' u = ' b ', ' c ' clusters'])
        axis([0 1 -1 clind])

        if ~isempty(colorlimits)
            if ~isempty(tmp)
                % make color bar
                %bfig = figure('Color', 'w');
                bax = axes('Position', [0.100    0.0800    0.40    0.02]);
                hold on;

                if any(allz > 0)
                    zh2 = zh./100;
                    for i = 2:size(h, 1), fill([zh2(i-1) zh2(i-1) zh2(i) zh2(i)], [0 1 1 0], h(i, :), 'EdgeColor', 'none'); end
                    plot([min(zh2) max(zh2)], [0 0], 'k-')
                    set(gca, 'YTickLabel', '');
                    %xlabel('Z-score', 'FontSize', 14)
                end
            end

            if any(allz < 0)
                %if ~exist('bfig'),
                %bfig = figure('Color', 'w'); hold on;
                %end
                zh2 = zhc./100;
                for i = 2:size(hc, 1), fill([zh2(i-1) zh2(i-1) zh2(i) zh2(i)], [0 1 1 0], hc(i, :), 'EdgeColor', 'none'); end
                plot([min(zh2) max(zh2)], [0 0], 'k-')
                set(gca, 'YTickLabel', '');
                %xlabel('Z-score', 'FontSize', 14)
            end

        end

        %set(gca, 'Position', [0.1300    0.1100    0.7750    0.8150])
        set(gca, 'FontSize', 16, 'Color', [1 1 1])
        %xlabel('Z-score', 'FontSize', 16)
        %set(gca, 'Position', [0.1000    0.4500    0.80    0.50])
        %ssize = get(0, 'ScreenSize');
        %set(mch, 'Position', ssize./[0.0370    0.0014    2.7948   15.2500])
    end

    drawnow

end



function display_underlay(mch, oimg, whsl, rc, V, textx, texty)
% mch: figure handle
% oimg: underlay image matrix 3-D
% whsl: which slices
% rc: rows/cols for subplots (integer number)
% V: volume structure with .mat field for underlay image
% textx, texty: x- and y-axis positions for text labels on plot

% tor: sept 2010: try to make consistent with spm flipping based on x
% coordinate
if sign(V.mat(1)) < 0
    for z = 1:size(oimg, 3)
        oimg(:, :, z) = flipud(oimg(:, :, z));
    end
end

if length(whsl) > 12, fsz = 14; else fsz = 18; end

%     whsl(whsl <=0) = [];
%     whsl(whsl >= size(oimg, 3)) = [];

    index = 1;
    for z = whsl
        set(0, 'CurrentFigure', mch);
        subplot(rc, rc, index);

        if ~isempty(oimg)
            imagesc(oimg(:, :, z)')
            set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
            hold on; axis image; axis off
            zmm = voxel2mm([1 1 z]', V.mat);
            text(textx, texty, ['z = ' sprintf('%3.0f', zmm(3))], 'Color', 'k', 'FontSize', fsz)
        else
            hold on
        end
        
        index = index + 1;
    end
    
end



function ph = plot_points(mch, XYZmm, rc, V, whsl, myc, myplab)

    XYZ = mm2voxel(XYZmm, V, 1)';     % suppress unique voxel output
    index = 1;
    phind = 1;

    for z = whsl
        set(0, 'CurrentFigure', mch);
        subplot(rc, rc, index);
        hold on
        myXYZ = XYZ(:, XYZ(3, :) == z);
        if ~isempty(myplab), myplz = myplab(XYZ(3, :) == z); end

        for i = 1:size(myXYZ, 2)
            ph(phind) = plot3(myXYZ(1, i), myXYZ(2, i), 100, [myc(1) '.'], 'MarkerFaceColor', myc(1), 'MarkerSize', 8);

            if ~isempty(myplab)
                text(myXYZ(1, i), myXYZ(2, i), 100, myplz{i}, 'Color', myc(1))
            end

            phind = phind + 1;
        end

        index = index + 1;
        %view(0, 90)
        %plot(0, 0, 'kd')
    end
end






function in = scn_mat_conform(in)
    %function in = scn_mat_conform(in)
    %
    % sets flipping to 0 (no flip) in SPM2 and adjusts mat file accordingly
    % input in spm-style mat file or struct with .mat or .M fields
    %

    global defaults
    spm_defaults();

    if isstruct(in)
        if isfield(in,'mat')
            in.mat = scn_mat_conform(in.mat);
        end

        if isfield(in,'M')
            in.M = scn_mat_conform(in.M);
        end

        return
    end



    if in(1) < 0
        disp('Orientation: Image has negative x voxel size, indicating radiological orientation.'); %''flipped'' in SPM.')
         disp('The x voxel size will be changed to positive for display. Montage_clusters program does not flip images for display.')
 
         in(1) = abs(in(1));     % voxel size
         in(1,4) = -(in(1,4));   % origin offset
    else
                disp('Orientation: Image has positive x voxel size, indicating neurological orientation, and will not be flipped for display.'); %''flipped'' in SPM.')
    end

    if isempty(defaults) || ~isfield(defaults, 'analyze')
        spm_defaults
    end

    switch spm('Ver')
        case 'SPM8'
            % no analyze field; flip is always 1

        case {'SPM5', 'SPM2', 'SPM99'}
            if defaults.analyze.flip
                disp('Warning: Setting defaults.analyze.flip to 0.  No flipping.')
                defaults.analyze.flip = 0;
            end
        otherwise
            warning('Unknown version of SPM. Results may be erratic.')
    end

end
    
