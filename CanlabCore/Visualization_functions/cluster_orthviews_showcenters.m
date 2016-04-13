function slices_fig_h = cluster_orthviews_showcenters(cl, myview, overlay, xhairs, sliceorder, bgcolor, varargin)
% :Usage:
% ::
%
%    slices_fig_h = cluster_orthviews_showcenters(cl, myview, [overlay], [xhairs], [order slices flag], [background color], [var args])
%
%    cluster_orthviews_showcenters(cl, 'coronal');
%    cluster_orthviews_showcenters(cl, 'sagittal');
%    cluster_orthviews_showcenters(cl, 'axial');
%
% used in cluster_orthviews_classes
%
% ..
%    tor wager, aug 2006
% ..

    if nargin < 3 || isempty(overlay), overlay = which('scalped_single_subj_T1.img'); end
    if nargin < 4, xhairs = 1; end
    if nargin < 5, sliceorder = 0; end
    
    whichorth = [];
    for i = 1:length(varargin)
       if ischar(varargin{i})
            switch varargin{i}

                case 'whichorth', whichorth = varargin{i+1};
                
                otherwise, warning('scn_tools:badInput', ['Unknown input string option:' varargin{i}]);
            end
        end
    end 
    
    % set up orthviews to be OK
    scn_export_spm_window('setup', overlay);
%     set(get(findobj('Tag','Graphics'), 'Children'), 'Color', 'white');
    set(findobj('Type', 'axes', 'Parent', findobj('Tag','Graphics')), 'Color', 'white');
    if xhairs, spm_orthviews('Xhairs', 'on'); end

    % Set window to be equal/constant
    [volInf_tmp, dat] = iimg_read_img(overlay, 2);
    dat = dat(volInf_tmp.wh_inmask);
    spm_orthviews('Window', whichorth, [prctile(dat, 2) prctile(dat, 98)]);

    myviews = {'axial' 'coronal' 'sagittal'};   % for selecting SPM window
    whview = find(strcmp(myviews, myview));
    if isempty(whview), error('myview must be axial, coronal, or sagittal.'); end

    if sliceorder
        % order slices from negative to positive along this dimension
        myviews2 = {'sagittal' 'coronal' 'axial' };  % for selecting coord
        whcoord = strmatch(myview, myviews2) ;
        cen = cat(1, cl.mm_center);
        [cen, wh] = unique(cen(:, whcoord));
        % remove duplicates and sort
        cl = cl(wh);

        % get text string base
        mystr = {'x = ' 'y = ' 'z = '};
        textbase = mystr{whcoord};
    end


    % get optimal number of axes
    num_clusters = length(cl);
    rc = ceil(sqrt(num_clusters));

    if isempty(whichorth)
        axh = get_orth_axishandles;
    else
        axh = get_orth_axishandles_whichorth(whichorth);
    end

    axh = axh(whview);

    slices_fig_h = figure; %create_figure(myview);
    set(slices_fig_h, 'Color', 'k');
    
    for i = 1:num_clusters
        newax(i) = subplot(rc, rc, i);
        axis off;
     end


    for i = 1:num_clusters
        %cluster_orthviews(cl{i}, colors2(classes(i)), 'overlay', overlay);
        spm_orthviews('Reposition', cl(i).mm_center);

        copyobj(get(axh, 'Children'), newax(i));
        axes(newax(i));
        axis image

        h = findobj(newax(i), 'Type', 'text');
        h= h(1); % kludgy fix for multiple matches ***
        
        if sliceorder
            % try to set a reasonable font size
            pos = get(newax(i), 'Position');
            height = pos(3);
            fs = round(height * 70);
            set(h, 'FontSize', fs)
            % set position of text (move down and right)
            pos = get(h, 'Position');
            pos(2) = 0; %pos(2) - fs./2;
            pos(1) = 0;
            set(h, 'Position', pos)
            % set text string based on cen
            set(h, 'String', [textbase num2str(cen(i))]);
        elseif ishandle(h)
            delete(h);
        end
    end

    % try to set a reasonable enlargement factor
    n = 1 + .15 * log(num_clusters ./ 2);
    n = max(n, 1); n = min(n, 3);
    enlarge_axes(gcf, n)

    % set background color to print as whatever it is on screen
    set(gcf,'InvertHardcopy', 'off');
end




function axish = get_orth_axishandles

    % Get figure and axis handles
    fh = findobj('Tag', 'Graphics');
    ch = get(fh, 'Children');
    for i= 1:length(ch)
        mytype = get(ch(i), 'Type');
        wh(i) = strcmp(mytype, 'axes');
    end
    axish = ch(find(wh));

    if isempty(axish), error('SPM figure orthviews do not exist'); end

    % get which axis is which
    for i = 1:length(axish)
        poss(i, :) = get(axish(i), 'Position');
    end

    % get rid of extra axes we may have created in the 4th quadrant
    other_axes = find(any(poss(:, 1) > .45 & poss(:, 2) < .2, 2));
    axish(other_axes) = [];
    poss(other_axes, :) = [];

    % sort into order:  axial, coronal, saggital
    ssum = sum(poss(:, 1:2), 2);
    [ssum, ind] = sort(ssum);
    axish = axish(ind);

end

function axish = get_orth_axishandles_whichorth(whichorth)
% Get the axis handles for the current orthviews
global st
for i = 1:length(st.vols), wh(i) = ~isempty(st.vols{i}); end
wh = find(wh); wh = wh(whichorth);
axish = cat(1, st.vols{wh}.ax{:});
axish = sort(cat(1, axish(:).ax));

end

