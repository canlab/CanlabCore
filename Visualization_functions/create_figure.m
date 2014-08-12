function f1 = create_figure(tagname, varargin)
    % f1 = create_figure(['tagname'], [subplotrows], [subplotcols], [do not clear flag])
    %
    % checks for old figure with tag of tagname
    % clears it if it exists, or creates new one if it doesn't
    
    if nargin < 1 || isempty(tagname)
        tagname = 'nmdsfig';
    end
    
    doclear = 1;    % clear if new or if old and existing
    if length(varargin) > 2 && varargin{3}
        % use same figure; do not clear
        doclear = 0;
    end
        
    old = findobj('Tag', tagname);
    old = old( strcmp( get(old, 'Type'), 'figure' ) );
    
    if ~isempty(old)

        if length(old) > 1
            % multiple figures with same tag!
            close(old(2:end))
            old = old(1);
        end
            
        if doclear, clf(old); end

        f1 = old;
        
    else
        % Or create new
        
        scnsize = get(0,'ScreenSize');

        xdim = min(scnsize(3)./2, 700);
        ydim = min(scnsize(4)./2, 700);

        f1 = figure('position',round([50 50 xdim ydim]),'color','white');
        set(f1, 'Tag', tagname, 'Name', tagname);
    end
    
    % activate this figure
    figure(f1);

        
    if doclear % true for new figs or cleared ones

        % Create subplots, if requested; set axis font sizes

        if length(varargin) > 0
            i = max(1, varargin{1});
            j = max(1, varargin{2});
        else
            i = 1;
            j = 1;
        end

        np = max(1, i * j);

        for k = 1:np
            axh(k) = subplot(i,j,k);
            cla;
            set(gca,'FontSize',18),hold on
        end
        axes(axh(1));

    end

    return