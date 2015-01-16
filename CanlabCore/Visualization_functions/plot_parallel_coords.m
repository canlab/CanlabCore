function plot_parallel_coords(dat)

    select_on = 'crossproducts';  % 'crossproducts' or 'variables'
    mcolor = [.2 .2 .2];
    facecolor = [.5 .5 .5];
    mtype = 'o';
    cutoff_percentile = 0;
    select_color_on_column = 1;
    
    % init figure and sizes
    % ---------------------------------------------
    create_figure('Parallel coords plot');

    [nobs, nvars] = size(dat);
    set(gca,'XLim',[0.5 nvars + .5], 'XTick', 1:nvars);


    % cross-products, if needed
    % ---------------------------------------------
    if strcmp(select_on, 'crossproducts')
        for i = 1:nvars - 1
            xp(:,i) = scale(dat(:,i)) .* scale(dat(:,i + 1));
        end
    end


    % initial plot
    % ---------------------------------------------

    run_plot;

    % color initial high values or xproducts
    % ---------------------------------------------
    
    cutoff_percentile = 66;
    facecolor = 'r';

    run_plot;

        % color last var high values or xproducts
    % ---------------------------------------------
    switch select_on
        case 'variables'
            select_color_on_column = nvars;
            xlabel('Red: high on v(1)  Blue: high on v(end)')
        case 'crossproducts'
            select_color_on_column = nvars - 1;
            xlabel('Red: high on xp(1,2)  Blue: high on xp(end-1:end)')
        otherwise 
            error('Unknown select_on!')
    end
    facecolor = 'b';

    run_plot;



    
    %
    %
    %
    % inline
    %
    %
    %
    
    function run_plot

        wh = find(dat(:, select_color_on_column) > prctile(dat(:, select_color_on_column), cutoff_percentile));

        if strcmp(select_on, 'crossproducts')
            wh = find(xp(:, select_color_on_column) > prctile(xp(:, select_color_on_column), cutoff_percentile));
        end

        for i = 1:nvars
            plot(i, dat(wh, i), mtype, 'Color', mcolor, 'MarkerFaceColor', facecolor)
        end
        
        for i = 1:length(wh)
            plot(dat(wh(i), :), '-', 'Color', facecolor);
        end
        
        drawnow

    end


end
