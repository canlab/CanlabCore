function [data_matrix12, xvals2, yvals2, data_matrix22, bestx, besty, bestz] = surf_plot_tor(data_matrix1, xvals, yvals, xname, yname, zname, varargin)
% Stylized surface plot of one or two surfaces
%
% :Usage:
% ::
%
%    [data_matrix12, xvals2, yvals2, data_matrix22] = surf_plot_tor(data_matrix1, xvals, yvals, xname, yname, zname, [data_matrix2])
%
% :Examples: See classify_search_script3...in meta-analysis classification
% ::
%
%    surf_plot_tor(corrc_mean, mya, mys, 'Activation feature cutoff', 'Sensitivity feature cutoff', 'Classification accuracy', worstcat)
%
% xvals are columns, yvals are rows!
%
% ..
%    Tor Wager, Sept. 07
% ..

    data_matrix22 = [];

    ystep = 1/5 .* range(yvals) ./ length(yvals); 
    xstep = 1/5 .* range(xvals) ./ length(xvals); 
    
    [X,Y] = meshgrid(xvals,yvals);
    
%     if size(X, 1) == size(data_matrix1, 2)
%         error('Wrong sizes...are the x and y inputs flipped?');
%     end
    
    yvals2 = 0:ystep:max(yvals);
    xvals2 = 0:xstep:max(xvals);

    [X2,Y2] = meshgrid(xvals2,yvals2);
    data_matrix12 = interp2(X,Y,data_matrix1,X2,Y2,'cubic');

    create_figure('Surface plot');

    han = surf(X2,Y2,data_matrix12);

    xlabel('Activation feature cutoff')
    ylabel('Sensitivity feature cutoff')
    zlabel('Classification accuracy');

    xlabel(xname);
    ylabel(yname);
    zlabel(zname);

    set(han,'EdgeColor','none')
    grid on

    % lines
    plot3(X2(end,:),Y2(end,:),data_matrix12(end,:),'r','LineWidth',2);
    plot3(X2(end,:),Y2(end,:),data_matrix12(end,:),'r','LineWidth',3);

    for i = 10:10:size(X2,1)
        plot3(X2(i,:),Y2(i,:),data_matrix12(i,:),'Color',[.5 .5 .5],'LineWidth',1);
    end

    view(135, 30)
    drawnow

    scn_export_papersetup(600);

    if length(varargin) > 0
        data_matrix2 = varargin{1};

        data_matrix22 = interp2(X,Y,data_matrix2,X2,Y2,'cubic');
        hold on
        han2 = surf(X2,Y2,data_matrix22);
        set(han2,'EdgeColor','none')
        set(han,'FaceAlpha',.7)

        % lines
        plot3(X2(end,:),Y2(end,:),data_matrix22(end,:),'b','LineWidth',3);
        plot3(X2(1,:),Y2(1,:),data_matrix22(1,:),'k','LineWidth',1);


        for i = 10:10:size(X2,1)
            plot3(X2(i,:),Y2(i,:),data_matrix22(i,:),'Color',[.5 .5 .5],'LineWidth',1);
        end


    end
    
    % calc max and put max point on map
    % --------------------------------------------------
    if length(varargin) > 0
        % average two maps
        mymap = .5 * (data_matrix12 + data_matrix22);
        fprintf('Using average of two maps to calculate max.');
    else
        mymap = data_matrix12;
    end
    
    [bestsum] = max(mymap(:));  
    [row, col] =  find(mymap == bestsum); 
bestx = xvals2(col);
besty = yvals2(row);
bestz = mymap(row, col);

fprintf('Maximum: x = %3.3f, y = %3.3f, z = %3.3f\n', bestx, besty, bestz);

plotz = bestz;
if length(varargin) > 0
    max1 = data_matrix12(row, col);
    max2 = data_matrix22(row, col);
    plotz = max([max1 max2]);
    
    bestz = [bestz max1 max2];
    fprintf('Height at overall max: Map 1: %3.3f, Map 2: %3.3f\n', max1, max2);
end

% plot
plot3(bestx, besty, plotz, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k')

mylowerbound = get(gca,'ZLim'); mylowerbound = mylowerbound(1);
plot3([bestx bestx], [besty besty], [plotz mylowerbound], 'k-', 'LineWidth', 2);

myxlowerbound = get(gca,'XLim'); %myxlowerbound = myxlowerbound(1);
plot3([myxlowerbound], [besty besty], [mylowerbound mylowerbound], 'k-', 'LineWidth', 2);   

myylowerbound = get(gca,'YLim'); %myylowerbound = myylowerbound(1);
plot3([bestx bestx], [myylowerbound], [mylowerbound mylowerbound], 'k-', 'LineWidth', 2); 

set(gca, 'ZLim', get(gca, 'ZLim'));
axis vis3d

 
    

end

