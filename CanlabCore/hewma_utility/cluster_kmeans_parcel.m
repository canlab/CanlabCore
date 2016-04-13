function [cl2,nclasses,colors] = cluster_kmeans_parcel(x,CLU,doplot,varargin)
% K-means clustering of voxels
%
% :Usage:
% ::
%
%     [cl2,nclasses,colors] = function cluster_kmeans_parcel(x,CLU,doplot,[overlay img])
%
% :Inputs:
%
%   **x:**
%        a voxels x data matrix
%
%   **CLU:**
%        a structure containing a list of XYZ coordinates and z-scores for voxels
%        (see clusters2CLU.m)
%
%   a plot flag
%
% :Outputs:
%
%   **cl2:**
%        a k-length cell array of clusters structures
%        each cell contains a clusters structure for one data class
%
% ..
%    tor wager
% ..

if length(varargin), ovl = varargin{1}; else, ovl = which('scalped_single_subj_T1.img');, end

doscale = 1;    % standardize each column of x before clustering 


if doplot
    
    tor_fig; view(135,30); axis vis3d

    % ----------------------------------------------------
    % make initial histogram 
    % ----------------------------------------------------


    nbins = [length(unique(x(:,1))) length(unique(x(:,2)))];
    [h,c] = hist3(x(:,1:2),nbins);
    h = log(h);

    %han = axes('Position',[.52 .08 .45 .38]);
    set(gca,'FontSize',16)
    barh = bar3(h);
    title('Log frequency')
    %set(gca,'ZLim',[0 max(h(:))]);
    set(barh,'FaceColor',[.5 .5 .5],'EdgeColor','none'); lighting phong; camlight right; [az,el] = view; lightangle(az,el);
    xlabel('Activation duration'),set(gca,'XTick',1:20:nbins(2),'XTickLabel',round(c{2}(1:20:nbins(2))));
    ylabel('Change Point'),set(gca,'YTick',1:20:nbins(1),'YTickLabel',round(c{1}(1:20:nbins(1))));
    
    
    % fix axis limits
    xl = get(gca,'XLim'); yl = get(gca,'YLim');
    set(gca,'ZLim',[0 max(h(:))+1]);
    set(gca,'XLim',xl)
    set(gca,'YLim',yl)
    set(gca,'DataAspectRatio',[1 1 .1])
end


% ----------------------------------------------------
% k-means clustering
% ----------------------------------------------------

nclasses = input('How many classes?');

err = 1; indx = 1;
while err
    try
        if doscale
            classes = kmeans(zscore(x), nclasses);   % ,'start','uniform');
        else
            classes = kmeans(x, nclasses); 
        end
        err = 0;
    catch
    end
    indx = indx + 1;
    if indx == 11, disp('kmeans: tried 10 times.  No solution.'); err = 0;, return, end
end

if doplot

    % ----------------------------------------------------
    % define colors and sort by class size
    % ----------------------------------------------------

    colors = {[1 0 0] [0 1 0] [0 0 1] [1 1 0] [1 .5 0] [0 1 1] [.5 1 0] [1 0 1]};
    while length(colors) < nclasses, colors = [colors {rand(1,3)}];,end

%     for i = 1:nclasses, nvox(i) = sum(classes==i);,end
%     [nvox,i] = sort(nvox,2,'descend');
% 
%     colors(i) = colors(1:length(i));

    nvars = size(x,2);

    for i = 1:nclasses,
        for j = 1:nvars
            meanvar{j}(i) = mean(x(classes==i,j));,
            stdvar{j}(i) = std(x(classes==i,j));,
        end
        nvox(i) = sum(classes==i);
        indx(i) = i;
    end
    % sort by frequency, descending
    [nvox,i] = sort(nvox,2,'descend');
    for j = 1:nvars
        meanvar{j} = meanvar{j}(i);
        stdvar{j} = stdvar{j}(i);
    end
    indx = indx(i);

    % ----------------------------------------------------
    % print table of classes
    % ----------------------------------------------------
    
    fprintf(1,'\nClass\tNum. Voxels\t');
    for j = 1:nvars, fprintf(1,'V%02d\t',j), end
    fprintf(1,'Color\t\t\n');
    for i =1:nclasses
        fprintf(1,'%3.0f\t%3.0f\t',indx(i),nvox(i));
        for j = 1:nvars, fprintf(1,'%3.2f (%3.2f)\t',meanvar{j}(i),stdvar{j}(i)), end      
        fprintf(1,'%3.2f\t%3.2f\t%3.2f\t\n',colors{i}(1),colors{i}(2),colors{i}(3));
    end
    fprintf(1,'\n');
    
    % colors are in sorted order. Re-order so that they correspond to original cl order 
    colors(indx) = colors(1:length(indx));

    % ----------------------------------------------------
    % bar color change
    % ----------------------------------------------------

    for cc = 1:nclasses
        hist_indic{cc} = zeros(size(h));
    end

    for ii = 1:size(h,1)
        for jj = 1:size(h,2)

            if h(ii,jj) > 0

                bincenter = [c{1}(ii) c{2}(jj)];

                % class for this bin
                d = distance(bincenter,x);
                wh = find(d == min(d));
                clas = classes(wh(1));

                % put in histogram indicator
                hist_indic{clas}(ii,jj) = 1;

            end
        end
    end


    % plot bars in color
    for cc = 1:nclasses

        htmp = h .* hist_indic{cc};

        hold on;
        hhtmp = bar3(htmp);
        set(hhtmp,'FaceColor',colors{cc},'EdgeColor','none')

    end

    hhtmp = bar3(zeros(size(h)));
    set(hhtmp,'FaceColor',[.8 .8 .8])

    lighting phong; camlight right; [az,el] = view; lightangle(az,el);
    xlabel('Activation duration'),set(gca,'XTick',1:20:nbins(2),'XTickLabel',round(c{2}(1:20:nbins(2))));
    ylabel('Change Point'),set(gca,'YTick',1:20:nbins(1),'YTickLabel',round(c{1}(1:20:nbins(1))));
    

    drawnow


    h1 = get(gca,'Children');

end     % if doplot



% ----------------------------------------------------
% 1 - D colored histograms
% ----------------------------------------------------

% ----------------------------------------------------
% re-plot histogram with color codes
% ----------------------------------------------------
% nplot = size(x,2);  tor_fig(nplot,1);
% for i = 1:nclasses
%     wh = find(classmap == i); range = [min(wh) max(wh)]; 
%     wh = find(x <= range(2) & x >= range(1));
%     hh = bar(x(wh),h(wh)); set(hh,'FaceColor',colors{i});
% end
% xlabel('Change point')
% ylabel('Number of voxels')







CLU.Z = classes';

% ----------------------------------------------------
% re-make separate clusters for each class
% and plot on brain
% ----------------------------------------------------

clear cl2 
for i = 1:nclasses
    CLUtmp = CLU;  
    wh = find(CLUtmp.Z == i);
    CLUtmp.XYZmm = CLUtmp.XYZmm(:,wh);
    CLUtmp.XYZ = CLUtmp.XYZ(:,wh);
    CLUtmp.Z = CLUtmp.Z(:,wh);
    %CLUtmp.Z = [CLUtmp.Z; x(wh,:)'];
    cl2{i} = tor_extract_rois([],CLUtmp,CLUtmp);
    
    if doplot
        if i == 1
            cluster_orthviews(cl2{i},colors(i),'overlay',ovl);
        else
            cluster_orthviews(cl2{i},colors(i),'add');
        end
    end
end

if doplot
    f1 = findobj('Tag','Graphics'); % spm fig window handle

    figure(f1)
    han = axes('Position',[.52 .08 .45 .38]);
    try
        copyobj(h1,han); view(135,30); axis vis3d

        set(gca,'FontSize',16)
        title('Log frequency')
        xlabel('Activation duration')
        ylabel('Change Point')

        xlabel('Activation duration'),set(gca,'XTick',1:20:nbins(2),'XTickLabel',round(c{2}(1:20:nbins(2))));
        ylabel('Change Point'),set(gca,'YTick',1:20:nbins(1),'YTickLabel',round(c{1}(1:20:nbins(1))));
        drawnow
    catch
        disp('Error copying object to SPM window.')
    end
end


return

