function cl = classify_viz_regions(indx,colors,ptaskcutoff,sizecutoff,maskimage,names)
% Visualize output of classify_bayes.m, or comparable output
%
% :Usage:
% ::
%
%     cl = classify_viz_regions(indx,colors,ptaskcutoff,sizecutoff,maskimage,names)
%
% indx should be voxels x images, with prob task (ptask) in non-zero
% elements
%
% :Examples:
% ::
%
%    cl = classify_viz_regions(indx,[],.6,5);
%    cl = classify_viz_regions(indx,[],.6,5,[],names);
%
% ..
%    tor wager
% ..

cl = [];

    if nargin < 2 || isempty(colors), colors = {[1 0 0] [0 1 0] [1 0 1] [1 1 0] [0 0 1]}; end
    if nargin < 3 || isempty(ptaskcutoff),ptaskcutoff = .4; end
    if nargin < 4 || isempty(sizecutoff),sizecutoff = 3; end
    if nargin < 5 || isempty(maskimage),maskimage = which('scalped_avg152T1_graymatter_smoothed.img'); end
    if nargin < 6, names = []; end
    
    [nvox,ntasks] = size(indx);

    while length(colors) < ntasks
        colors{end+1} = rand(1,3);
    end
    
    volInfo = iimg_read_img(maskimage,2);

    if (nvox ~= volInfo.nvox) &&  (nvox ~= volInfo.n_inmask)
        disp('Cannot plot because # vars does not match # voxels in mask.'); 
        return
    end

    addstr = 'noadd';
    
    % light colors, first threshold
    indxlite = indx;
    indxlite(indxlite >= ptaskcutoff) = 0;
    
    
    for i = 1:ntasks

        cl{i} = iimg_indx2clusters(indxlite(:,i),volInfo);

        if ~isempty(cl{i})
            sz = cat(1,cl{i}.numVox);
            sz = sz < sizecutoff;
            cl{i}(sz) = [];

            mycolor = {(.7 .* colors{i})}; %^+ ([1 1 1].*.5)};
            cluster_orthviews(cl{i},mycolor,'solid',addstr);
            addstr = 'add';
        end
    end
    
    % solid colors, higher threshold
    indx(indx < ptaskcutoff) = 0;
    
    for i = 1:ntasks

        cl{i} = iimg_indx2clusters(indx(:,i),volInfo,ptaskcutoff,sizecutoff);

        if ~isempty(cl{i})
            sz = cat(1,cl{i}.numVox);
            sz = sz < sizecutoff;
            cl{i}(sz) = [];

            cluster_orthviews(cl{i},colors(i),'solid',addstr);
            addstr = 'add';
        end
        
    end

    if ~isempty(names), makelegend(names,colors,1); end

    return
