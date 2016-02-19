function [ds, g, mystd, d, d2, c, c2, mi, b, eigv, eigval] = compare_subjects(varargin)
% This function compares a set of images to one another and does some diagnostics on the similarity among images.
% - It returns multivariate distances and dissimilarities among images
% - It works on the GLOBAL signal after standardizing each image (case 1) or the REGIONAL values in each cluster (case 2) 
% - You can also enter a reference image, in which case each image will be correlated with the ref.
%
%
% :Usage:
% ::
%
%     function [ds, g, mystd, d, d2, c, c2, mi, b, eigv, eigval] = compare_subjects([img files or clusters], [mask], ...
%                                        [plot flag], [title on figure], [standardize flag], [text labels], [ref image])
%
% :Inputs:
%
%     a list of image names to compare
%
%     OR
%
%     a clusters structure, with data to compare
%     in timeseries field
%
% If a mask is entered, only voxels in the mask (e.g., with value of 1) will be used.
% You can use this option to specify brain-only or gray-matter only voxels
%
% textlab: optional text labels for each image, can be empty []
%
% If a ref image is entered, each image will be correlated with the ref,
% and values will be saved for the correlation (plot 2 will show these values)
% Useful for comparing anatomical imgs with template, etc.
%
% :Outputs:  from correls with ref image are in variable "c"
%
%   **ds:**
%        multivariate distance (sim. to Mahalanobis) for each image
%        ds is a matrix of squared distances, case numbers, and
%        expected chi2 values (in columns in this order) rows are cases
%
%   **g:**
%        global value for each image
%
%   **d:**
%        global distance from mean image
%        distance, or dissimilarity, is the average absolute deviation between images
%
%   **d2:**
%        matrix of distances among all images
%
%   **c:**
%        correlation between real valued voxels and mean image
%
%   **c2:**
%        correlations among all images (treating voxels as cases)
%
%   **mi:**
%        mutual information between images, with hist2.m
%
%   **b:**
%        principal component scores on correlation matrix for eigenvalues > 1
%
%   **eigv:**
%        eigenvectors
%
%   **eigval:**
%        eigenvalues
%
% :Examples:
% ::
%
%    % Compare normalized anatomcals with standard brain
%    P = get_filename2(['sub*\Anatomy\nscalped_ft1.img']);
%    [ds, g, mystd, d, d2, c, c2, mi] = compare_subjects(P, which('brain_avg152T1.img'), 1, 'intext_countloc', 1, [], which('avg152T1.img'));
%
% ..
%    by Tor Wager
% ..


    doplot = 1; 
    dostd = 0;
    do_mi = 1;
    do_corr = 1;
    do_chi2 = 1;
    textlab = [];
    if isempty(varargin)
        error('no inputs.');
    end

    if length(varargin) > 2, doplot = varargin{3}; end
    if length(varargin) > 3, mytitle = varargin{4}; end
    if length(varargin) > 4, dostd = varargin{5}; end
    if length(varargin) > 5, textlab = varargin{6}; end
    if length(varargin) > 6, refimg = varargin{7}; end
    
    for i=8:length(varargin)
        if(ischar(varargin{i}))
            switch(varargin{i})
                case 'no_mi'
                    do_mi = 0;
                case 'no_corr'
                    do_corr = 0;
                case 'no_chi2'
                    do_chi2 = 0;
            end
        end
    end




    if isstruct(varargin{1})
        clusters = varargin{1};
        error('No method implemented yet.');
    else
        hP = varargin{1};

        disp(' compare subjects.m Running ');
        
        % Check for images whose dims do not match
        % -----------------------------------------------------------------
        V = spm_vol(hP);
        [bad, whbad] = iimg_check_volinfo(V(1), V);

        hP(whbad, :) = [];
        V(whbad) = [];
        
        if any(bad)
            disp('Warning!!! Some image dims do not match.  Removing subjects:');
            fprintf('%3.0f ', whbad);
            fprintf('\n');

            if ~isempty(textlab),
                fprintf('%s ', textlab{whbad});
                fprintf('\n');
                textlab(whbad) = [];
            end

        end

        % -----------------------------------------------------------------

        fprintf(1, '\nLoading volumes.\t');

        v = spm_read_vols(V);
        num_images = size(v, 4);

        if length(varargin) > 1
            mP = varargin{2};
            fprintf(1, 'masking volumes.\t')
            % ----------------------------------------------------------------------------------
            % * load mask, and mask all subjects' contrast values
            % so that we show only those voxels existing for all subjects.
            % ----------------------------------------------------------------------------------

            vm = spm_read_vols(spm_vol(mP));
            tmp = size(v);
            if any(size(vm) - tmp(1:3)),
                fprintf(1, '(reslicing mask).\t')
                [tmp, mP] = reslice_imgs(deblank(hP(1,:)), mP);
                vm = spm_read_vols(spm_vol(mP));
            end

            vm = (vm ~= 0);
            vm = +vm;
            vm(vm == 0) = NaN;

            %attempting something less memory-intensive - commented out line
            %required 3 times the size of v in memory... approx 1.2 gigs
            %v = v .* repmat(vm, [1 1 1 size(v, 4)]);
            for i=1:num_images
                v(:,:,:,i) = v(:,:,:,i) .* vm;
            end

            if length(varargin) > 6
                rimg = spm_read_vols(spm_vol(refimg));
                tmp = size(v);
                if any(size(rimg) - tmp(1:3)),
                    fprintf(1, '(reslicing ref image).\t')
                    [tmp, refimg] = reslice_imgs(deblank(hP(1,:)), refimg);
                    rimg = spm_read_vols(spm_vol(refimg));
                end
                rimg(rimg == 0) = NaN;
            end
        end

        % ----------------------------------------------------------------------------------
        % * Standardize, if asked for
        % ----------------------------------------------------------------------------------
        if dostd
            fprintf('\nScaling images to mean 0 and var 1.\t')
        end

        for i = 1:num_images
            tmp = v(:,:,:,i);
            tmp = tmp(:);
            mystd(i) = std(tmp(~isnan(tmp)));
            g(i) = mean(tmp(~isnan(tmp) & tmp ~= 0));

            if dostd
                v(:,:,:,i) = (v(:,:,:,i) - g(i)) ./ mystd(i);
            end
        end
        clear tmp

        % Metrics
        [d, d2, g, mystd] = get_dist(v);

        c = [];
        c2 = [];
        mi = [];
        if exist('rimg', 'var')
            [c, c2, mi] = get_correl(v, do_corr, do_mi, rimg);
        else
            [c, c2, mi] = get_correl(v, do_corr, do_mi);
        end

        if(do_chi2)
            [b, eigv, eigval] = pc(c2, doplot);
            [ds, S, p] = multivar_dist(b);
        else
            ds = [];
            S = [];
            p = [];
            b = [];
            eigv = [];
            eigval = [];
        end

        % plotting stuff

        if doplot
            plot_results();
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inline Functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function plot_results()
        figure('Color', 'w'); subplot(2, 2, 1);
        bar(g);
        xlabel('Image');
        ylabel('Global for each image');
        hold on;
        plot([0 num_images], [mean(g) mean(g)], 'r-');
        if length(varargin) > 3
            title(mytitle);
        end

        subplot(2, 2, 2);
        bar(mi);
        xlabel('Image');
        if length(varargin) < 7
            varargin{7} = 'refimg';
        end

        if exist('rimg', 'var')
            ylabel('Mutual info with reference');
            title(varargin{7});
        else
            ylabel('Mutual info with mean');
        end

        subplot(2, 2, 3);
        imagesc(d2);
        colormap hot;
        xlabel('Image');
        ylabel('Image');
        title('Avg abs dist');

        if isempty(textlab)
            for j = 1:size(b, 1)
                textlab{j} = num2str(j);
            end
        end

        % mds-like (pca version) on similarities (correlations)
        subplot(2, 2, 4); hold on
        if size(b, 2) > 1
            plot(b(:,1), b(:,2), 'Color', 'w');
            xlabel('Component 1');
            ylabel('Component 2');
            for j = 1:size(b, 1)
                text(b(j, 1), b(j, 2), textlab{j}, 'Color', 'b');
            end
        else
            plot(ones(1, length(b)), b, 'Color', 'w');
            xlabel('Component 1');
            for j = 1:size(b, 1)
                text(1, b(j), textlab{j}, 'Color', 'b');
            end
            set(gca, 'XTick', [-1 1]);
        end

        title('MDS of global image values');
    end
end




function [d, d2, g, mystd] = get_dist(v)
    % get the distances from the average
    fprintf('getting distances from mean.\t')
    gmn = mean(v, 4);
    d2 = zeros(size(v, 4), size(v, 4));

    for i = 1:size(v, 4)
        dd = abs((v(:,:,:,i) - gmn));
        d(i) = nanmean(dd(:));
        gg = v(:,:,:,i);
        gg = gg(~isnan(gg));
        g(i) = mean(gg(:));
        mystd(i) = std(gg(:));

        % get the distances from all other images
        d2(i, i) = 0;
        for j = (i+1):size(v, 4)
            dd = abs((v(:,:,:,i) - v(:,:,:,j)));
            d2(i, j) = mean(dd(~isnan(dd)));
        end
    end

    %fill in lower triangle
    d2 = d2 + d2';
    if any(isnan(g)), disp('Warning!  Some images have no valid values.'); end
end

function [c, c2, mi] = get_correl(v, do_corr, do_mi, varargin)
    % get correlations with the average and with all
    fprintf(1, 'getting correlations.\t')

    if length(varargin) > 0  %#ok  % we have a ref image instead of the grand mean
        rimg = varargin{1};
        gv = rimg(:);
    end

    if ~exist('gv', 'var') || isempty(gv)
        gmn = mean(v, 4);
        gv = gmn(:);
    end

    szv = size(v);
    v = reshape(v, prod(szv(1:3)), szv(4));

    % eliminate NaNs
    wh = union(find(isnan(gv)), find(any(isnan(v), 2)));
    if ~isempty(wh)
        gv(wh,:) = [];
        v(wh,:) = [];
    end

    c = [];
    mi = [];
    for i = 1:size(v, 2),
        if(do_corr)
            cc = corrcoef(gv, v(:,i));
            c(i) = cc(1, 2);
        end

        % mutual information
        if(do_mi)
            fprintf(1, 'MI.')
            [H, mi(i)] = hist2(gv, v(:,i), 64);
        end
    end

    c2 = corrcoef(v);
end


function [b, v, d] = pc(a, doplot)
    % a is original matrix, b is principal components, v is eigenvectors
    % (weights on columns, which = weights on voxels)

    [v, d] = eig(a);
    b = (pinv(v) * a')' ./ repmat((diag(d)').^.5, size(a, 1), 1);
    b = fliplr(b);
    v = fliplr(v);

    num = min(10, sum(diag(d) >= 1));
    b = b(:,1:num);
    v = v(:,1:num);
    origd = diag(d);
    d = diag(d)';
    d = fliplr(d);

    if doplot
        figure('Color', 'w');
        bar(real(d.^.5));
        cumprct = cumsum(d ./ sum(d));
        for i = 1:length(d)
            str = sprintf('%3.2f', cumprct(i));
            text(i-.2, d(i)^.5+.2, str);
        end
        title('Eigenvalues');
    end

    d = d(1:num);

    if num == 0
        warning('No eigenvalues above 1!');
        disp(origd);
    end
end


