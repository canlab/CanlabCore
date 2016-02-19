function [ds,g,mystd,d,d2,c,c2,mi,b,eigv,eigval] = compare_subjects256(varargin)
% This function compares a set of images to one another and does some diagnostics on the similarity among images.
% - It returns multivariate distances and dissimilarities among images
% - It works on the GLOBAL signal after standardizing each image (case 1) or the REGIONAL values in each cluster (case 2) 
% - You can also enter a reference image, in which case each image will be correlated with the ref.
%
% :Usage:
% ::
%
%     function [ds,g,mystd,d,d2,c,c2,mi,b,eigv,eigval] = compare_subjects256([img files or clusters],[mask], ...
%                                    [plot flag],[title on figure],[standardize flag],[text labels],[ref image])
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
% :Outputs: from correls with ref image are in variable "c"
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
%    [ds,g,mystd,d,d2,c,c2,mi] = compare_subjects256(P,which('brain_avg152T1.img'),1,'intext_countloc',1,[],which('avg152T1.img'));
%
% ..
%    by Tor Wager
% ..

doplot = 1; dostd = 0; textlab = [];
if length(varargin) > 2, doplot = varargin{3};  end
if length(varargin) > 3, mytitle = varargin{4};  end
if length(varargin) > 4, dostd = varargin{5};  end
if length(varargin) > 5, textlab = varargin{6};  end
if length(varargin) > 6, refimg = varargin{7};  end

if length(varargin) == 0
    disp('no inputs.'), return
end

if isstruct(varargin{1}), 
    clusters = varargin{1};
    disp('No method implemented yet.')
    return
    
else
    hP = varargin{1};
    
    mypwd = pwd;

    disp([' compare subjects.m Running '])

    fprintf(1,'\nLoading volumes.\t')

    v = spm_read_vols(spm_vol(hP));
    % for i = 1:size(hP,1), v(:,:,:,i) = spm_read_vols(spm_vol(deblank(hP(i,:)))); end
    
    if length(varargin) > 1, mP = varargin{2};  
        fprintf(1,'masking volumes.\t')
        % ----------------------------------------------------------------------------------
        % * load mask, and mask all subjects' contrast values
        % so that we show only those voxels existing for all subjects.
        % ----------------------------------------------------------------------------------

        vm = spm_read_vols(spm_vol(mP));
        tmp = size(v);
        if any(size(vm) - tmp(1:3)),
            fprintf(1,'(reslicing mask).\t')
            [tmp,mP] = reslice_imgs(deblank(hP(1,:)),mP);
            vm = spm_read_vols(spm_vol(mP));
        end
        
        vm = (vm ~= 0); vm = real(vm); vm(vm == 0) = NaN;
        v = v .* repmat(vm,[1 1 1 size(v,4)]);
        
        if length(varargin) > 6,
            rimg = spm_read_vols(spm_vol(refimg));
            tmp = size(v);
            if any(size(rimg) - tmp(1:3)),
                fprintf(1,'(reslicing ref image).\t')
                [tmp,refimg] = reslice_imgs(deblank(hP(1,:)),refimg);
                rimg = spm_read_vols(spm_vol(refimg));
            end
        end
        
        rimg(rimg == 0) = NaN;
            
    end

    % ----------------------------------------------------------------------------------
    % * Standardize, if asked for
    % ----------------------------------------------------------------------------------
    if dostd
        fprintf(1,'\nScaling images to mean 0 and var 1.\t')
        for i = 1:size(v,4)
            tmp = v(:,:,:,i); tmp = tmp(:); mystd = std(tmp(~isnan(tmp)));
            v(:,:,:,i) = (v(:,:,:,i) - mean(tmp(~isnan(tmp)))) ./ mystd;
            
            %chk = v(:,:,:,i); chk = chk(:); chk = chk(~isnan(chk)); mean(chk), std(chk)
        end
    end
    
    % Metrics
    
    [d,d2,g,mystd] = get_dist(v);
    
    if exist('rimg') == 1
        [c,c2,dummy,mi] = get_correl(v,rimg);
    else
        [c,c2,dummy,mi] = get_correl(v);
    end
    
    [b,eigv,eigval] = pc(c2,doplot);
    ds = multivar_dist(b);
    
    % plotting stuff
    
    if length(varargin) > 3, title(mytitle), end
      
    if doplot,
        
        figure('Color','w'); subplot 221
        bar(g), xlabel('Image'),ylabel('Global for each image')
        hold on; plot([0 size(v,4)],[mean(g) mean(g)],'r-')
        if length(varargin) > 3, title(mytitle), end
    
        subplot 222
        bar(mi), xlabel('Image'),
        if exist('rimg')==1, ylabel(['Mutual info with reference']),title(varargin{7}), else, ylabel('Mutual info with mean'),end
    
        subplot 223, imagesc(d2), colormap hot, xlabel('Image'),ylabel('Image'),title('Avg abs dist')
    
        if isempty(textlab), for i = 1:size(b,1), textlab{i} = num2str(i);  end, end
        
        % mds-like (pca version) on similarities (correlations)
        subplot 224, hold on
        if size(b,2) > 1
            plot(b(:,1),b(:,2),'Color','w'); xlabel('Component 1');ylabel('Component 2')
            for i = 1:size(b,1), text(b(i,1),b(i,2),textlab{i},'Color','b'), end
        else
            plot(ones(1,length(b)),b,'Color','w'); xlabel('Component 1');
            for i = 1:size(b,1), text(1,b(i),textlab{i},'Color','b'), end
            set(gca,'XTick',[-1 1])
        end
        
        title('MDS of global image values')
    
        %figure; hold on; plot(b(:,1),b(:,2),'Color','w'); xlabel('Component 1');ylabel('Component 2')
        %for i = 1:size(b,1), text(b(i,1),b(i,2),textlab{i},'Color','b','FontWeight','b'), end
        %title([mytitle ':MDS'],'FontSize',14)
        %figure; bar(mystd),title('Variance of images')
    end

end

return
    
       


function [d,d2,g,mystd] = get_dist(v)

% get the distances from the average
    fprintf(1,'getting distances from mean.\t')
    gmn = mean(v,4);
    for i = 1:size(v,4)
        dd = abs((v(:,:,:,i) - gmn)); dd = dd(:);
        d(i) = mean(dd(~isnan(dd)));
        gg = v(:,:,:,i);  gg = gg(~isnan(gg));
        g(i) = mean(gg(:));
        mystd(i) = std(gg(:));
        
        % get the distances from all other images
        for j = 1:size(v,4)
            if i ~= j
                dd = abs((v(:,:,:,i) - v(:,:,:,j))); dd = dd(:);
                d2(i,j) = mean(dd(~isnan(dd)));
            else
                d2(i,j) = 0;
            end
        end
        
    end
return
    
function [c,c2,vv,mi] = get_correl(v,varargin)

% get correlations with the average and with all
    fprintf(1,'getting correlations.\t')
    
    if length(varargin) > 0,    % we have a ref image instead of the grand mean
        rimg = varargin{1}; 
        gv = rimg(:);
    else
        gmn = mean(v,4); gv = gmn(:);
    end
    
    for i = 1:size(v,4)
        v2 = v(:,:,:,i); vv(:,i) = v2(:);
    end
    
    % eliminate NaNs
    wh = find(isnan(gv));
    if ~isempty(wh), gv(wh,:) = []; vv(wh,:) = [];  end
    wh = find(any(isnan(vv),2));
    if ~isempty(wh), gv(wh,:) = []; vv(wh,:) = [];  end
    
    for i = 1:size(v,4),
        cc = corrcoef(gv,vv(:,i));
        c(i) = cc(1,2);
        
        % plot all points
        %figure('Color','w');plot(vv(:,i),gv,'b.');xlabel('Ref'),ylabel(['Image ' num2str(i)])
        
        % mutual information
        fprintf(1,'MI.')
        [H,mi(i),H2] = hist2(gv,vv(:,i),256);
        
    end
    
    c2 = corrcoef(vv);
    
return


function [b,v,d] = pc(a,doplot)
% a is original matrix, b is principal components, v is eigenvectors 
% (weights on columns, which = weights on voxels)

[v,d]=eig(a);
b = (pinv(v) * a')' ./ repmat((diag(d)').^.5,size(a,1),1);
b = fliplr(b); v = fliplr(v); 

num = min(10,sum(diag(d) >= 1));
b = b(:,1:num); v = v(:,1:num); 
origd = diag(d);
d = diag(d)'; d= fliplr(d); 

if doplot
    figure('Color','w'); bar(real(d.^.5));
    cumprct = cumsum(d ./ sum(d));
    for i = 1:length(d),str=sprintf('%3.2f',cumprct(i)); text(i-.2,d(i)^.5+.2,str),end
    title('Eigenvalues')
end

d = d(1:num);

if num == 0, warning('No eigenvalues above 1!');  origd, end

%figure;plot(b,'r'),hold on;plot(a,'k'), hold on; plot(mean(a,2),'g--'),legend({'eig' 'orig' 'avg'})

return


