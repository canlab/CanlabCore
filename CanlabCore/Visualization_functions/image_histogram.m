function h = image_histogram(P,varargin)
% :Usage:
% ::
%
%    h = image_histogram(P,[method(string)],[range])
%
% eliminates 0, NaN voxels from either image
% red then blue
%
% :Inputs:
%
%   **P:**
%        is one or two images in string array
%
%   **Methods:**
%        'def'
%
%   **range:**
%        optional 3rd input, range of values to include in histogram
%
% :Example:
% ::
%
%    h = image_histogram('p_Omnibus.img','def',[0 1-eps]);
%
% ..
%    tor wager
% ..

V = spm_vol(P); %v = spm_read_vols(V);

v = iimg_read_vols(V);


if length(varargin) > 0, meth = varargin{1};  else, meth = 'def'; end
if length(varargin) > 1, range = varargin{2};  else, range = [-Inf Inf]; end

for f = 1:size(P,1), [dum,ff{f}] = fileparts(P(f,:)); end


%p2 = which('/Users/tor/Documents/tor_scripts/3DheadUtility/canonical_brains/scalped_avg152T1_graymatter.img')
%[p2,p2new] = reslice_imgs(P(1,:),p2,1);
%maskV = spm_read_vols(spm_vol(p2new)); maskV = maskV(:);

switch meth
    case 'def'
        
        disp('Plotting histograms')
        
        v1 = squeeze(v(:,:,:,1));
        if size(v,4) > 1
            v2 = squeeze(v(:,:,:,2));
        end
        
        v1 = v1(:); 
        if size(v,4) > 1, v2 = v2(:); ,end
        
        if size(v,4) > 1
            wh = isnan(v1) | isnan(v2) | v1 == 0 | v2 == 0 ...
                | v1 < range(1) | v1 > range(2) | v2 < range(1) | v2 > range(2); %| ~maskV;
        else
            wh = isnan(v1) | v1 == 0 | v1 < range(1) | v1 > range(2); 
        end
        
        v1(wh) = []; 
        if size(v,4) > 1
            v2(wh) = [];
        else
            v2 = [];
        end
        
        [h,xx] = hist([v1;v2],max(10,length(v1)./300));
        h1 = hist(v1,xx); 
        
        if size(v,4) > 1, h2 = hist(v2,xx); end
        
    otherwise
        error('unknown action string')
end



figure('Color','w'); 

if size(v,4) > 1, plot(xx,h1,'r',xx,h2,'b'); 
else
        plot(xx,h1,'r');
end
    
set(gca,'FontSize',16)

legend(ff)


return
