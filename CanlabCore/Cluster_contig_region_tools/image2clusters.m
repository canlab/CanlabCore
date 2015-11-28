function varargout = image2clusters(varargin)
% Menu-driven function for getting clusters from an image file (e.g., a
% t-image)
%
% :Usage:
% ::
%
%    cl = image2clusters([overlay image name])
%
% Can also return clusters active in two contrasts, sorted by increases in
% both, decreases in both, inc in first, dec in first
% useful for testing whether something is both activated and correlated!
% e.g., see active_plus_corr_scatterplot_plugin
%
% % :Example:
% ::
%
%    [pospos,negneg,posneg,negpos] = image2clusters(overlay)
%
% ..
%    NOTE - documented usage does not appear to correspond to function
%    behavior - comment added by Jared 7/13/06
% ..


overlay = which('scalped_single_subj_T1.img');
if length(varargin) > 0 && ~isempty(varargin{1}), overlay = varargin{1};, end
if isempty(overlay), warning('Cannot find overlay image on path!');, end


dd = spm_get(1,'*img','Select image file to get clusters from.',pwd);

ist = spm_input(['Is the image a t-image (1 or 0)? '],[],'y/n',[1 0],1);


if ist,
    thr = spm_input(['Enter p-value threshold (uncorrected): ']);
    df = spm_input(['Enter degrees of freedom (N - params estimated): ']);
    
    u = tinv(1-thr,df);
else
    u = spm_input(['Enter threshold value for image: ']);
end

kthr = spm_input(['Enter extent threshold in voxels: ']);

%fx = spm_input(['Get which effects?  Type pos, neg, or both: '],[],'s','both');
fx = spm_input(['Get which effects? '],[],'b',{'pos' 'neg' 'both'},{'pos' 'neg' 'both'}); fx = fx{1};

cross = spm_input(['Mask with another t-image (e.g., correlation t-values): '],[],'y/n',[1 0],0);

[P2,P,sigmat,sigmatneg] = threshold_imgs(dd,u,kthr,fx);
disp(['Thresholding image: created ' num2str(P2)])
    
if ~cross, 
    % go ahead
    
    cl = mask2clusters(P2);

    doplot = spm_input(['Display orthviews? '],[],'y/n',[1 0],0);

    if doplot
        switch fx
            case 'both'
                cluster_orthviews(cl,'bivalent','overlay',overlay);
            otherwise

                color = spm_input(['Enter vector of colors [r g b], e.g., [1 0 0] for red: ']);
                cluster_orthviews(cl,{color},'overlay',overlay);
        end
    end

    varargout{1} = cl;
    
else
    % cross significant regions in first image with 2nd image
    dd2 = spm_get(1,'*img','Select second t-image (e.g., correlation t-vals).',pwd);
    thr = spm_input(['Enter p-value threshold (uncorrected): ']);
    df = spm_input(['Enter degrees of freedom (N - params estimated): ']);
    u = tinv(1-thr,df);
    kthr = spm_input(['Enter extent threshold in voxels: ']);
    %fx = spm_input(['Get which effects?  Type pos, neg, or both: '],[],'s','both');
    fx = spm_input(['Get which effects? '],[],'b',{'pos' 'neg' 'both'},{'pos' 'neg' 'both'}); fx = fx{1};

    [P3,tmp,sigmat,sigmatneg] = threshold_imgs(dd2,u,kthr,fx);
    disp(['Thresholding image: created ' num2str(P3)])
    
    V1 = spm_vol(P2);   % first image
    v1 = spm_read_vols(V1);
    
    V2 = spm_vol(P3);   % second image
    v2 = spm_read_vols(V2); 
    
    pospos = mask2clusters(v1>0 & v2>0, V1.mat); if ~isempty(pospos),pospos(1).title = ['Positive in ' P2 ' and ' P3];,end
    posneg = mask2clusters(v1>0 & v2<0, V1.mat); if ~isempty(posneg),posneg(1).title = ['Positive:' P2 ' and negative:' P3];,end
    negneg = mask2clusters(v1<0 & v2<0, V1.mat); if ~isempty(negneg),negneg(1).title = ['Negative in ' P2 ' and ' P3];,end
    negpos = mask2clusters(v1<0 & v2>0, V1.mat); if ~isempty(negpos),negpos(1).title = ['negative:' P2 ' and positive:' P3];,end

     cluster_orthviews(pospos,{[1 0 0]},'overlay',overlay);
     cluster_orthviews(negneg,{[0 0 1]},'overlay',overlay,'add');
     cluster_orthviews(posneg,{[1 .5 0]},'overlay',overlay,'add');
     cluster_orthviews(negpos,{[0 .5 1]},'overlay',overlay,'add');
     
     % make legend
     names = {'Activated and + correlation' 'Deactivated and - correlation' 'Activated and - correlation' 'Deactivated and + correlation'};
    colors = [{[1 0 0]} {[0 0 1]} {[1 .5 0]} {[0 .5 1]}]
    
    curfig=gcf;
    figure;set(gcf,'Color','w');set(gca,'FontSize',16); 
    makelegend(names,colors,1);
    figure(curfig);

    cl = pospos; 
    cl = merge_clusters(cl,negneg);
    cl = merge_clusters(cl,posneg);
    cl = merge_clusters(cl,negpos);

    varargout{1} = cl;
    varargout{2} = pospos;
    %varargout{3} = negneg;
    varargout{4} = posneg;
    varargout{5} = negpos;
end

return

