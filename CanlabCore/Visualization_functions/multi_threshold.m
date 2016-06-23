function cl = multi_threshold(P,type,df,varargin)
% ::
%
%    cl = multi_threshold(P,type,df,[overlay image name])
%       This function wants a statistic image.
%       F contrast:  df = xSPM.df;
%
%       type = 'F' or 'T'or 'none'
%
% ..
%    really good function.
%    tor wager
% ..

% ..
%    defaults
%       Colors:     postive - Red Yellow Orange
%                   negative -  Blue, Light Blue, Aqua
%       Threshold:  [.001 .005 .05] uncorrected
%       K Sizes:    [5 5 10]
%   
%   example
%
%   o2 = multi_threshold(out.t, 'o2', o2, 'thresh', [.005 .01 .05], 'sizethresh', [1 1 1]);
%   o2 = montage(o2, 'coronal', 'slice_range', [-20 20], 'onerow');
%   o2 = addblobs(o2, region(out.t));
% ..

red = [1 0 0]; % Colors
yellow = [1 1 0];
orange = [.9 .5 0];
orange2 = [1 .6 .1];
orange3 = [1 .7 .3];
buff = [.7 .6 .4];
dkred = [.7 .3 .1];

blue = [0 0 1];
ltblue = [0 .5 .9];
ltblue2 = [.2 .7 1];
aqua = [0 .5 .5];
dkblue = [0 .1 .7];
dkblue2 = [0 .2 .5];

% low thresh colors
colors = {yellow orange3 orange2};
colors2= {blue ltblue ltblue2};

% thresholds

thr = [.001 .005 .05];   % Inf = FDR, thr should be increasing p-values
sizes = [5 5 10];     % minimum sizes
direction = 'both';  % 'pos','neg','both'

% masking

maskimg = [];       % empty, or specify image name
%which('scalped_avg152T1_graymatter_smoothed.img');
%maskimg = '/Users/tor/Desktop/Opioid_Placebo2/mean_images/avg_thresh_1-2.img';

% overlay

ovl = [];
if length(varargin) > 0, ovl = varargin{1};,end

% read image
V = spm_vol(P); v = spm_read_vols(V);

% mask image, if necessary
if ~isempty(maskimg)
    disp(['Masking with: ' maskimg])
    VV = spm_vol(maskimg); v2 = spm_read_vols(VV); 
    v = v .* (v2 > 0);
end

% convert to p-values
switch type
    case 'F'
        p = 1 - fcdf(v(:),df(1),df(2));
        direction = 'pos';
        
    case 'T'
        v(isnan(v)) = 0;
        p = 1 - tcdf(v(:),df);
        
    case 'none'
        direction = 'pos';  % for p-image
        % do nothing
        p = v(:);
        %thr = [3 2 1.2];
        colors={[1 1 0] [.9 .5 0] [.7 0 0]};   % pos
        
    otherwise
end

cumwh = logical(zeros(size(p))); % define cumulative significant voxels
cumwh2 = logical(zeros(size(p))); % define cumulative significant voxels

for i = 1:length(thr)
    
    % find significant voxels - only those sig at this thresh, but not
    % higher ones tested previously
    
    if isinf(thr(i)),
        tmp = p; tmp(v==0 | isnan(v))=[];   % eliminate out-of-analysis voxels
        pt = FDR(tmp,.05); 
        
        switch type
            case 'F', 
            disp(['Height thresh: F = ' num2str(finv(1-pt,df(1),df(2))) ', p = ' num2str(pt)])
            case 'T',
             disp(['Height thresh: T = ' num2str(tinv(1-pt,df(1))) ', p = ' num2str(pt)])
        end
        
        
    else
        pt = thr(i);
    end
    if isempty(pt), pt = -Inf; end
    
    switch direction    % applies only for t-values!!
        case 'pos',
            wh = logical(p <= pt); doboth = 0;
        case 'neg'
            wh = logical(p >= 1-pt); doboth = 0;
        case 'both'
            wh = logical(p <= pt);
            wh2 = logical(p >= 1-pt); doboth = 1;
    end
    
    % special for "none" - just threshold -- this would be for t-image, not
    % p-image!
    %if strcmp(type,'none'), wh = logical(p >= pt); doboth = 0;,end
    
    % positive response, or neg only response
    mask{i} = zeros(V.dim(1:3));
    mask{i}(wh) = 1;
    mask{i}(cumwh) = 0; % voxels previously sig get 0
    
    cumwh = cumwh | wh; % update cumulative which sig
            
    % write an image of it
    warning off
    V.fname = 'tmp_mask.img';
    spm_write_vol(V,mask{i});
    warning on
    
    % get clusters
    cl{i} = mask2clusters(V.fname);
    
    
    % negative response (do both)
    if doboth
        mask2{i} = zeros(V.dim(1:3));
        mask2{i}(wh2) = 1;
        mask2{i}(cumwh2) = 0; % voxels previously sig get 0
    
        cumwh2 = cumwh2 | wh2; % update cumulative which sig
        
        
        % write an image of it
        warning off
        V.fname = 'tmp_mask.img';
        spm_write_vol(V,mask2{i});
        warning on
        
        % get clusters
        cl2{i} = mask2clusters(V.fname);
    end
    
    !rm tmp_mask.img
    !rm tmp_mask.hdr
end


% size threshold
for i = 1:length(cl)
    if ~isempty(cl{i}), wh = find(cat(1,cl{i}.numVox) < sizes(i));,else, wh = [];,end
    if ~isempty(wh), cl{i}(wh) = []; end
end

if doboth
    for i = 1:length(cl2)
        if ~isempty(cl2{i}), wh = find(cat(1,cl2{i}.numVox) < sizes(i));,else, wh = [];,end
        if ~isempty(wh), cl2{i}(wh) = []; end
    end
end


% add together, if both
if doboth
    for i = 1:length(thr)   % for legend
        if isempty(cl{i}), saveme(i) = 0; else, saveme(i) = 1;,end
    end
    len = sum(saveme);   % for legend
    
    thr = [thr thr];
    colors = [colors colors2];
    cl = [cl cl2];
end

% remove empties 
for i = 1:length(thr)
    if isempty(cl{i}), saveme(i) = 0; else, saveme(i) = 1;,end
end
saveme = logical(saveme);
cl = cl(saveme);
colors = colors(saveme);
thr = thr(saveme);

if isempty(cl), return, end

% now image clusters
if ~isempty(cl{1}),  cluster_orthviews(cl{1},colors(1),'overlay',ovl);, end

for i = 2:length(cl)
    if ~isempty(cl{i}),cluster_orthviews(cl{i},colors(i),'add');, end
end

% make legend string and legend
for i = 1:length(cl)
    lstr = [];
    if doboth, 
        if i <= len, lstr = ['+ '];,
        else, lstr = ['- '];,
        end
    end
    if isinf(thr(i)), 
        lstr = [lstr 'p < .05 FDR'];,
    else,
        lstr = [lstr 'p < ' num2str(thr(i))];
    end
    legstr{i} = lstr;
end

h = axes('Position',[.55 .2 .25 .2]); hold on; 
set(gca,'FontSize',24)
for i = 1:length(colors),hh(i)=plot(0,0,'Color',colors{i},'LineWidth',10);,end
axis off
legend(legstr)

% check for cases in which we DO NOT want montage
if strcmp(type,'none'), return, end

tmp =[];
for i = 1:length(cl), tmp = [tmp cat(2,cl{i}.XYZ)];,end
tmp = unique(tmp(3,:));
if length(tmp) > 30, stopme = input('More than 30 slices in montage.  Make montage figure? (1/0) ');,if ~stopme, return,end, end

% montage
if length(cl) == 6
    montage_clusters(ovl,cl{1},cl{2},cl{3},cl{4},cl{5},cl{6},colors);
elseif length(cl) == 5
     montage_clusters(ovl,cl{1},cl{2},cl{3},cl{4},cl{5},colors);
elseif length(cl) == 4
     montage_clusters(ovl,cl{1},cl{2},cl{3},cl{4},colors);       
elseif length(cl) == 3
    montage_clusters(ovl,cl{1},cl{2},cl{3},colors);
elseif length(cl) == 2
    montage_clusters(ovl,cl{1},cl{2},colors);
elseif length(cl) == 1
    montage_clusters(ovl,cl{1},colors);
end

cl{1}(1).thr = thr; 
cl{1}(1).colors = colors; 

% montage medial slices
tmp = unique(tmp(1,:));
if length(tmp) > 30, stopme = input('More than 30 slices in medial montage.  Make montage figure? (1/0) ');,if ~stopme, return,end, end

if length(cl) == 6
    montage_clusters_medial(ovl,cl{1},cl{2},cl{3},cl{4},cl{5},cl{6},colors);
elseif length(cl) == 5
     montage_clusters_medial(ovl,cl{1},cl{2},cl{3},cl{4},cl{5},colors);
elseif length(cl) == 4
     montage_clusters_medial(ovl,cl{1},cl{2},cl{3},cl{4},colors);       
elseif length(cl) == 3
    montage_clusters_medial(ovl,cl{1},cl{2},cl{3},colors);
elseif length(cl) == 2
    montage_clusters_medial(ovl,cl{1},cl{2},colors);
elseif length(cl) == 1
    montage_clusters_medial(ovl,cl{1},colors);
end

cl{1}(1).thr = thr; 
cl{1}(1).colors = colors; 

return
