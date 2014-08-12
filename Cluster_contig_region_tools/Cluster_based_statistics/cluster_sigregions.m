function [cl,clpos,clneg,clrpos,clrneg] = cluster_sigregions(cl,pthr,varargin)
% [cl,clpos,clneg,clrpos,clrneg] = cluster_ttest(clusters,pthr,[behavioral regressor])
% tor wager  Oct 9, 2003
%
% does t-tests on cluster timeseries and all voxels
% which should contain group data (e.g., con* data)
% 
% if optional behavioral regressor is entered,
% computes cluster extent nonparametric p-values, and reports
% extent npm p-values for intercept and behavioral reg.
% 
% writes output to cluster_ttest_output.txt

fid = fopen('cluster_ttest_output.txt','a');

if length(varargin) > 0, beh = varargin{1};,end

for i = 1:length(cl)
    
% df = n - k    primary threshold
tthr1 = abs(tinv(pthr,size(cl(i).all_data,1) - size(cl(i).ttest.t,1)));

% ----------------------------------------
% * do nonparametric permutation 
% ----------------------------------------

if ~isfield(cl(i),'npm_ttest'),cl(i).npm_ttest = [];,end
if isfield(cl(i).npm_ttest,'tmax') & isfield(cl(i).npm_ttest,'textmax')
    tthr = cl(i).npm_ttest.tthr ;
    rthr = cl(i).npm_ttest.rthr ;
    tnegthr = cl(i).npm_ttest.tnegthr ;
    rnegthr = cl(i).npm_ttest.rnegthr ;
    tmax = cl(i).npm_ttest.tmax ;
    tmin = cl(i).npm_ttest.tmin ;
    rmax = cl(i).npm_ttest.rmax ;
    rmin = cl(i).npm_ttest.rmin ;
    textmax = cl(i).npm_ttest.textmax ;
    rextmax = cl(i).npm_ttest.rextmax ;
else
    
[tthr,rthr,tnegthr,rnegthr,tmax,rmax,tmin,rmin,textmax,rextmax] = ...
    npm_ttest(cl(i).all_data,1000,beh,pthr,cl(i).XYZ);
close

cl(i).npm_ttest.tthr = tthr;
cl(i).npm_ttest.rthr = rthr;
cl(i).npm_ttest.tnegthr = tnegthr;
cl(i).npm_ttest.rnegthr = rnegthr;
cl(i).npm_ttest.tmax = tmax;
cl(i).npm_ttest.rmax = rmax;
cl(i).npm_ttest.textmax = textmax;
cl(i).npm_ttest.rextmax = rextmax;
cl(i).npm_ttest.tmin = tmin;
cl(i).npm_ttest.rmin = rmin;

end

fprintf(fid,'\nCl.\tVoxels\tThreshold\tMask\tData\t\n');
fprintf(fid,'\tx\ty\tz\tVoxels\tCluster corr. p.\tMax. t\tCorrected p\tDirection of effect\tr\n');
fprintf(fid,'%3.0f\t%3.0f\t%3.4f\t%s\t%s\t\n',i,cl(i).numVox,pthr,cl(i).P,cl(i).imP(1,:));

% ----------------------------------------
% * get sig clusters
% ----------------------------------------

clpos(i).from_cluster = i;
clpos(i).M = cl(i).M;
clpos(i).voxSize = cl(i).voxSize;
clpos(i).threshold = tthr1;
clpos(i).title = 'Positive effects';

clrpos(i).from_cluster = i;
clrpos(i).M = cl(i).M;
clrpos(i).voxSize = cl(i).voxSize;
clrpos(i).threshold = tthr1;
clrpos(i).title = 'Positive correlations';

sig = cl(i).ttest.t > tthr1;

if any(sig(1,:))
    clpos(i).XYZ = cl(i).XYZ(:,find(sig(1,:))); 
    clpos(i).XYZmm = cl(i).XYZmm(:,find(sig(1,:)));
    clpos(i).Z = cl(i).Z(:,find(sig(1,:))); 
    clpos(i).t = cl(i).ttest.t(:,find(sig(1,:))); 
    gopos(i) = 1;
    
    % print table entry
    wcl = spm_clusters(clpos(i).XYZ);
    clpos(i).clusters = wcl;
    for j = 1:max(wcl)
     
        wh = find(wcl==j);
        
        % corrected p
    
        xyz = mean(clpos(i).XYZmm(:,wh),2)';
        if length(xyz) == 1, xyz = clpos(i).XYZmm(:,wh)';,end
        nv = length(wh);
        
        if length(wh)>1,max_t = max(clpos(i).t(:,wh)'); max_t = max_t(1);,
        else, max_t = clpos(i).t(1,wh);
        end
        
        clpos(i).center = xyz;
    
        extcor_p = 1 - (sum(textmax <= nv) ./ length(textmax));
        max_p = 1 - (sum(tmax <= max_t) ./ length(tmax));

        fprintf(fid,'\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.4f\t%3.2f\t%3.4f\t+\n', ...
            xyz(1),xyz(2),xyz(3),nv,extcor_p,max_t,max_p)
    end
        
else
    gopos(i) = 0;
end

if any(sig(2,:))
    clrpos(i).XYZ = cl(i).XYZ(:,find(sig(2,:))); 
    clrpos(i).XYZmm = cl(i).XYZmm(:,find(sig(2,:)));
    clrpos(i).Z = cl(i).Z(:,find(sig(2,:))); 
    clrpos(i).t = cl(i).ttest.t(:,find(sig(2,:))); 
    gorpos(i) = 1;
    
    % print table entry
    wcl = spm_clusters(clrpos(i).XYZ);
    clrpos(i).clusters = wcl;
    for j = 1:max(wcl)
     
        wh = find(wcl==j);
        
        % corrected p
    
        xyz = mean(clrpos(i).XYZmm(:,wh),2)';
        if length(xyz) == 1, xyz = clrpos(i).XYZmm(:,wh)';,end
        nv = length(wh);
        
        tmp = cl(i).all_data(:,find(sig(2,:)));
        tmp = tmp(:,wh);
        tmp = corrcoef([tmp beh]);
        tmp = tmp(end,1:end-1);
        max_r = max(tmp);
        if max_r > .99, warning('problem?'),keyboard,end
    
        if length(wh)>1,max_t = max(clrpos(i).t(:,wh)'); max_t = max_t(2);,
        else, max_t = clrpos(i).t(2,wh);
        end
        clrpos(i).center = xyz;
    
        extcor_p = 1 - sum(rextmax <= nv) ./ length(rextmax);
        max_p = 1 - sum(rmax <= max_t) ./ length(rmax);

        fprintf(fid,'\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.4f\t%3.2f\t%3.4f\t+\t%3.2f\n', ...
            xyz(1),xyz(2),xyz(3),nv,extcor_p,max_t,max_p,max_r);
    end
    
else
    gorpos(i) = 0;
end

% negative

clneg(i).from_cluster = i;
clneg(i).M = cl(i).M;
clneg(i).voxSize = cl(i).voxSize;
clneg(i).threshold = tthr1;
clneg(i).title = 'Negative effects';

clrneg(i).from_cluster = i;
clrneg(i).M = cl(i).M;
clrneg(i).voxSize = cl(i).voxSize;
clrneg(i).threshold = tthr1;
clrneg(i).title = 'Negative correlations';

sig = cl(i).ttest.t < -tthr1;

if any(sig(1,:))
    clneg(i).XYZ = cl(i).XYZ(:,find(sig(1,:))); 
    clneg(i).XYZmm = cl(i).XYZmm(:,find(sig(1,:)));
    clneg(i).Z = cl(i).Z(:,find(sig(1,:))); 
    clneg(i).t = cl(i).ttest.t(:,find(sig(1,:))); 
    goneg(i) = 1;
        
    % print table entry
    wcl = spm_clusters(clneg(i).XYZ);
    clneg(i).clusters = wcl;
    for j = 1:max(wcl)
     
        wh = find(wcl==j);
        
        % corrected p
    
        xyz = mean(clneg(i).XYZmm(:,wh),2)';
        if length(xyz) == 1, xyz = clpos(i).XYZmm(:,wh)';,end
        nv = length(wh);
        
        if length(wh)>1,max_t = min(clneg(i).t(:,wh)'); max_t = max_t(1);,
        else, max_t = clneg(i).t(1,wh);
        end
        
        clneg(i).center = xyz;
    
        extcor_p = 1 - sum(textmax <= nv) ./ length(textmax);
        max_p = 1 - sum(tmin >= max_t) ./ length(tmin);

        fprintf(fid,'\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.4f\t%3.2f\t%3.4f\t-\n', ...
            xyz(1),xyz(2),xyz(3),nv,extcor_p,max_t,max_p);
    end
        
else
    goneg(i) = 0;
end

if any(sig(2,:))
    clrneg(i).XYZ = cl(i).XYZ(:,find(sig(2,:))); 
    clrneg(i).XYZmm = cl(i).XYZmm(:,find(sig(2,:)));
    clrneg(i).Z = cl(i).Z(:,find(sig(2,:))); 
    clrneg(i).t = cl(i).ttest.t(:,find(sig(2,:))); 
    gorneg(i) = 1;
    
    % print table entry
    wcl = spm_clusters(clrneg(i).XYZ);
    clrneg(i).clusters = wcl;
    for j = 1:max(wcl)
     
        wh = find(wcl==j);
        
        % corrected p
    
        xyz = mean(clrneg(i).XYZmm(:,wh),2)';
        if length(xyz) == 1, xyz = clrneg(i).XYZmm(:,wh)';,end
        nv = length(wh);
   
        if length(wh)>1,max_t = min(clrneg(i).t(:,wh)'); max_t = max_t(2);,
        else, max_t = clrneg(i).t(2,wh);
        end
        
        tmp = cl(i).all_data(:,find(sig(2,:)));
        tmp = tmp(:,wh);
        tmp = corrcoef([tmp beh]);
        tmp = tmp(end,1:end-1);
        max_r = min(tmp);
        if max_r > .99, warning('problem?'),keyboard,end
    
        clrneg(i).center = xyz;
    
        extcor_p = 1 - sum(rextmax <= nv) ./ length(rextmax);
        max_p = 1 - sum(rmin >= max_t) ./ length(rmin);

        fprintf(fid,'\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.4f\t%3.2f\t%3.4f\t-\t%3.2f\t\n', ...
            xyz(1),xyz(2),xyz(3),nv,extcor_p,max_t,max_p,max_r);
    end
    
else
    gorneg(i) = 0;
end

end % loop through clusters

fprintf(fid,'\n')


% Make montage
try
    
if ~isempty(clpos) & ~isempty(clneg)
    montage_clusters([],cl(unique(cat(2,[clpos.from_cluster clneg.from_cluster]))),clpos,clneg,{'k' 'r' 'b'},'nooverlap')
    %hh=get(gcf,'Children');
    %for i = 1:length(hh),axes(hh(i));, h = findobj(gca,'FaceColor',[0 0 0]);, set(h,'FaceAlpha',.5) ,end
elseif ~isempty(clpos)
    montage_clusters([],cl(cat(2,clpos.from_cluster)),clpos,{'k' 'r'},'nooverlap')
elseif ~isempty(clneg)
    montage_clusters([],cl,clneg,{'k' 'b'},'nooverlap')
else
    disp('No main contrast effects')
end

if ~isempty(clrpos) & ~isempty(clrneg)
    montage_clusters([],cl(unique(cat(2,[clrpos.from_cluster clrneg.from_cluster]))),clrpos,clrneg,{'k' 'r' 'b'},'nooverlap')
elseif ~isempty(clrpos)
    montage_clusters([],cl(cat(2,clrpos.from_cluster)),clrpos,{'k' 'r'},'nooverlap')
elseif ~isempty(clrneg)
    montage_clusters([],cl(cat(2,clrneg.from_cluster)),clrneg,{'k' 'b'},'nooverlap')
else
    disp('No correlation effects')
end

catch
    disp('error with montage')
end

fclose(fid);

return



