function clusters = cluster_ttest(clusters,varargin)
% clusters = cluster_ttest(clusters,[covariate])
% tor wager  Oct 9, 2003
%
% does t-tests on cluster timeseries and all voxels
% which should contain group data (e.g., con* data)
% 
% 
fid = fopen('cluster_ttest_output.txt','a');

fprintf(fid,'\nT-test results for clusters\n');
fprintf(fid,'Cluster\tOverall t\tOverall p\tNum vox\tMax t\tMin t\tBonf_thresh\tBonf_posvox\tBonf_negvox');

fprintf(fid,'\tNPM_posthresh\tNPM_negthresh\tNPM_posvox\tNPM_negvox\tpos_NPM_corp\tneg_NPM_corp');

if length(varargin) > 0, 
    beh = varargin{1};,
    if size(beh,2) ~= 1, beh = beh';,end
    %X = [beh ones(length(beh),1)];
    X = beh - mean(beh);
    X(isnan(beh),:) = [];
    clusters(1).covariate = beh;
    fprintf(fid,'\tOverall r\tt\tp\tMax r\tMin r\tBonf_r_pos\tBonf_r_neg');
    fprintf(fid,'\tNPM_r+thresh\tNPM_r-thresh\tNPM_r_posvox\tNPM_r_negvox\tposr_NPM_corp\tnegr_NPM_corp');
    
else,beh = [];,X = []; % define for eliminating NaN voxels when no beh vector
end



fprintf(fid,'\n');

for i = 1:length(clusters)
    
    
    if ~isempty(beh)
        
            % simple regression
            [b,dev,stat] = glmfit(X,clusters(i).timeseries);
            if ~isfield(stat,'p'),stat.p = NaN;,stat.t = NaN;,end
            stats.tstat = stat.t;
            stats.p = stat.p;
            stats.beta = stat.beta;
            stats.df = stat.dfe;
            tmp = corrcoef(X,clusters(i).timeseries);
            roverall = tmp(1,2);
            
        else
            % ttest
            [h,p,ci,stats]=ttest(clusters(i).timeseries);
            stats.p = p;
        end
        
    clusters(i).overall_ttest = stats;
    
    fprintf(fid,'%3.0f\t%3.2f\t%3.4f\t%3.0f',i,stats.tstat(1),stats.p(1),clusters(i).numVox);
    
    % Bonferroni threshold
    bfp = .05 ./ clusters(i).numVox;
    bf = tinv(1-(bfp),stats.df);
    
    clear h, clear p, clear t, clear r,
    for j = 1:size(clusters(i).all_data,2)
        
        if any(isnan(clusters(i).all_data(:,j)))
            p(:,j) = zeros(size(X,2)+1,1) * NaN;,t(:,j) = zeros(size(X,2)+1,1) * NaN;
            h(:,j) = zeros(size(X,2)+1,1); r(j) = NaN;
        elseif ~isempty(beh)
        
            % simple regression
            [b,dev,stat] = glmfit(X,clusters(i).all_data(:,j));
            if ~isfield(stat,'p'),stat.p = zeros(size(X,2)+1,1) * NaN;,stat.t = zeros(size(X,2)+1,1) * NaN;,end
            h(:,j) = stat.p < bfp;
            p(:,j) = stat.p;
            t(:,j) = stat.t;
            tmp = corrcoef(X,clusters(i).all_data(:,j));
            r(j) = tmp(1,2);
            
        else
            
            [h(j),p(j),ci,stats]=ttest(clusters(i).all_data(:,j),0,bfp);
            t(j) = stats.tstat;
        end
    
    
    end
    
    mns = nanmean(clusters(i).all_data);
    pos = h(1,:) & (mns > 0);
    neg = h(1,:) & (mns < 0);
    
    if any(pos), clusters(i).ttest.poscenter = mean(clusters(i).XYZmm(:,find(pos)),2)';, end
    if any(neg), clusters(i).ttest.negcenter = mean(clusters(i).XYZmm(:,find(neg)),2)';, end
    
    if ~isempty(beh) & size(h,1) > 1
        poscor = h(2,:) & (t(2,:) > 0);
        negcor = h(2,:) & (t(2,:) < 0);
        clusters(i).ttest.poscor = poscor;
        clusters(i).ttest.negcor = negcor;
        if any(poscor), clusters(i).ttest.posrcenter = mean(clusters(i).XYZmm(:,find(poscor)),2)';, end
        if any(negcor), clusters(i).ttest.negrcenter = mean(clusters(i).XYZmm(:,find(negcor)),2)';, end
    end
    
    tm = max(t(1,:));
    
    % Max t\tBonf_thresh\tBonf_posvox\tBonf_negvox')
    fprintf(fid,'\t%3.2f\t%3.2f\t%3.2f\t%3.0f\t%3.0f',tm,min(t(1,:)),bf,sum(pos),sum(neg));

    % NPM simulation
    [tthr,rthr,tnegthr,rnegthr,tmax,rmax,tmin,rmin,textmax,rextmax] = ...
        npm_ttest(clusters(i).all_data,1000,beh,.005,clusters(i).XYZ);
    close
    
    tcorp = 1 - (sum(tmax <= tm) ./ length(tmax));
    tcorpneg = 1 - (sum(tmin >= tm) ./ length(tmin));
    clusters(i).npm_ttest.tcorp = tcorp;
    clusters(i).npm_ttest.tcorpneg = tcorpneg;
    clusters(i).npm_ttest.tmax = tmax;
    clusters(i).npm_ttest.tmin = tmin;
    
    if any(any(~isnan(t)))
    
    fprintf(fid,'\t%3.2f\t%3.2f\t%3.0f\t%3.0f\t%3.4f\t%3.4f', ...
        tthr,tnegthr,sum(t(1,:)>tthr),sum(t(1,:)<tnegthr),tcorp,tcorpneg);
    
    
    clusters(i).npm_ttest.tthr = tthr;
    clusters(i).npm_ttest.rthr = rthr;
    clusters(i).npm_ttest.tnegthr = tnegthr;
    clusters(i).npm_ttest.rnegthr = rnegthr;
    clusters(i).npm_ttest.tmax = tmax;
    clusters(i).npm_ttest.tmin = tmin;
    clusters(i).npm_ttest.rmax = rmax;
    clusters(i).npm_ttest.rmin = rmin;
    clusters(i).npm_ttest.textmax = textmax;
    clusters(i).npm_ttest.rextmax = rextmax;
    
    if any(t(1,:)>tthr), clusters(i).npm_ttest.poscenter = mean(clusters(i).XYZmm(:,find(t(1,:)>tthr)),2)';, end
    if ~isempty(rthr), if any(t(2,:)>rthr), clusters(i).npm_ttest.posrcenter = mean(clusters(i).XYZmm(:,find(t(2,:)>rthr)),2)';, end, end
    if any(t(1,:)<tnegthr), clusters(i).npm_ttest.negcenter = mean(clusters(i).XYZmm(:,find(t(1,:)<tnegthr)),2)';, end
    if ~isempty(rnegthr), if any(t(2,:)<rnegthr), clusters(i).npm_ttest.negrcenter = mean(clusters(i).XYZmm(:,find(t(2,:)<rnegthr)),2)';, end,end
    
    if ~isempty(beh)
        %fprintf(fid,'\tOverall r\tt\tp\tMax r\tMin r\tBonf_r_pos\tBonf_r_neg\t')
        fprintf(fid,'\t%3.2f\t%3.2f\t%3.4f\t%3.2f\t%3.2f\t%3.0f\t%3.0f', ...
            roverall,stats.tstat(end),stats.p(end),max(r),min(r),sum(poscor),sum(negcor));
        
        if ~isempty(rthr)
            rm = max(t(2,:));
            rcorp = 1 - (sum(rmax <= rm) ./ length(rmax));
            rcorpneg = 1 - (sum(rmin >= tm) ./ length(rmin));
            clusters(i).npm_ttest.rcorp = rcorp;
            clusters(i).npm_ttest.rcorpneg = rcorpneg;
            clusters(i).npm_ttest.rmax = rmax;
            clusters(i).npm_ttest.rmin = rmin;
    
            fprintf(fid,'\t%3.2f\t%3.2f\t%3.0f\t%3.0f\t%3.4f\t%3.4f\t', ...
                rthr,rnegthr,sum(t(2,:)>rthr),sum(t(2,:)<rnegthr),rcorp,rcorpneg);
        end
    end
    fprintf(fid,'\n');
    
    clusters(i).ttest.bf_threshold = bf;
    clusters(i).ttest.numSig = sum(h);
    clusters(i).ttest.sig = h;
    clusters(i).ttest.p = p;
    clusters(i).ttest.t = t;
    
    else    % no results because of NaNs
        fprintf(fid,'\n')
    end

end

fclose(fid)

return

    
        
    
    