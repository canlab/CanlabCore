% example: get s
%
% load ROIfile
% for each voxel:
%   - remove session means
%   - HP filter
%   - convert to percent change
%   - fit model
%   - save fits and residual


/Users/tor/Documents/Tor_Documents/CurrentExperiments/VNL/ROIFILES
load ../vnl_condf.mat

hrf = hrf ./ max(hrf);
[LinearX] = getPredictors(vnl_condf,hrf);
LinearX = LinearX - repmat(mean(LinearX),size(LinearX,1),1);

xtxitx = pinv(LinearX);

D = dir('*vis.mat');
b_by_ymean = [];

for s = 1:length(D)

    % load the file
    load(D(s).name)
    eval(['ROI = ' D(s).name(1:end-4)])
    
    clear b, clear sigma
    
for i = 1:size(ROI.ts.indiv,2), 
    
        if length(ROI.adjustedy) == 6020
			[y,ymean] = voistat('adjusty',ROI.ts.indiv(:,i),41,.5,7);
		else
			[y,ymean] = voistat('adjusty',ROI.ts.indiv(:,i),41,.5,8);
		end
        
    b(:,i) = xtxitx(:,1:length(y)) * y; 
    sigma(i) = std(y - LinearX(1:length(y),:) * b(:,i));,
    
end

%mean(b')
%sigma

%btmp = xtxitx * ROI.adjustedy
%figure; plot(ROI.adjustedy); hold on; plot(LinearX * btmp,'r')

    B(s,:) = mean(b');
    SIGMA(s) = mean(sigma);
    SNR(s,:) = mean(b ./ repmat(sigma,6,1),2)';
    
end


save vnl_snr_out B SIGMA SNR b_by_ymean D

B(1:2,:) = [];
sigma(1:2) = [];
SNR(1:2,:) = [];

% mean B, 1 2 5 6 10 11 flashes
% 0.8258    1.1419    1.2200    1.2461    0.8393    0.9499

% mean SNR
% 0.4724    0.6557    0.7084    0.7278    0.5038    0.5574

% mean sigma
% 1.6684    percent change


% to get betas from UNSCALED data, to assess fits of betas by mean
% level...does the signal scale with baseline levels?


b_by_ymean = [];

for s = 3:length(D)

    % load the file
    load(D(s).name)
    eval(['ROI = ' D(s).name(1:end-4)])
    
    clear b, clear sigma
    
for i = 1:size(ROI.ts.indiv,2), 
    
    %   ALTER VOISTAT TO GET RID OF SCALING BEFORE DOING THIS
        if length(ROI.adjustedy) == 6020
			[y,ymean] = voistat('adjusty',ROI.ts.indiv(:,i),41,.5,7);
		else
			[y,ymean] = voistat('adjusty',ROI.ts.indiv(:,i),41,.5,8);
		end
        
    b(:,i) = xtxitx(:,1:length(y)) * y; 
    sigma(i) = std(y - LinearX(1:length(y),:) * b(:,i));,
    
    b_by_ymean(end+1,:) = [b(:,i)' ymean s];
    
end

end

figure;gscatter(b_by_ymean(:,7),b_by_ymean(:,1),b_by_ymean(:,8),'rgbycmk','o+xv^s');
t1 = unique(b_by_ymean(:,8));
for i = t1'
    
    tmp = b_by_ymean(b_by_ymean(:,8) == i,[1 7]);
    t2 = corrcoef(tmp);  b_by_mean_r(i) = t2(1,2);
    bb = polyfit(tmp(:,2),tmp(:,1),1); h = refcurve(bb);
    set(h,'Color',rand(1,3))
end

set(gca,'FontSize',16); xlabel('Mean scanner signal'),ylabel('Estimated response magnitude'),title('Estimated magnitude as a function of mean signal')





