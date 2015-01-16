function D = cluster_discrim(cl,beh,varargin)
% D = cluster_discrim(cl,beh,[covs of no interest])
%
%x = cat(2,pain_matrix_early(:).timeseries);
%x = cat(2,pain_matrix_peak(:).timeseries);

diary Cluster_Discrim_Output.txt

disp('cluster_discrim.m: Predictions of behavior with contrast scores');
fprintf(1,'\n-----------------------------------------------------------------\n');

disp(['Input clusters: ' inputname(1)])
disp('Loading data from cl.timeseries');


% Get values from clusters
% ---------------------------------------------
x = cat(2,cl(:).timeseries);
if any(isnan(sum(x))), 
    whnan = find(isnan(sum(x)));
    cl(whnan) = []; x = cat(2,cl(:).timeseries);
    disp(['Found NaNs! Removing clusters ' num2str(whnan)]); 
end



% Remove covariates of no interest
% ---------------------------------------------
covs = [];
if length(varargin) > 0
    covs = varargin{1};
    D.covs = covs;
    disp('Found covariate(s) of no interest; Removing them from behavior and brain data.')
        % partialcorr gets adjusted x, y, and correls; here just use it to
        % get adjusted x
        beh  = partialcor([beh covs],ones(size(beh,1),1),1);
        
    for i = 1:size(x,2)
        x(:,i) = partialcor([x(:,i) covs],ones(size(x,1),1),1);
    end
end


D.data = x;

figure;imagesc(x); colorbar; title('Data for discriminant analysis'); drawnow

D.xyzmm = cat(1,cl(:).mm_center);

y = beh;
D.beh = beh;

y2 = mediansplit(y);
D.beh_medsplit = y2;

if isfield(cl,'imnames'),D.imnames = cl(1).imnames;,end

% Print univariate correlations with regions and location (t-test)
% ---------------------------------------------
disp('Correlations with behavior, after removing covariates of no interest.')
fprintf(1,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n','Cluster','Name','x','y','z','V','r', 'p', 'IRLSr', 'IRLSp','Zavg','Inc/Dec');
warning off
fprintf(1,'\n');
for i = 1:size(x,2)
    [dummy,dummy,r,p,rrob,prob]  = partialcor([beh covs],x(:,i),1);
    
     % print output
     fprintf(1,'Cluster %3.0f\t',i);
     if isfield(cl,'shorttitle'), fprintf(1,'%s\t',cl(i).shorttitle),end
     if isfield(cl,'mm_center'),fprintf(1,'%3.0f\t%3.0f\t%3.0f\t',cl(i).mm_center);,end
     if isfield(cl,'numVox'),fprintf(1,'%3.0f\t',cl(i).numVox);,end            
     fprintf(1,'%3.2f\t%3.2f\t%3.2f\t%3.2f\t',[r p rrob prob])
     
     % t-test
     [hh,p,ci,stats] = ttest(cl(i).timeseries);
     strs = {'Decrease' '--' 'Increase'};
     if ~hh, strs = strs{2};,elseif stats.tstat >0,strs=strs{3};,else,strs = strs{1};end
     fprintf(1,'%3.2f\t%s\t',stats.tstat,strs)
     
     
     fprintf(1,'\n');
end           
warning on    

% --------------------------------------
% get eigenvalues, choose # dims
% --------------------------------------
disp('Computing significant principal components')
[v,score,e] = princomp(x);
[pc,D.pca_npm] = pca_npm(x,1000);
whsig = find(cumsum(D.pca_npm.sig) >= 1:length(D.pca_npm.sig));
ndim = length(whsig) + 1;

try, saveas(gcf,'Cluster_eigenplot','tif');,end

% --------------------------------------
% predict behavior with component scores
% --------------------------------------

fprintf(1,'\n-----------------------------------------------------------------\n');
fprintf(1,'Stepwise regression: Predictions of behavior with component scores');
fprintf(1,'\n-----------------------------------------------------------------\n');
        
D.PCR = stepwise_tor(score(:,1:ndim),y);

sig = find(D.PCR.inmodel);
D.PCR.sscore = score(:,sig);    % scores for significant components

% plot all significant (now done in cluster_discrim_montage)
%figure('Color','w');
%for i = 1:size(D.PCR.sscore,2)
%    subplot(1,size(D.PCR.sscore,2),i);
%    discrim_plot(D.PCR.sscore(:,i),y,0);
%end

D.PCR.seigv = v(:,sig);         % eigenvectors for sig components
D.PCR.seig = e(sig);            % significant eigenvalues

tmpc = corrcoef([D.PCR.sscore x]);
tmpc(1:size(D.PCR.sscore,2),:) = []; tmpc = tmpc(:,1:size(D.PCR.sscore,2));
D.PCR.compcorr = tmpc;
D.PCR.compcorr_descrip = 'Regions x sig. comps, correlation with component values';
D.PCR.whichcomps = sign(D.PCR.compcorr) .* real(abs(D.PCR.compcorr) > .4);

tmpc(abs(tmpc) < .4) = 0;
D.PCR.threshcomps = tmpc;

% --------------------------------------
% manova to get categorical canonical discrim functions
% this gives us overall p-value and discrim functions
% --------------------------------------
[manovasig,manovap,D.manova] = manova1(score(:,1:ndim),y2);
D.manova.p = manovap;
fprintf(1,'\n-----------------------------------------------------------------\n');
fprintf(1,'Manova to discriminate high from low behavior');
fprintf(1,'\n-----------------------------------------------------------------\n');
fprintf(1,'Wilks lambda (%3.0f,%3.0f) = %3.2f\tChi-sq (%3.0f) = %3.2f\t, p = %3.4f\t\n',D.manova.dfB,D.manova.dfW,D.manova.lambda,D.manova.chisqdf,D.manova.chisq,D.manova.p);
%discrim_plot(D.manova.canon(:,1),y);
fprintf(1,'Eigenvalues\t')
fprintf(1,'%3.2f\t',D.manova.eigenval')
fprintf(1,'\n')
fprintf(1,'Eigenvector 1\t')
fprintf(1,'%3.2f\t',D.manova.eigenvec(:,1)')
fprintf(1,'\n')

% --------------------------------------
% visualize clusters on brain
% --------------------------------------

cluster_discrim_montage(cl,D.PCR.threshcomps,D.PCR.sscore,D.beh);

% --------------------------------------
% manova on individual regions??
% --------------------------------------


diary off


return







