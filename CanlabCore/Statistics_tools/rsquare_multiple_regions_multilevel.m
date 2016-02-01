function stats = rsquare_multiple_regions_multilevel(Y, X, varargin)
% Predict outcome var (y) using data from multiple vars (X, e.g., brain regions)
% Test R-square against permuted data
%
% :Usage:
% ::
%
%     stats = rsquare_multiple_regions_multilevel(Y, X, varargin)
%
% :Inputs: Variable args
%
%   **'colors':**
%        followed by colors cell ({'r' 'g' 'b'}) for each col. of X
%
%   **'nperms:**
%        then num perms
%
% :Examples:
% ::
%
%    % SETUP data
%    cl = [clpos_data clneg_data];
%
%    cl(cat(1, cl.numVox) < min_cluster_size) = [];
% 
%    % get brain data cell
%    for i = 1:size(cl(1).all_data, 2) 
%        for c = 1:length(cl)
%            data{i}(:,c) = cl(c).all_data(:, i); 
%        end
%    end
%
%    % NMDS ANALYSIS
%    OUT = [];
%
%    OUT.ridge = matrix_direct_effects_ridge(data);
%    D = OUT.ridge.mean; D(find(eye(size(D)))) = 1;
%    D = (D' + D) ./ 2;
%    OUT.ridge.D = (1 - D) ./ 2;
%    [OUT.stats_mds.GroupSpace,OUT.stats_mds.obs,OUT.stats_mds.implied_dissim] = shepardplot(OUT.ridge.D,[]);
%
%    OUT.stats_mds = nmdsfig_tools('cluster_solution',OUT.stats_mds, OUT.stats_mds.GroupSpace, 2:max_networks, nperms, []);
%    OUT.stats_mds.colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};
%    create_figure('nmdsfig');
%
%    OUT.stats_mds.names = [];
%    nmdsfig(OUT.stats_mds.GroupSpace,'classes',OUT.stats_mds.ClusterSolution.classes,'names',OUT.stats_mds.names,'sig',OUT.ridge.fdrsig);
%    hh = nmdsfig_fill(OUT.stats_mds);
%    axis image; axis equal
%
%
%    % Multiple regions predict behavior
%
%    % Design matrix with cluster averages
%    classes = OUT.stats_mds.ClusterSolution.classes;
%    clear X
%    for i = 1:length(data)
%        for j = 1:max(classes)
%            X{i}(:, j) = nanmean(data{i}(:, classes == j), 2); 
%        end
%    end
%
%    OUT.stats_regression = rsquare_multiple_regions_multilevel(acc_x_dist, X, 'colors', OUT.stats_mds.colors, 'nperms', 100);
%
% ..
%    Tor Wager, June 2008
%
%    var args not all done yet (in development)
% ..


nperms = 100;


for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                
                    % functional commands
                case 'colors', colors = varargin{i+1};

                case 'nperms', nperms = varargin{i + 1};
                    
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
end
    
%% initial model fits, all together, then separately to get order
stats = glmfit_multilevel(Y, X, ones(length(X), 1), 'weighted');

nvars = size(X{1}, 2);

for i = 1:nvars
    for s = 1:length(X)
        Xr{s} = [X{s}(:, i)]; % reduced model, with only this regressor.
    end

    stats_ind = glmfit_multilevel(Y, Xr, ones(length(X), 1), 'weighted');
    stats.separate_reg_pvals(i) = stats_ind.p(2);

end

%%% **** Notes: should do best 1, then best combo of 2, best 3...etc ****
%%% permute additional param after best combo from prev. step

%% Predicting using variable # regions
[tmp, myorder] = sort(stats.separate_reg_pvals); % eliminate intercept
stats.myorder = myorder;

for i = 1:length(myorder)

    for s = 1:length(X)

        Xr{s} = [ones(size(X{s}, 1), 1) X{s}(:, myorder(1:i))]; % reduced model, with regressors ordered by predictive accuracy.

        [B,BINT,R,RINT,sstats] = regress(Y{s}, Xr{s});

        r2(s, i) = sstats(1);

    end

end

create_figure('R-squared'); barplot_columns(r2, 'R-squared', [], 'nofig', 'noind');
xlabel('Number of networks in model');
ylabel('Average R^2 value');

if exist('colors', 'var') && ~isempty(colors)

    for i = 1:size(r2, 2)
        hbar(i) = bar(i, nanmean(r2(:, i)), 'FaceColor', colors{myorder(i)}(1));
    end

end

stats.observed_r2 = r2;

drawnow

%% permutation

clear permindx
for i = 1:length(X)
    permindx{i} = permute_setupperms(size(X{i}, 1), nperms);
    
    %Xrp{s} = [ones(size(X{s}, 1), 1) X{s}(:, myorder(1:i))]; % reduced model, with regressors ordered by predictive accuracy.
    
end

[r2meanp, b1mean] = deal(zeros(nperms, length(myorder)));

    
%%
for p = 1:nperms
    
    fprintf('%3.0f ', p);

    if mod(p, 20) == 0, fprintf(1, '\n'); end

    [r2, bb] = deal(zeros(length(X), length(myorder)));
    
    for i = 1:length(myorder)

        for s = 1:length(X)

            % Method 1: Permute all variables (permute data)
            Xr{s} = [ones(size(X{s}, 1), 1) X{s}(:, myorder(1:i))]; % reduced model, with regressors ordered by predictive accuracy.

            Yp{s} = Y{s}(permindx{s}(p, :)', :); % permute

% %             % Method 2: Permute only the i:endth regressors
% %             % part to not permute: conditional on this being in model
% %             Xsnp = X{s}(:, myorder(1:i-1));
% %             
% %             % part to permute: if everything else were random...
% %             Xsp = X{s}(:, myorder(i:end));
% %             Xsp = Xsp(permindx{s}(p, :)', :); % permute
% %             Xr{s} = [ones(size(X{s}, 1), 1) Xsnp Xsp]; % reduced model, with regressors ordered by predictive accuracy.
% %             Yp{s} = Y{s};

            [B,BINT,R,RINT,sstats] = regress(Yp{s}, Xr{s});

            r2(s, i) = sstats(1);

            bb(s, i) = B(2);
            
            % %         Xrp_ml{s} = [X{s}(permindx{s}(p, :)', myorder(1:i))]; % reduced model, with regressors ordered by predictive accuracy.

        end

        % %     stats = glmfit_multilevel(Y, Xrp_ml, ones(length(X), 1), 'weighted', 'noverbose');
        % %     multilev_pvals{i}(p, :) = stats.p;

    end

    r2meanp(p, :) = nanmean(r2);
    
    b1mean(p, :) = nanmean(bb);
    
end
    
fprintf(1, '\n');

bar((1:nvars) + .2, mean(r2meanp), 'FaceColor', [.2 .2 .2]);

%% Z-scores and p-values based on perm dist
fp = r2meanp > repmat(nanmean(stats.observed_r2, 1), nperms, 1);

stats.permute.pvals = sum(double(fp), 1) ./ nperms;

yval = nanmean(stats.observed_r2) + .15 * nanmean(stats.observed_r2);

wh = find(stats.permute.pvals <= .001);
if ~isempty(wh), for i = 1:length(wh), text(wh(i), yval(wh(i)), '***', 'FontSize', 18); end, end

wh = find(stats.permute.pvals <= .01 & stats.permute.pvals > .001);
if ~isempty(wh), for i = 1:length(wh), text(wh(i), yval(wh(i)), '**', 'FontSize', 18); end, end

wh = find(stats.permute.pvals <= .05 & stats.permute.pvals > .01);
if ~isempty(wh), for i = 1:length(wh), text(wh(i), yval(wh(i)), '*', 'FontSize', 18); end, end

%%
disp('Added .permute field to stats output strucuture.');
stats.permute.nperms = nperms;
stats.permute.r2meanp = r2meanp;
stats.permute.r2mean_ub = prctile(stats.permute.r2meanp, 95);
stats.permute.permindx = permindx;

stats.permute.b1mean = b1mean;
