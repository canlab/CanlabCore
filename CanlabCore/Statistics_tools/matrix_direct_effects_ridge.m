function out = matrix_direct_effects_ridge(data, varargin)
% :Usage:
% ::
%
%     [out] = matrix_direct_effects_ridge(data, [optional args]);
%
% Take an n x v matrix of data for each of N subjects.
% data are entered in a square matrix within each cell of data,
% data{1}, {2}, etc.
%
% Uses ridge regression to assess linear slopes for each variable on
% each other one.
%
% :Optional Inputs:
%
%   **'k':**
%        followed by ridge parameter value.  Default = 0.
%
%   **'scale' or 'zscore':**
%        which z-scores within-subjects
%
% Future: Consider weighting based on multicolinearity
%
% Note: Matrix is asymmetrical! Rows predicting columns...
%
% :Examples:
%
% Take mediation brain results clusters cell clpos_data{...}
% Concatenate data matrix and run to get connections.
% ::
%
%    for i = 1:length(clpos_data), data{i} = cat(2, clpos_data{i}.timeseries); end
%    out = matrix_direct_effects_ridge(data);
%    [out.GroupSpace,out.obs_dist,out.implied_dissim] = shepardplot(out.dissim,[]);
%    create_figure('nmdsfig');
%    nmdsfig(out.GroupSpace,'classes',ones(out.k, 1),'names',[],'sig',out.fdrsig, 'legend',{'Pos' 'Neg'},'sizescale',[4 16]);
%
% Example 2:
% Include mediation predictor and outcome
% ::
%
%    for i = 1:length(clpos_data), data{i} = [SETUP.data.X{i} SETUP.data.Y{i} data{i}]; end
%    out = matrix_direct_effects_ridge(data);
%
% Example of MDS and plotting direct relationships:
% ::
%
%    OUT.ridge = matrix_direct_effects_ridge(data);
%    D = OUT.ridge.mean; D(find(eye(size(D)))) = 1;
%    D = (D' + D) ./ 2;
%    OUT.ridge.D = (1 - D) ./ 2;
%    [OUT.stats_mds.GroupSpace,OUT.stats_mds.obs,OUT.stats_mds.implied_dissim] = shepardplot(OUT.ridge.D,[]);
%    OUT.stats_mds = nmdsfig_tools('cluster_solution',OUT.stats_mds, OUT.stats_mds.GroupSpace, 2:10, 1000, []);
%    OUT.stats_mds.colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};
%    create_figure('nmdsfig');
%    nmdsfig(OUT.stats_mds.GroupSpace,'classes',OUT.stats_mds.ClusterSolution.classes,'names',OUT.stats_mds.names,'sig',OUT.ridge.fdrsig);
%    hh = nmdsfig_fill(OUT.stats_mds);
%
% Prediction using multiple regions:
% ::
%
%    stats = rsquare_multiple_regions_multilevel(Y, X, varargin)
%
% Visualization of networks:
% ::
%
%    classes = OUT.stats_mds.ClusterSolution.classes;
%    % create_figure('Surfaces');
%    % shan = addbrain;
%    shan = mediation_brain_surface_figs({}, {});
%    scolors = {[1 0 0] [0 1 0] [0 0 1] [1 1 0] [0 1 1] [1 0 1] [0 0 0] [1 0 0] [0 1 0] [0 0 1] [1 1 0] [0 1 1] [1 0 1] [0 0 0] };
%    for j = 1:max(classes)
%    wh = (classes == j); sum(wh)
%    montage_clusters([], cl(wh), scolors(j))
%    % create_figure('Surfaces', 1, 1, 1);
%    cluster_surf(cl(wh), 5, scolors(j), shan);
%    end
%
% :See Also: xcorr_multisubject.m
%
% ..
%    Tor Wager, March 2008
% ..

kval = 0; % 0 is no ridge reg; other k > 0 implements ridge
doscale = 0;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case 'k', kval = varargin{i+1};
            case {'zscore', 'scale'}, doscale = 1;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if ~iscell(data)
    % single-level
    
    error('Not implemented yet.');
    
else
    % multi-level
    N = length(data);
    k = size(data{1}, 2); % number of variables
    
    tthr = tinv(1 - (.05) ./ 2, N - 1);
    
    betamatrix = zeros(k, k, N);
    
    for i = 1:N % For each subject
        
        fprintf(' %03d', i);
        
        b = build_matrix(data{i}, kval, doscale);
        
        betamatrix(:,:,i) = b;
        
    end
    
    fprintf(' Done.\n');
end

%b = nanmean(betamatrix, 3);
[b, t, sig, out] = ttest3d(betamatrix);

out.dissim = (1 - out.mean) ./ 2;
out.dissim(find(eye(size(out.dissim)))) = 0; % make diagonal 0
out.dissim = (out.dissim' + out.dissim) ./ 2; % make symmetrical

out.ridge_param = kval;
out.N = N;
out.k = k;

end



function b = build_matrix(sdata, kval, doscale)

k = size(sdata, 2); % number of variables
b = zeros(k);

for v = 1:k  % For each variable
    
    y = sdata(:, v);
    X = sdata;
    X(:, v) = [];
    
    [nanvec, X] = nanremove(X'); X = X';
    
    % standardize
    if doscale
        y = scale(y);
        X = scale(X);
    end
    
    % intercept
    if isempty(X) % we have only one predictor
        X = ones(size(y)); % just intercept
    else
        X(:,end + 1) = 1;
    end
    
    %fprintf(' %03d %3.2f\n', i, cond(X));
    
    if rank(X) < size(X, 2)
        % rank deficient; skip
        b(:, v) = NaN;
    else
        
        bsubi =  ridge(y, X, kval);
        bsubi = bsubi(1:end-1);              % exclude intercept
        bsubi = naninsert(nanvec, bsubi);
        
        % insert space for v in matrix of betas
        indx = zeros(size(nanvec, 1) + 1, 1);
        indx(v) = 1;
        bsubi = naninsert(indx, bsubi);
        
        b(:, v) = bsubi;
        
    end
    
end  % variable within matrix

end




% b(:, end) = []; % remove intercept parameters
%
% bpop = nanmean(b);
% t = bpop ./ ste(b);
%
% wh = find(t > tthr);
% disp('Positively related to X:');
% disp(wh);
% whn = find(t < -tthr);
% disp('Negatively related to X:');
% disp(whn);
%
% allt(indx, :) = t;
%
% end
%
% create_figure('t-vals by ridge param'); plot(kvals, mean(abs(allt), 2))
%
%
%







% %
% %
% %
% %
% %
% %
% %
% % [rows,cols] = find(tril(sig));
% %
% % n = size(sig,1);
% % indx = 1:n;
% %
% % n_sig_pairs = length(rows);
% %
% % % outputs
% % direct = ones(1, n_sig_pairs);
% % mediated = zeros(1, n_sig_pairs);
% %
% % mediators = zeros(1, n);
% %
% % t1 = clock;
% % fprintf(1, '\nTesting for mediators of %d significant bivariate relationships: %6d', n_sig_pairs, 1);
% %
% % for i = 1:n_sig_pairs
% %
% %     fprintf(1, '\b\b\b\b\b\b%6d', i);
% %
% %     this_pair = [rows(i) cols(i)];
% %
% %     % x and y, the relationship potentially mediated by another
% %     x = data(:, rows(i));
% %     y = data(:, cols(i));
% %
% %     % potential mediators
% %     indx_other_regions = indx(~(indx == rows(i) | indx == cols(i)));
% %     m_candidates = data(:, indx_other_regions);
% %
% %     % search for mediators; see mediation.m
% %     med_results = mediation_search('M', x, y, m_candidates, 'noverbose');
% %
% %     sigfx = med_results.pvals < .05;
% %     sig_mediation = sigfx(:,1) & sigfx(:,2) & sigfx(:,5);
% %
% %     % preserve signs for later use; store as 1 or -1 in matrix
% %     mediation_sign = sign(med_results.paths(:,5));
% %     direct_sign = sign(med_results.paths(1,4));     % should all be same!
% %
% %     complete_mediation = sig_mediation & ~sigfx(:,3);
% %
% %     % collect final output
% %     if any(complete_mediation)
% %         direct(i) = 0;
% %
% %         mediators(complete_mediation) = mediation_sign(complete_mediation);
% %     else
% %         direct(i) = direct_sign(1);
% %     end
% %
% %     if any(sig_mediation)
% %         mediated(i) = 1;
% %     end
% %
% %
% % end
% %
% % fprintf(1, '\nDone in %d s\n', etime(clock,t1))
% %
% % % collect output matrices
% % mediated_mtx = zeros(n);
% % direct_mtx = zeros(n);
% %
% % for i = 1:n_sig_pairs
% %     mediated_mtx(rows(i), cols(i)) = mediated(i);
% %     direct_mtx(rows(i), cols(i)) = direct(i);
% %
% %     mediated_mtx(cols(i), rows(i)) = mediated(i);
% %     direct_mtx(cols(i), rows(i)) = direct(i);
% % end
% %
% %
% % end

