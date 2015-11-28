function robustfit(fmri_model_obj, fmri_data_obj, varargin)
% Robust fit for a model object to data object
%
% :Usage:
% ::
%
%     robustfit(fmri_model_obj, fmri_data_obj, [optional args])
%
% :Features:
%    spatial smoothing of weights at 12 mm FWHM
%    ridge regression ***not yet***
%
% :Preproc scaling:
%   1. Remove covariates using ridge reg; ridge trace for full model
%   2. scale to % signal change across time (cols) OR rank time points (for
%      w/i ss predictions??) AND/OR rank or center rows (images; for 'shape'
%      analysis
%
% :Example:
%    %sig across time, rank across rows: relative % sig change
%
%    Different models of noise lead to different ideas about optimal preproc
%    If large diffs in nuisance scaling in BOLD across individuals, ranking cols may
%    be good idea. but then individual diffs in overall activity will be removed...
%
% :Options:
%
%   **tune:**
%        tuning const for robust reg
%
%   **iter:**
%        'maxiterations', robust reg /WLS iterations. 1 = OLS only!
%
%   **smooth:**
%        'spatial_smooth_fwhm', 0 or smoothing kernel for weights
%
%   **nosmooth:**
%        spatial_smooth_fwhm = 0;
%
%   **stats:**
%        'calculate_stats', calculate_stats = 1; IN DEVELOPMENT
%
%   **noresiduals:**
%        write_residuals = 0;
%
%   **noplots:**
%        save_plots = 0;

% ..
%    Defaults/constants
% ..

tune = 4.6850;  % from Matlab's implementation
bisquare_fcn = @(r) (abs(r)<1) .* (1 - r.^2).^2;

maxiterations = 12;
spatial_smooth_fwhm = 12;   % if 0, do not smooth at all

calculate_stats = 0;
write_residuals = 1;
save_plots = 1;

% Rotate nuisance regressors only to PCA space to avoid colinearity in wfit
fmri_model_obj = rotate_to_pca(fmri_model_obj);

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'tune', tune = varargin{i+1};
            case {'iter', 'maxiterations'}, maxiterations = varargin{i+1};
            case {'smooth', 'spatial_smooth_fwhm'}, spatial_smooth_fwhm = varargin{i+1};
                
            case 'nosmooth', spatial_smooth_fwhm = 0;
            case {'stats', 'calculate_stats'}, calculate_stats = 1;
            case 'noresiduals', write_residuals = 0;
            case 'noplots', save_plots = 0;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

if maxiterations == 1, disp('OLS only: no robust iterations'); end

% display analysis info
% ----------------------------------------------------

SETUP = struct('fmri_data', [], 'fmri_model', [], 'algorithm_info', []);

SETUP.fmri_data = struct('source_notes', fmri_data_obj.source_notes, ...
    'image_names', fmri_data_obj.image_names, ...
    'mask', fmri_data_obj.mask, ...
    'mask_descrip', fmri_data_obj.mask_descrip);

SETUP.fmri_data.history = fmri_data_obj.history;

SETUP.fmri_model = fmri_model_obj;

SETUP.algorithm = struct('name', 'robust regression', ...
    'tune', tune, ...
    'maxiterations', maxiterations, ...
    'spatial_smoothing_fwhm', spatial_smooth_fwhm);

% do not do this...for some reason, bug makes SETUP file 600 MB
% SETUP.algorithm.robust_weight_function = bisquare_fcn;

t1 = tic;
fprintf('Saving SETUP.mat');
save SETUP SETUP
t2 = toc(t1);
fprintf(' %3.0f s\n', t2);


% Check data and dims
% ---------------------------------------------------------------------
t1 = tic;

% Plot and save output
if save_plots
    save_plots_inlinefcn;
end

X = fmri_model_obj.xX.X;
Y = fmri_data_obj.X;

fmri_data_obj.X = []; % save space

fprintf('Checking data');
badvox = any(isnan(Y)) | all(Y == 0);
anybadvox = any(badvox); % use later as well
if anybadvox
    fprintf('Voxels with bad values: %3.0f. (all zeros or some NaNs)\n', sum(badvox));
    Y(:, badvox) = [];
    % we will have to insert these later.
end

% We should not need this, as NaNs are illegal in fmri_data objects.
% fprintf('Removing NaNs\n');
% [wasnan, X, Y] = nanremove(X,Y);

[n, p] = size(X);
[n2, v] = size(Y);
if n ~= n2, error('model and data dims do not match'), end

t2 = toc(t1);
fprintf(' %3.0f s\n', t2);


% Initialize vars
% ---------------------------------------------------------------------

[ols_s, s] = deal(zeros(1, v));

tol_for_convergence = .001 * min([nanstd(X) nanstd(Y)]);  % ...times the smallest beta

% Initial OLS model
% ---------------------------------------------------------------------
t1 = tic;

fprintf('Fitting OLS for %3.0f voxels', v);

pinvX = pinv(X);

b = pinvX * Y;

fprintf('Getting residual error\n');

dfe = n - p;
r = Y - X*b;

for i = 1:v
    ols_s(i) = sqrt((r(:, i)' * r(:, i)) ./ dfe); % same as norm(y-X*b) / sqrt(dfe); % same as sqrt((r' * r) ./ (n - p)
end

% estimate leverages for robust reweighting
H = X * pinvX;  % Hat matrix; diags are leverage
h = diag(H);
adjfactor = 1 ./ sqrt(1-h);
adjfactor = repmat(adjfactor, 1, v);

t2 = toc(t1);
fprintf(' %3.0f s\n', t2);

% -------------------------------------------------------------------
% Iterate robust regression
% -------------------------------------------------------------------

iter = 1;                % Set counter
b0 = b;
max_dev = NaN;

fprintf('Iteration ');

while (iter < maxiterations)  && (iter == 1 || max(sum((b - b0).^2)) > tol_for_convergence)
    
    t1 = tic;
    
    fprintf('%3.0f\n', iter);
    
    %fprintf('\b\b\b\b%3.0f ', iter);
    iter = iter + 1;
    
    % -------------------------------------------------------------------
    % Get weights for robustfit
    % -------------------------------------------------------------------
    
    % adjust residuals
    % -------------------------------------------------------------------
    fprintf('  Weighting residuals');
    
    radj = r .* adjfactor;
    
    for i = 1:v
        s(i) = madsigma(radj(:, i), p);
    end
    s(s == 0) = 1;
    
    radj = radj ./ repmat(s * tune, n, 1);
    w = bisquare_fcn(radj);
    
    
    erase_string('  Weighting residuals');
    
    % Takes too long.  Just do at end.
    % %     % Smooth weights
    % %     % -------------------------------------------------------------------
    % %
    % %     if spatial_smooth_fwhm
    % %
    % %         fprintf('  Smoothing weights in 3D space: %04d', 0);
    % %
    % %         % w is n x v, pass in transposed:
    % %         w = smooth_3d(w', fmri_data_obj.mask.volInfo, spatial_smooth_fwhm, badvox)';
    % %
    % %         erase_string(sprintf('  Smoothing weights in 3D space: %04d', 0));
    % %
    % %     end
    
    
    % -------------------------------------------------------------------
    % Update betas and residuals
    % -------------------------------------------------------------------
    
    fprintf('  Updating betas for all voxels');
    
    %     for i = 1:v
    %         if i > 10, warning off; end  % if weights are 0, this will give a warning  % need to check out ***
    %
    %         b(:, i) = wfit(Y(:, i), X, w(:, i));
    %     end
    
    b = wfit_allvoxels(Y, X, w, size(X, 2), v);
    
    warning on
    
    erase_string('  Updating betas for all voxels');
    
    fprintf('  Getting residuals');
    
    r = Y - X * b;         % residuals
    
    erase_string('  Getting residuals');
    
    % save these
    max_dev(iter) = max(sum((b - b0).^2));
    
    t2 = toc(t1);
    fprintf(' %3.0f s\n', t2);
    
end % Loop

% -------------------------------------------------------------------
% -------------------------------------------------------------------
%
% Final model - calculate and save results
%
% -------------------------------------------------------------------
% -------------------------------------------------------------------

% Smooth weights
% -------------------------------------------------------------------

if maxiterations > 1 && spatial_smooth_fwhm
    
    fprintf('Smoothing weights in 3D space: %04d', 0);
    
    % w is n x v, pass in transposed:
    w = smooth_3d(w', fmri_data_obj.mask.volInfo, spatial_smooth_fwhm, badvox)';
    
    erase_string(sprintf('Smoothing weights in 3D space: %04d', 0));
    
    fprintf('Recalculating final b and r');
    b = wfit_allvoxels(Y, X, w, size(X, 2), v);
    r = Y - X * b;         % residuals
end


SETUP.algorithm.max_dev = max_dev;
SETUP.algorithm.max_dev_descrip = 'Max deviation from old to new betas in interation; check convergence.';

% Write mask, betas, residuals, weights
% ---------------------------------------------------------------------

vI = fmri_data_obj.mask.volInfo;

% write mask
iimg_reconstruct_vols(fmri_data_obj.mask.dat, vI, 'outname', 'mask.img');


% write betas and weights
write_image(b', 'regression_betas.img', vI, badvox);

if maxiterations > 1
    write_image(w', 'robust_reg_weights.img', vI, badvox);
end

if write_residuals
    write_image(r', 'residuals.img', vI, badvox);
end

% Stats, if asked for
% ---------------------------------------------------------------------

fprintf('Final stats for all voxels\n');

if calculate_stats
    
    % we already have r and b
    
    % Standard error
    MSE = single(zeros(1,v));
    v = single(zeros(n, v));
    
    for i = 1:v
        W = diag(w(:, i));                    % Weight matrix, W = V^-1, inverse of cov.matrix
        
        MSE(i) = r(:,i)' * W * r(:,i);
        
        v(:, i) = diag(inv(X' * W * X)) .* MSE(i);   % variances
        
    end
    
    % could smooth MSE and dfe
    % could replace dfe with robust one??
    
    MSE = MSE ./ dfe;
    t = b ./ sqrt(v);
    
    p = 2 * ( 1 - tcdf(abs(t),dfe) );
    
    write_image(t, 'regressors_t.img', vI, badvox)
    write_image(p, 'regressors_p.img', vI, badvox)
    
end

% Write fitted responses (HRFs)
% ---------------------------------------------------------------------

% For each condition and each image...
% Save HRFs(condition-specififc fits)
% Estimate htw
fprintf('Writing 4-D hrf images and height, time to peak, width, AUC images\n');

% list all bf
wh_bf = [1 2 1 2];

c = fmri_model_obj.xX.cond_assignments;

for i = 1:size(c, 2)
    t1 = tic;
    
    fprintf('   Condition %3.0f, %s\n', i, fmri_model_obj.xX.name{i});
    
    cbeta = b(c(:, i), :);
    cfit = fmri_model_obj.xBF(wh_bf(i)).bf * cbeta;
    
    % estimate h, t, w
    % ---------------------------------------------------------------------
    
    fprintf('Estimating h, t, w, AUC');
    
    % peak in first 20 sec of response, or length of response if shorter
    % period is modeled
    hconstraint = min( round(20 ./ fmri_model_obj.xBF(wh_bf(i)).dt), size(fmri_model_obj.xBF(wh_bf(i)), 1) );
    [h, t, wid, auc] = deal(zeros(v, 1));
    
    for j = 1:v
        [h(j, 1), t(j, 1), wid(j, 1), w_times, halfh, auc(j, 1)] = fir2htw2(cfit(:, j), hconstraint, 0);
    end
    
    erase_string('Estimating h, t, w, AUC');
    
    t2 = toc(t1);
    fprintf('      h, t, w, AUC estimated in %3.0f s\n', t2);
    
    % write images
    % ---------------------------------------------------------------------
    t1 = tic;
    
    write_image(h, sprintf('amplitude_cond%03d_%s.img', i, fmri_model_obj.xX.name{i}), vI, badvox);
    write_image(t, sprintf('time_to_peak_cond%03d_%s.img', i, fmri_model_obj.xX.name{i}), vI, badvox);
    write_image(wid, sprintf('duration_cond%03d_%s.img', i, fmri_model_obj.xX.name{i}), vI, badvox);
    write_image(auc, sprintf('AUC_cond%03d_%s.img', i, fmri_model_obj.xX.name{i}), vI, badvox);
    
    % downsample
    cfit_at_TR = cfit(1:16:end, :);
    
    %write hrf images
    write_image(cfit_at_TR', sprintf('hrf_cond%03d_%s.img', i, fmri_model_obj.xX.name{i}), vI, badvox);
    
    t2 = toc(t1);
    fprintf('      Images written in %3.0f s\n', t2);
end

% Ridge reg

% Satterthwaite stuff
% See Module8_GLS_timeseries

% % % H = X * inv(X' * W * X) * X' * W;
% % % trace(H)
% % %
% % %
% % % % This could be matrix-ized
% % %
% % % sw = sqrt(w(:, i));
% % % [r c] = size(X);
% % % yw = Y(:, i) .* sw;
% % % xw = X .* sw(:,ones(1,c));
% % % [Q,R]=qr(xw,0);
% % % b = R\(Q'*yw);
% % %
%
% v = repmat(diag(invxwx),1,v) .* repmat(MSE,k,1);
%
% t = b ./ sqrt(v);
%
% if nargout > 1, p = 2 * ( 1 - tcdf(abs(t),dfe) ); end
%
% Wi{i} = diag(W(:, i));              %
%
%             invxvx{i} = inv(X' * Wi{i} * X);       % Save these for later, for speed; don't need to re-use
%             bforming{i} = invxvx{i} * X' * Wi{i};
%  % R = Wi.^.5 * (eye(n) - X * invxvx * X' * Wi{i});          % Residual inducing matrix
%             R = Wi{i} .^ .5 * (eye(n) - X * bforming{i});               % Residual inducing matrix
%
%             Q = R * inv(Wi{i});                                         % Q = RV
%             dfe(i) = (trace(Q).^2)./trace(Q * Q);                       % Satterthwaite approximation for degrees of freedom


% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
%
% Inline functions
%
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------


    function save_plots_inlinefcn
        
        savedir = fullfile(pwd, 'output_images');
        
        plot(fmri_model_obj);
        saveplots(fmri_model_obj, savedir);
        
        plot(fmri_data_obj)
        saveplots(fmri_model_obj, savedir);
        
    end



end % function



function [b,R,xw] = wfit(y,x,w)

% weighted least squares fit
% one voxel

sw = sqrt(w);
[r c] = size(x);
yw = y .* sw;
xw = x .* sw(:,ones(1,c));
[Q,R]=qr(xw,0);
b = R\(Q'*yw);

end



function [b,R,xw] = wfit_allvoxels(y, x, w, k, v)

% ***why are some weights > 1???

% weighted least squares fit
% all voxels

sw = w .^ .5;

yw = y .* sw;

b = zeros(k, v, 'single');

for i = 1:v
    xw = x .* sw(:,i * ones(1,k));
    
    [Q,R]=qr(xw, 0);
    
    b(:, i) = R\(Q'*yw(:, i));
end

end


function s = madsigma(r,p)

% Compute sigma estimate using MAD of residuals

m = median(r);
rs = sort(abs(r-m));
if (abs(m) > rs(end))
    % Unexpectedly all residuals are very small
    rs = sort(abs(r));
end
s = median(rs(p:end)) / 0.6745;
if (s==0), s = .5*mean(rs); end

end

%
% function [t,p,b,v] = weighted_glmfit(X,Y,wts,dfe)
%
% v = size(Y,2);
%
%
% W = diag(wts);                    % Weight matrix
% invxwx = inv(X'*W*X);
% bform = invxwx * X'* W;         % beta-forming matrix.  hat = X * bform
% %
% % % rows are columns of design (X), cols are Y variables
% b = bform*Y;
%
% k = size(b,1);
%
% r = Y - X * b;         % residuals
%
% % % standard error
% MSE = zeros(1,v);
% for i=1:v
%     W = diag(w(:, i));                    % Weight matrix
%     MSE(i) = r(:,i)'*W*r(:,i);
% end, MSE = MSE/dfe;
%
% v = repmat(diag(invxwx),1,v) .* repmat(MSE,k,1);
%
% t = b ./ sqrt(v);
%
% if nargout > 1, p = 2 * ( 1 - tcdf(abs(t),dfe) ); end
%
% end
%
%
% function erase_string(str1)
% fprintf(1,repmat('\b',1,length(str1))); % erase string
% end



function yout = zeroinsert(wasbad, y)
% yout = zeroinsert(wasbad, y)
%
% Re-insert removed CASES (rows) and fill with zeros
% wasbad is indicator for removed cases, of size size(yout).  y is data.
%
% if you enter y', inserts VOXELS (cols).  here, pass in v x n matrix, y'
% to fill empty/removed voxels

% See nanremove.m and naninsert.m

wh = find(wasbad);

if isempty(wh), yout = y; return, end

yout = zeros(length(wasbad), size(y, 2));

yout(1:wh(1) - 1, :) = y(1:wh(1) - 1, :);

for i = 2:length(wh)
    ystart = wh(i-1) + 2 - i; % NaN index value - num previous removed; wh(i-1) + 1 - i + 1;
    yend = wh(i) - i; % wh(i) - 1 - i + 1
    
    yout(wh(i-1) + 1 : wh(i) - 1, :) = y(ystart : yend, :);
end
% last segment
i = i+1;
if isempty(i), i = 2; end
ystart = wh(end) + 2 - i;

if ystart <= size(y, 1)
    yout(wh(end)+1:end, :) = y(ystart:size(y, 2), :);
end

end



function write_image(x, myname, vI, badvox)

t1 = tic;

% Replace bad vox for image writing

if any(badvox)
    x = zeroinsert(badvox, x);
end

mypath = fullfile(pwd, myname);

imgv = image_vector('X', x, 'volInfo', vI, 'filename', myname, 'fullpath', mypath);
write(imgv);

t2 = toc(t1);
fprintf('Done in %3.0f s\n', t2);

end



function w = smooth_3d(w, volInfo, sfwhm, varargin)
% Smooth 3-D images stored in columns with FWHM in mm of sfwhm
%
% function w = smooth_3d(w, volInfo, sfwhm, [badvox])
%
% Take v x n matrix w, and smooth columns (images), returning v x n matrix
% again
% Optional: if v is smaller than the original image because some voxels
% were removed, enter logical vector badvox, and the missing voxels
% will be filled with zeros.


% NOTE: horribly slow and memory intensive:
%wvol = iimg_reconstruct_vols(w', fmri_data_obj.mask.volInfo);

if ~isempty(varargin), badvox = varargin{1}; end

% transpose, and ...
% insert zeros back into bad vox
if any(badvox)
    w = zeroinsert(badvox, w);
end

n = size(w, 2);

try
    
    for i = 1:n
        fprintf('\b\b\b\b%04d', i);
        
        wvol = iimg_reconstruct_vols(w(:, i), volInfo);
        
        spm_smooth(wvol, wvol, [sfwhm sfwhm sfwhm]);
        
        % go back to vector
        wvec = wvol(volInfo.wh_inmask);
        
        w(:, i) = wvec;
        
    end
    
catch
    disp('Error with volume reconstruction. Mask in fmri_data obj not properly defined?');
    disp('Stopped in debugger so you can check...');
    keyboard
end

if any(badvox)
    
    w(badvox, :) = [];
    
end

end % function



