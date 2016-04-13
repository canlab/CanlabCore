function [X,mP,spmP] = tor_get_physio(varargin)
% :Usage:
% ::
%
%     [X,mP,spmP] = tor_get_physio([mP],[spmP],[nvoxels],[doortho])
% arguments are optional, but you must enter them in this order.
%
% Tor Wager 10/21/02
%
% Get nuisance covariates likely to be related to physiological noise and head motion
% The algorithm:
%
% The program extracts raw/preprocessed image data from the ventricles (CSF space), as
% defined by a mask denoting which voxels are CSF for that subject.  
% Either all voxels or a randomly selected subset [nvoxels] is subjected to
% principal components analysis, to determine regular patters of drift over time
% and across voxels.  Those patterns are expected to be related to global signal drift,
% head movement, and physiological noise, and are assumed to be UNrelated to the task
% of interest, by virtue of the fact that they occur in the ventricles.
% 
% PCA is done twice on the timeseries' of CSF voxels.  The first time, PCA is done
% on the sums of squared values (not the correlations) of voxel timeserieses across
% the entire experiment, mean-centered based on the whole experiment.  Most of the
% coherent variation in this case is expected to be due to head movement and changes
% in shims/gradients/etc. from run to run.  The SS values are used because we want to
% weight the voxels with the highest variation most heavily, as they are presumably
% picking up most of this signal.  The first 3 eigenvariates (canonical timeseries)
% are saved.
%
% Following, a separate, second PCA is done on the correlation matrix of data
% within each session.  Session data for each voxel are mean-centered and scaled
% relative to the session (variance of each voxel = 1).  We do this because 
% physiological noise-related signals may produce periodic signals of different
% magnitudes in different voxels, and we want to extract the most coherent signals
% we can within each session.  So these eigenvariates are expected to reflect
% primarily noise related to physiology (heart rate, respiration).  Up to 5 eigenvariates
% for each session are saved (nothing with eigenvalue < 1 is saved).
%
% Next, the CSF-related nuisance covariates (eigenvariates from PCA) are combined
% with existing nuisance covariates and intercept columns from the existing 
% design matrix (SPMcfg xX).  The proportion of variance in each predictor of interest
% explained by this nuisance basis set is calculated using regression, and the
% nuisance covariates are orthogonalized with respect to each predictor of interest.
% There are good and bad results of this step.  The bad is that any signal that 
% tracks the predictors is attributed to the task, not to noise, even if it's actually
% caused by physiological artifact.  So the orthogonalized basis set does not
% protect you from physiology or movement-related false positives.  However,
% the nuisance covariates are also unlikely to reduce power in estimating you effects
% of interest.  More importantly, it avoids false positives created when one 
% predictor (A) is more highly correlated with the nuisance covariates than another
% (B).  In practice, betas for A will tend to be smaller than B, given the same
% actual response to both, and a random effects analysis on A-B will produce
% false positive activations.  Orthogonalization of the nuisance set precludes this.
%
% :Inputs:
%
%   **mP:**
%        CSF mask image file.  *_seg3.img output from SPM is appropriate
%        should be in same space and have same dims as functionals
%        but automatic reslicing is done if necessary.
%
%   **spmP:**
%        name (full path name preferred) of SPMcfg.mat file to use
%        This contains the design matrix and raw/preproc image file names to use.
% 
%   **nvoxels:**
%        Number of CSF voxels to use in PCA analysis
%        More than 100 can be very slow and memory intensive.
%        Fewer than 100 voxels loads a different way, and may be slower.
%        Best is probably between 100 - 1000.  800 runs pretty fast.
%
%   **doortho:**
%        Orthogonalize nuisance covariates with respect to regs of interest
%        This assumes that any signal that covaries with the task is, in fact,
%        due to the task, so it gives you some bias towards finding positive results.
%        However, the alternative is that nuisance covariates may soak up variance
%        related to the task, and you'll miss activations.
%        In addition, if some regressors are more colinear with the nuisance set,
%        you can create false "activations" when comparing these regressors to other
%        ones.  This problem exists whether or not we choose to model nuisance 
%        covariates.  One solution is to use the ortho when doing random effects analyses,
%        as the sign and magnitude of nuisance-related activations would not be expected to be
%        the same across subjects unless the variance was really task-related.
%        Default is 1, or "yes, do orthogonalization."
%
% for functions called, see this .m file.
%
% :Examples:
% ::
%
%    % get filenames for SPMcfg files and CSF mask for each subject
%    cd C:\Tor_Documents\CurrentExperiments\intext2\RESULTS\model1
%    spmP = get_filename('sub*','SPMcfg.mat');
%    cd C:\Tor_Documents\CurrentExperiments\intext2\
%    mP = get_filename('sub*','anatomy/nscalped_f*seg3.img');
%    % Now run:
%    for i = 1:size(mP,1)  
%        tor_get_physio(mP(i,:),spmP(i,:),300);  % 300 voxels
%        pause(10); close all
%    end

% :Functions called:
%   - spm functions: spm_get, etc.
%   - timeseries2.m	(for < 100 voxels)
%   - read_hdr.m	(big-little endian dependent; validate for your data)
%   - timeseries3.m	(for > 100 voxels; uses SPM's image reading)
%   - reslice_imgs.m
%   - mask2voxel.m 	(only if ind2sub.m from Matlab is not found)


CSFprob = .95;      % this is the value a voxel in the mask img must have to be considered
                    % an 'on' value.
                    % if using SPM segmentation output (e.g., *_seg3) for the mask,
                    % this is something like the prob. of being in CSF
                    % low thresholds will result in HUGE eigenvalue problems
                    % and memory difficulties.

mypwd = pwd;
t1 = clock;

% ----------------------------------------------------------------------------------
% * set up input arguments
% ----------------------------------------------------------------------------------

if length(varargin) > 0
    mP = varargin{1};
else
    % mask filename
    mP = spm_get(1,'*img','Select CSF mask for this subject.');
end
if isempty(mP), mP = spm_get(1,'*img','Select CSF mask for this subject.');, end

if length(varargin) > 1    
    spmP = varargin{2};
else
    spmP = spm_get(Inf,'SPMcfg.mat','Choose SPMcfg.mat file for this subject or Done to skip.');
end

if length(varargin) > 2, srand = varargin{3};, else, srand = 0;, end

if length(varargin) > 3, doortho = varargin{4};, else, doortho = 1;, end

% ----------------------------------------------------------------------------------
% * cd to SPM results directory so we can write output file there
% ----------------------------------------------------------------------------------
d = fileparts(spmP);
eval(['cd ' d])

diary physio_nuisance_covariates.out

% ----------------------------------------------------------------------------------
% * load SPMcfg.mat file, which contains all relevant info except mask image
% ----------------------------------------------------------------------------------

        load(spmP)
        nsess = length(Sess);
        P = str2mat(VY.fname);
	if ~(exist(deblank(P(1,:))) == 2)
		disp(['Looking for: ' P(1,:)])
		disp(['Can''t find original img files!! Please specify.'])
		P = spm_get(Inf,'*.img','Select raw image files.');
		VY = spm_vol(P);
		disp(['VY has been modified!!! Using: ' P(1,:)])
	end
        for i = 1:length(Sess), nimgs(i) = size(Sess{i}.row,2);,end
        
%  [nsess] = spm_input_ui('Enter number of runs, or 0 to choose SPM.mat file [recommended]',.1,'i',[],1);


% ----------------------------------------------------------------------------------
% * reslice the CSF mask if necessary, to be in space of functionals!
% ----------------------------------------------------------------------------------

mV = spm_vol(mP);
pV = spm_vol(P(1,:));
if any(mV.dim(1:3) - pV.dim(1:3)) | any(any(pV.mat(1:3,1:3) - mV.mat(1:3,1:3)))
    [d f e] = fileparts(mP); reslice_imgs(P(1,:),mP,0);
    mP = fullfile(d,['r' f e]);
    disp(['Resliced mask to space of functionals: ' mP])
    mV = spm_vol(mP);
end

% ----------------------------------------------------------------------------------
% * load and check mask
% ----------------------------------------------------------------------------------

mv = spm_read_vols(mV); mv = mv > CSFprob;
fprintf(1,'\nCSF mask has %3.0f voxels out of %3.0f total.\t',sum(mv(:)),prod(size(mv)))
try
    % if we have the right toolbox...
    [x y z] = ind2sub(size(mv),find(mv)); [XYZ] = [x y z];
catch
    disp('elmat toolbox function ind2sub not found; using SLOWER version mask2voxel.m')
    XYZ = mask2voxel(mv);
end

tmp = sum(sum(mv)); 
if ~findobj('Tag','Graphics'), spm fmri; figure(findobj('Tag','Graphics'));,end
if ~findobj('Tag','Interactive'), spm fmri;,end
%imagesc(mv(:,:,find(tmp == max(tmp)))); colormap gray; warning off; title([mP]), warning on,;drawnow
spm_check_registration(str2mat(mP,P(1,:)));

if sum(tmp) > 100, 
    disp(['More than 100 CSF voxels in mask - this may be computationally intensive!'])
    % if srand is not entered as an input argument, prompt.
    if ~srand, 
        figure(findobj('Tag','Interactive'))
        srand = spm_input_ui('Enter n vox to use, or 0 for all',.1,'i',[],1);, 
    end
    if srand > sum(tmp), 
	disp(['More voxels requested than available at threshold ' num2str(CSFprob) ': using ' num2str(max(tmp))])
	srand = sum(tmp);
    end

else
    srand = 0;
end

if srand
    wh = rand(size(XYZ,1),1) * size(XYZ,1);
    wh2 = sort(wh); wh2 = wh2(1:srand);
    for i = 1:srand, XYZ2(i,:) = XYZ(find(wh == wh2(i)),:);, end
    XYZ = XYZ2;
end

% ----------------------------------------------------------------------------------
% * Load images and mask with CSF
% ----------------------------------------------------------------------------------
% M1 is unscaled timeseries over all sessions; eigs computed over whole experiment
% M2 is timeseries scaled within each session; (cell array)
%   eigs are computed within session

fprintf(1,'\nEigenvariates based on %3.0f voxels \t',size(XYZ,1))
fprintf(1,'\nLoading volumes, extracting voxels, and scaling \t')

wh = [0 cumsum(nimgs)];
M1 = [];
M2 = [];

if size(XYZ,1) < 101
    % this is super slow for large n!
    ts = timeseries2('multi',P,struct('coords',XYZ));
    M1 = ts.indiv;
    
    % adjust M2 to mean 0 var 1
    for j = 1:size(M1,2)
        for i = 1:nsess
            ind = wh(i)+1:wh(i+1);
            M2{i}(:,j) = (M1(ind,j) - mean(M1(ind,j))) ./ std(M1(ind,j));
        end
    end
    
    
    
else
    


for i = 1:nsess
    subP = P(wh(i)+1:wh(i+1),:);
    fprintf(1,'.')
    ts = timeseries3(XYZ,subP);
    
    % save 2 matrices, one for session-mean centered and scaled and one uncentered (until overall mean is known)
    M1 = [M1;ts.all_data];
    
    for j = 1:size(ts.all_data,2)
        % center and scale to variance = 1
        ts.all_data(:,j) = (ts.all_data(:,j) - mean(ts.all_data(:,j))) ./ std(ts.all_data(:,j));
    end
    M2{i} = [ts.all_data];
end

end

fprintf(1,'Done.\n')

% ----------------------------------------------------------------------------------
% * Mean-center first (unscaled) timeseries and check for NaN or Inf values
% ----------------------------------------------------------------------------------
for i = 1:size(M1,2)
    M1(:,i) = M1(:,i) - mean(M1(:,i));
end

ex = any(isnan(M1) | isinf(M1)) | all(M1 == 0);
if any(ex)
	disp(['WARNING! NaN or Inf values for ' num2str(sum(ex)) ' voxels in timeseries!!  Mis-registration of funct and anatomy?'])
	M1(:,find(ex)) = [];
end

figure;imagesc(M1); title('Overall timeseries (y) for all voxels (x)')


% ----------------------------------------------------------------------------------
% * Find Principal Components overall
% ----------------------------------------------------------------------------------
fprintf(1,'\nPrincipal components overall ') 
[eigvec,eigval] = eig(M1'*M1);
eigvec = eigvec(:,end-2:end);
X = M1 * eigvec;
for i = 1:size(X,2), X(:,i) = (X(:,i) - mean(X(:,i))) ./ std(X(:,i));, end
figure;subplot 221; plot(eigvec); title('Weights on original variables (=voxels) == eigenvectors')
subplot 222; plot(X),title('Components (X * v)')
legend({'Comp 3' 'Comp 2' 'Comp 1'},0)
    subplot 223; bar(diag(eigval));, title('Scree plot for eigenvalues')
    
    xx = abs(fft(X)); xxx = (1:size(xx,1)) ./ (xX.RT * size(xx,1));
    subplot 224; plot(xxx(1:round(length(xx)./2)),xx(1:round(length(xx)./2),:));, title('FFT of components')
    xlabel('Frequency (Hz)')
drawnow

% ----------------------------------------------------------------------------------
% * Find Principal Components for each session
% ----------------------------------------------------------------------------------
wh = [0 cumsum(nimgs)];

for i = 1:nsess
    fprintf(1,'\nPrincipal components for session %1.0f',i) 

    ex = any(isnan(M2{i}) | isinf(M2{i}));
    if any(ex)	
	disp('')
	disp(['WARNING! NaN or Inf values for ' num2str(sum(ex)) ' voxels in sess ' num2str(i) '!!  Mis-registration of funct and anatomy?'])
	M2{i}(:,find(ex)) = [];
    end

    [eigvec,eigval] = eig(corrcoef(M2{i}));
    fprintf(1,'\t%3.0f eigenvalues > 1\t',sum(diag(eigval) > 1))
        
    num2save = min(sum(diag(eigval)>1),5);    % save at most 5 eigenvalues from this session
    eigvec = eigvec(:,end-num2save+1:end);
    X2 = M2{i} * eigvec;
    for j = 1:size(X2,2), X2(:,j) = (X2(:,j) - mean(X2(:,j))) ./ std(X2(:,j));, end
    
    figure;subplot 221; plot(eigvec); title(['Weights (eigenvectors) for session ' num2str(i)])
    subplot 222; plot(X2),title('Components (X * v)')
    subplot 223; bar(diag(eigval));, title('Scree plot for eigenvalues')
    
    xx = abs(fft(X2)); xxx = (1:size(xx,1)) ./ (xX.RT * size(xx,1));
    subplot 224; plot(xxx(1:round(length(xx)./2)),xx(1:round(length(xx)./2),:));, title('FFT of components')
    xlabel('Frequency (Hz)')
    
    drawnow
    
    % pad with zeros to get in the right session
    zbef = zeros(wh(i),size(X2,2));
    X2 = [zbef; X2];
    zaft = zeros(size(X,1) - size(X2,1),size(X2,2));
    X2 = [X2; zaft];
    
    X = [X X2];
end
    
X(:,end+1:end+length(xX.iB)) = xX.X(:,xX.iB); % add the nuisance covariates already in xX



% ----------------------------------------------------------------------------------
% * Correlate / Regress design vectors on components
% ----------------------------------------------------------------------------------
px = X * pinv(X);
fprintf(1,'\n Variation in predictors explained by nuisance covariates before orthogonalization')
fprintf(1,'\n Differences among predictors could create false activations without orthogonalization.')

for i = 1:length(xX.iC)
    r = xX.X(:,xX.iC(i)) - px * xX.X(:,xX.iC(i));       % residuals
    pve = 1 - ((r' * r) ./ (xX.X(:,xX.iC(i))' * xX.X(:,xX.iC(i)))); % percentage of variance explained
    fprintf(1,'\n%s\t%3.2f%%',xX.Xnames{xX.iC(i)},100*pve)
end



% ----------------------------------------------------------------------------------
% * Orthogonalize components (?) 
% This isn't a good idea from the standpoint of misattributing noise variance to signal
% but it will prevent artifactual activations based on differential correlations with
% nuisance subspace among conditions.  Lesser of two evils?
% ----------------------------------------------------------------------------------

if doortho
	fprintf(1,'\n Orthogonalizing nuisance covariates wrt model and scaling')
	mX = xX.X(xX.iC); mX(:,end+1) = 1;
	px = mX * pinv(mX);

	for i = 1:size(X,2)
    	X(:,i) = X(:,i) - px * X(:,i);                      % residuals
    	if ~(std(X(:,i)) == 0)
    		X(:,i) = (X(:,i) - mean(X(:,i))) ./ std(X(:,i));    % re-scale
    	end
	end

	figure; imagesc(X); colormap gray; title('Found and orthogonalized nuisance covariates')
	Xn = X; save Nuisance_covariates Xn
else
	figure; imagesc(X); colormap gray; title('Found nuisance covariates')
	Xn = X; save Nuisance_covariates_ortho Xn
end



% modify SPMcfg.mat
% just don't do this, it gets messy.  
% Better to just add them by loading Nuisance_covariates
% and entering as user-spec regressors.
% add_nuisance_to_SPMcfg(Xn);


fprintf(1,'\n Total running time is %3.2f s.\n',etime(clock,t1))
eval(['cd ' mypwd])

diary off
return

    
