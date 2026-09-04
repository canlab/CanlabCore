function sFIR_estimation(subIDs)
% A wrapper function for estimating voxel-wise smooth FIR model fits for each subject 
% in a multi-subject group, by Stephan Geuter.  It uses the hrf_fit method for fmri_data objects,
% which calls the HRF_est_Toolbox2 functions.  Future extensions could:
% - standardize the use of obj.images_per_session here, allowing flexibility in
% number of images per run
% - increase documentation
% - provide a complete walkthrough with data needed to run in our data
% repository


% folders
rdir   = '/work/ics/data/projects/wagerlab/labdata/projects/YOURPROJECT/'; % root directory
odir   = fullfile(rdir,'/data/onsets/'); % directory with onset files
adir   = fullfile(rdir,'/analysis/'); % directory for analysis outputs
imgdir = fullfile(rdir,'data','mri'); % directory hold subdirectories with fmri data/nifti

%%%

% design spec
TR        = 1.3; % TR in seconds
T         = TR*23; % total duration of estimated (s)FIR response
fl        = 'fl_D_cat_cue_on_sFIR'; % folder name for this analysis (sub directory of 'adir')

nimgs     = 1848; % n total images per subject
nrun      = 4; % n runs per subject
img_run   = nimgs/nrun; 


% filters
imgfilt   = 'swrarun_r*.nii';  % filter for nifti files of each subject
onsetfilt = 'cat_cue_on_stick.mat'; % filter for onset file for each subject

% analysis mask
maskfile  = fullfile(rdir,'masks','BasalGanglia','BG_T50.nii'); % fullpath to mask image
% standard mask: maskfile = which('brainmask.nii');


%%%

if ~isdir(fullfile(adir,fl)), mkdir(fullfile(adir,fl)); end


% loop subjects
nrun = numel(subIDs);
for crun = 1:nrun

    % current subject and folder
    fprintf(1,'\nnow running subject %d...\n\n',subIDs(crun));
    fldir = fullfile(adir,fl,sprintf('sub%04d',subIDs(crun)));
    if ~isdir(fldir), mkdir(fldir); end;
    cd(fldir);
    
    % filter all image files
    clear nimgs imgs
    imgs  = filenames(fullfile(imgdir,sprintf('sub%04d',subIDs(crun)),'run*',imgfilt));
    regnames = {};
    C = {};
    
    % filter onset file
    for o = 1:length(onsetfilt)
        
        % get onsets from subject specific onset directory
        ons  = filenames(sprintf('%ssub%04d/%s',odir,subIDs(crun),onsetfilt{o}),'char');
        load(ons); 
        
        % select regressors
        sel = find(cellregexp(names,'\w'));
        
        for k=1:length(sel)
            
            % make design cell array for HRF toolbox
            C{end+1} = zeros(nimgs,1);
            
            C{end}(ceil(onsets{sel(k)})) = 1;
            regnames{end+1} = names{sel(k)};
            %         plot(C{k}+k,'color',hsv2rgb([1/numel(osel)*k,1,1])); hold on
            %         fprintf(1,'\n%s\t\t%d',names{k},sum(C{k}));
        end
    end
    
    % load data within mask
    d     = fmri_data(imgs, maskfile);
    
    % HP filter
    [y, pI] = hpfilter(d.dat', TR, 160, repmat(img_run,nrun,1));
    
    % add intercept back
    y = y + pI * d.dat';
    
    % scale SPM style
    v  = [0 1:nrun-1] * 462;
    gs = [];
    for r = 1:nrun
        gs(r) = nanmean( nanmean( double(y(1+v(r):img_run+v(r),:)),2));
    end
    gs = repmat(gs,img_run,1); gs = gs(:)';
    gs = repmat(gs,size(d.dat,1),1);
    d.dat = y' * 100 ./ gs;
    
    % fit sFIR model. writes images in pwd (fldir)
    [params_obj hrf_obj] = hrf_fit(d, TR, C, T, 'FIR', 1);

    % save results in .mat files
    fn = fullfile(fldir,'sFIR_results.mat');
    save(fn,'params_obj','hrf_obj','onsetfilt','regnames','C','onsets');
    
end



end