% [clusters, cl_in_cell_struct] = extract_raw_data(EXPT, clusters, varargin)
%
% tor wager
% Extracts raw image data for each subject in each cluster
% averaging over cluster voxels (saved in clusters.all_data)
% and averaging over subjects (saved in clusters.timeseries)
%
% Options: High-pass filtering, spline detrending, principal components
%
% :Inputs:
%
%   **EXPT:**
%        as defined with get_expt_info.m
%   **clusters:**
%        as defined with tor_extract_rois.m
%        Uses EXPT.FILES.im_files to get raw image names
%
%        Optional: a string matrix of contrast or t-images to define
%        individually significant regions.  Empty input skips this step.
%
%        For high-pass filtering (recommended), use the 'hpfilter' option, and
%        enter HPDESIGN.TR, HPDESIGN.HP, and HPDESIGN.spersess = [number of images in each session], e.g., [220 220 220 220].
%        For session intercept fitting only, use EXPT.HP = Inf;
%        requires hpfilter.m
%
% :Optional Inputs:
%
% case 'extract_from', extract_imgs = varargin{i + 1}; varargin{i + 1} = [];
%
% case 'define_ind_rois', imgs = varargin{i + 1}; varargin{i + 1} = [];
% case 'subjects', subjidxs = varargin{i + 1};
%
% case 'noraw', doraw_data = 0;
% case {'dotrimts', 'windsorize'}, dotrimts = 0;
%
% case 'dospline', dospline = 1;
% case 'doprincomp', doprincomp = 1;
% case 'dospline', dospline = 1;
% case 'hpfilter', HPDESIGN = varargin{i + 1};
% case 'multiple_volumes', multiple_volumes = 1; % skip checking if files exist,
% because "exist" can't read comma-appended vols
%
% case 'load_and_continue', load_and_continue = 1; doraw_data = 1;
% case '4D', dim4 = 1; %use this for 4D data
%
% Example: Do 1st subject last, save results afterwards
% cl = extract_raw_data(EXPT, cl, 'subjects', [2:7 9:length(EXPT.FILES.im_files) 1]);
% save conjunction_cl_raw cl
%
% :Examples:
% ::
%
%    % Extract data from amygdala and put in convenient format:
%    % SETUP.data.M contains data images from a mediation analysis directory
%    % Amygdala mask images were created with SPM Anatomy toolbox
%    cm = mask2clusters('/Users/tor/Documents/matlab_code/3DheadUtility/SPM2_brains/ROI_spmanatomy_CMamy_MNI.img');
%    bl = mask2clusters('/Users/tor/Documents/matlab_code/3DheadUtility/SPM2_brains/ROI_spmanatomy_BLamy_MNI.img');
%    [cl, clcell] = extract_raw_data([], [cm bl], 'extract_from', SETUP.data.M, 'noraw');
%   amy = mediation_multilev_reformat_cl(clcell)
%
% ..
%    Modified March 9, 2008, Tor Wager; remove automatic spline detrending if no HP filter
%    Mod. May 23, 2008, Tor Wager; re-format input options
% ..

function [clusters, varargout] = extract_raw_data(EXPT, clusters, varargin)
    global defaults;

    dotrimts = 0;   % Windsorize voxel-by-voxel
    dospline = 0;   % spline detrend, if no HP filter
    doprincomp = 0; % get principal components from clusters
    doraw_data = 1; % save raw data in output
    multiple_volumes = 0; %can check if files exist, because none are comma-appended vols
    imgs = [];
    load_and_continue = 0; % pick up where we left off
    dim4 = 0; %assume not 4D images

    % --------------------------------------------------------
    % Set up inputs
    % --------------------------------------------------------

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % reserved keywords
                case 'extract_from', extract_imgs = varargin{i + 1}; varargin{i + 1} = [];

                case 'define_ind_rois', imgs = varargin{i + 1}; varargin{i + 1} = [];
                case 'subjects', subjidxs = varargin{i + 1};

                case 'noraw', doraw_data = 0;
                case {'dotrimts', 'windsorize'}, dotrimts = 0;

                case 'dospline', dospline = 1;
                case 'doprincomp', doprincomp = 1;
                case 'dospline', dospline = 1;
                case 'hpfilter', HPDESIGN = varargin{i + 1};
                case 'multiple_volumes', multiple_volumes = 1;

                case 'load_and_continue', load_and_continue = 1; doraw_data = 1;
                    
                case '4D', dim4 = 1;

                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end

    % Set up images to extract from
    if ~exist('extract_imgs', 'var')
        if isfield(EXPT, 'FILES') && isfield(EXPT.FILES, 'im_files')
            extract_imgs = EXPT.FILES.im_files;
        else
            disp('Cannot find either image names in extract_from input or EXPT.FILES.im_files.');
            error('Enter image names in cell array, ''extract_imgs'' followed by cells of image names (i.e., one per subject).');
        end
    end

    % Set up save file name
    if doraw_data && ~isfield(clusters, 'descrip')
        clusters(1).descrip = input('Enter descriptive name for these clusters (no spaces or special chars): ', 's');
    end

    % Set up raw data
    if doraw_data
        fname = [clusters(1).descrip '_raw_data.mat'];
    end

    % Set up continue
    if load_and_continue
        if exist(fname, 'file')
            disp(['Loading saved data: ' fname]);
        else
            fname = spm_get(0, '*mat', 'Choose raw_data.mat file to load or done to start fresh');
        end

        if ~isempty(fname) && exist(fname, 'file')
            load(fname)
        else
            error(['Cannot find saved file:' fname])
        end

    end

    % Set up subject indices
    N = length(extract_imgs);

    if ~exist('subjidxs', 'var'), subjidxs = 1:N; end


    % --------------------------------------------------------
    % try to get stuff for high-pass filtering
    % --------------------------------------------------------
    dohp = 0;
    TR = [];
    HP = Inf;
    spersess = [];

    if exist('HPDESIGN', 'var') && ~isempty(HPDESIGN)
        if isfield(HPDESIGN, 'TR'), TR = HPDESIGN.TR; end
        if isfield(HPDESIGN, 'HP'), HP = HPDESIGN.HP; end
        if isfield(HPDESIGN, 'spersess')
            spersess = HPDESIGN.spersess;
        elseif isfield(HPDESIGN, 'FIR') && isfield(HPDESIGN.FIR, 'nruns')
            spersess = HPDESIGN.FIR.nruns;
        elseif isfield(HPDESIGN, 'DX') && isfield(HPDESIGN.DX, 'nsess')
            spersess = HPDESIGN.DX.nruns;
        end
    end

    if ~isempty(TR) && ~isempty(HP) && ~isempty(spersess)
        dohp = 1;
        fprintf('Doing session intercept removal and high-pass filtering using %3.2f TR, %3.2f filter (s).\n', TR, HP);

        % set up filtering matrices for fast processing
        y = ones(sum(spersess), 1);   % dummy data
        [y, I, S] = hpfilter(y, TR, HP, spersess);
    else
        if dospline
            fprintf('HP filter info missing, using spline detrending every 100 images instead.\n');
        else
            fprintf('HP filter info missing (no HP); spline detrend option off.\n');
        end
    end

    if ~doprincomp
        fprintf('\t(principal components option off.)\n')
    end

    % --------------------------------------------------------
    % checks
    % --------------------------------------------------------

    switch spm('Ver')
        case 'SPM2'
            % spm_defaults is a script
            disp('WARNING: spm defaults not set for spm2. Make sure your defaults are set correctly');

        case 'SPM5'
            % spm_defaults is a function
            if(isempty(defaults))
                spm_defaults();
            end
    end

    % check clusters to see if SPM mat file info in clusters matches that of
    % data to extract.  Adjust voxel coordinates in clusters if necessary.
    if multiple_volumes == 0
        fprintf('Checking whether files exist...')
        for i = subjidxs
            for j = 1:size(extract_imgs{i}, 1)
                iname = deblank(extract_imgs{i}(j, :));
                if ~exist(iname, 'file')
                    disp(['Image does not exist! : ' iname]);
                    disp('Exiting extract_raw_data.m');
                    return
                end
            end
        end
    elseif multiple_volumes == 1
        fprintf('Checking whether files exist...')
        for i = subjidxs
            for j = 1:size(extract_imgs{i}, 1)
                iname = deblank(extract_imgs{i}(j, :));
                comma = find(iname(:)==',');
                iname = iname(1:comma-1);
                if ~exist(iname, 'file')
                    disp(['Image does not exist! : ' iname]);
                    disp('Exiting extract_raw_data.m');
                    return
                end
            end
        end
    end
    
    fprintf('Checking coordinates...')
    V = spm_vol(extract_imgs{1}(1, :));
    [chk, clusters] = check_spm_mat(clusters(1).M, V(1).mat, clusters);
    if chk, fprintf('\nAdjusted vox coords to extracted image space using XYZmm; new XYZ and XYZmm created.\nASSUMING ALL IMAGES HAVE SAME DIM/SPACE AS FIRST!\n'), end
    fprintf('\n')


    
    
    % get max length, set up output
    % -------------------------------------------------------- 
    nsubjects = max(subjidxs);
    nimgs = zeros(1, nsubjects);
    
    for i = subjidxs
        % get rid of empty image names
        wh = all(extract_imgs{i} == ' ', 2);
        extract_imgs{i}(logical(wh), :) = [];
        
        if ~dim4
            nimgs(i) = size(extract_imgs{i}, 1);
        else
            nimgs(i) = size(expand_4d_filenames(extract_imgs{i}), 1);
        end
        
    end
    maximgs = max(nimgs);
    if ~all(nimgs(subjidxs) == maximgs)
        fprintf('Warning! Not all image sets have same number of images\n');
        disp(nimgs(subjidxs))
    end
    
    % if we're not loading from file, initialize output vars
    if ~exist('ts', 'var')
        for j = 1:length(clusters)
            ts{j} = NaN * zeros(maximgs, nsubjects);
            
            if(doraw_data) && ~exist('raw', 'var')
                raw{j} = NaN * zeros(maximgs, size(clusters(j).XYZ, 2), nsubjects);
            end
        end
    end
    
    % --------------------------------------------------------
    % loop through subjects, extract, and process
    % --------------------------------------------------------


        
    for i = subjidxs
        t1 = clock;
        fprintf('Subject %3.0f  Extracting from %3.0f images...', i, nimgs(i))

        no_errors = 1;

        % Main data extraction
        % --------------------------------------------------------
        subjcl = tor_extract_rois(extract_imgs{i}, clusters);

        fprintf('%3.0f s. ', etime(clock, t1));

        %if nimgs(i) > 0, no_errors = 0; end
        
        if nimgs(i) ~= maximgs
            % Pad
            for j = 1:length(subjcl)
                subjcl(j).all_data = padwithnan(subjcl(j).all_data, ts{j}, 1);
            end
        end
        
% %         % check to see if this is the same length as first subject
% %         if i ~= 1 && nimgs(i) > size(ts{1}, 1)
% %             fprintf('Warning! Timeseries for subject %d is too long!!!\n', i);
% %             no_errors = 0;
% %         elseif i ~= 1 && nimgs(i) < size(ts{1}, 1)
% %             fprintf('Warning! Timeseries for subject %d is too short!!!\n', i);
% %             no_errors = 0;
% %             
% %         end

        if no_errors
            fprintf(' Processing: ')
            if dotrimts, fprintf(' Windsorizing. '); end

            t1 = clock;

            for j = 1:length(subjcl)
                % index so 1st is cluster, matrix of subjects
                if(doraw_data)
                    raw{j}(:, :, i) = subjcl(j).all_data;
                end
                subjcl(j).timeseries = [];

                % Windsorize and re-average
                tmp = subjcl(j).all_data;
                subjcl(j).all_data = [];

                for k = 1:size(tmp, 2)
                    if dotrimts
                        tmp(:, k) = trimts(tmp(:, k), 3, [], 1);     % last arg does sess break removal and spike correct
                    end

                    if dohp
                        tmp(:, k) = hpfilter(tmp(:, k), [], S, spersess, I);  % high-pass
                    elseif dospline
                        tmp(:, k) = splineDetrend(tmp(:, k));     % spline detrend every 100 images - remove LF drift
                    end
                end
                ts{j}(:, i) = nanmean(tmp', 1)';

                % eigenvectors and principal components
                if doprincomp
                    warning off                                 % divide by 0 warning if 1 voxel
                    [eigv, eigscore, eigval] = princomp(tmp);
                    ne = min(3, max(size(eigval)));              % number to save
                    eigs{j}(:, i, 1:ne) = eigscore(:, 1:ne);       % time x subjects x vectors
                    warning on
                else
                    eigs{j} = [];
                end

                % save names
                if i == 1 && j == 1
                    clusters(j).imnames = subjcl(j).imnames;
                end

                clusters(j).imP{i} =  extract_imgs{i};
                clusters(j).no_errors(i) = no_errors;
            end

            fprintf('%3.0f s. \n', etime(clock, t1));
        else
            fprintf('\tError in subject %d\n', i);
            for j = 1:length(subjcl)
                % bad subject, wrong # imgs
                if(doraw_data), raw{j}(:, :, i) = NaN; end
                ts{j}(:, i) = NaN;
                eigs{j}(:, i, :) = NaN;
                clusters(j).imP{i} =  extract_imgs{i};
                clusters(j).no_errors(i) = no_errors;
            end
        end

        if(doraw_data)
            save(fname, 'ts', 'raw', 'subjcl', 'eigs');
        end
        % %             save(fname, 'ts', 'subjcl', 'eigs');
        % %         end
    end % loop through subjects


    % --------------------------------------------------------
    % attach data to main clusters
    % --------------------------------------------------------

    for j = 1:length(clusters)
        ts{j}(ts{j} == 0) = NaN;

        clusters(j).all_data = ts{j};                   % one column per subject, average over voxels
        clusters(j).timeseries = nanmean(ts{j}')';      % average across subjects
        if doprincomp, clusters(j).prin_comps = eigs{j}; end
        if(doraw_data), clusters(j).raw_data = raw{j}; end                 % time x voxels x subjects
    end

    % --------------------------------------------------------
    % extract individual peaks, if contrast/tmap input is entered
    % --------------------------------------------------------

    if ~isempty(imgs)
        try
            clusters = extract_indiv_peak_data(clusters, imgs);
        catch
            disp('Error using extract_indiv_peak_data')
        end
    end
    
    % --------------------------------------------------------
    % Reformat to cell structure with all data (clpos_data format) if asked
    % for
    % --------------------------------------------------------

    if nargout > 0
        N = size(clusters(1).all_data, 2);
        clusters_cell = cell(1, N);


        for i = 1:N
            for j = 1:length(clusters)
                clusters_cell{i}(j) = clusters(j);
                clusters_cell{i}(j).imP = clusters_cell{i}(j).imP{i};
                clusters_cell{i}(j).timeseries = clusters_cell{i}(j).all_data(:, i);

                if isfield(clusters_cell{i}(j), 'raw_data')
                    clusters_cell{i}(j).all_data = squeeze(clusters_cell{i}(j).raw_data(:, :, i));

                else
                    clusters_cell{i}(j).all_data = [];
                end

            end

            if isfield(clusters_cell{i}, 'raw_data')
                clusters_cell{i} = rmfield(clusters_cell{i}, 'raw_data');
            end


        end

        varargout{1} = clusters_cell;
    end


end




function [chk, clusters] = check_spm_mat(mat1, mat2, clusters)
    %check_spm_mat(mat1, mat2, clusters)
    % mat1 is from clusters, mat2 is functional (imgs to extract)


    chk = mat1 - mat2;
    chk = chk(:);
    chk = chk(1:end-1);             % eliminate SPM scale factor
    chk = any(chk);

    if chk

        % we need to get the correct voxel coords from mat2 (funct) into
        % clusters, keeping the mm_coordinates the same!
        VOL.M = mat2;

        for i = 1:length(clusters)

            clusters(i).XYZ = mm2voxel(clusters(i).XYZmm, VOL)'; % functional img space, cluster mm coordinates

            clusters(i).Z = ones(1, size(clusters(i).XYZ, 2));
            clusters(i).XYZmm = voxel2mm(clusters(i).XYZ, mat2);

            clusters(i).M = mat2;

            clusters(i).voxSize = diag(clusters(i).M(1:3, 1:3)');

            % skip this, and clusters will have mm list of different length than
            % voxel list (XYZ).  data will be from voxel list.
            %SPM.XYZmm = voxel2mm(SPM.XYZ, VOL.M);    % slow, but gives unique voxels
        end
    end
end


