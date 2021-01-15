function spm_htw_from_fit(varargin)
% When fitting a 1st-level SPM model with multiple basis functions,
% reconstructs fitted response at each voxel beta images and basis functions, and
% estimates height (amplitude), time-to-peak, and response width (duration and half-max)
% For each condition. Saves maps of these across voxels, and contrasts
% across height, time-to-peak, and width maps if contrasts are specified.
% This provides contrast images to take to 2nd-level analyses for group statistics when
% using multiple basis functions at the 1st level.
%
% :Usage:
% ::
%
%     function spm_htw_from_fit(varargin)
%
% :NOTES:
%
% This function loads the SPM.mat file in the current directory and uses
% the basis set specified in the loaded SPM structure.
%
% : Required Inputs:
% 	None - run within an SPM first-level analysis directory
% 
% :Optional Inputs:
%
%   **'amplitudes'**:
%       do amplitudes only; no contrasts
%
%   **'contrasts'**:
%       do contrasts only; no amplitudes
%
%   **'noamplitudes'**
%        skip create amping images (combination of betas
%        across basis functions)
%
%   **'nocontrasts'**
%        skip creation of contrast images
%        (running contrast images assumes that amplitude images are already
%        created)
%
%   **'all'**
%        will run both the amplitudes and contrasts sections
%
%          - The second way uses CANlab HTW code to estimate height, time to peak,
%            width, and area under the curve (see Lindquist & Wager 2007 for
%            simulations using a version of this code).
%            It requires SCANlab specific functions, in SCN_Core_Support
%            (unlike the deriv. boost).
%            To turn this OFF, enter 'nohtw' as an optional argument
%
%   **'startend'**
%        followed by starting and ending values in seconds for amplitude
%        estimation window (for HTW estimation only).
%        If you do not enter this, it will show you a plot and ask you to pick
%        these values.
%        If you enter them here as inputs, you can skip the interactive step and
%        loop over subjects.
%
%   **'condition_numbers'**
%        followed by which index numbers in your list should
%        be used to calculate h, t, w from.  You should use this if you
%        are entering regressors of no interest, besides the intercepts.
%
% :Important for Contrasts:
%    disp('Using contrasts entered as F-contrasts. Assuming the first contrast vector in each F-contrast '
%
%    disp('is a contrast of interest across the CANONICAL basis function regressors.')
%
% :Output:
% ::
% Reconstructed amplitude, time-to-peak, and duration (width) images for each event type
% e.g., for Event type (condition) 001:
% htw_amplitude_001.nii	
% htw_time_to_peak_001.nii
% htw_width_001.nii
% htw_area_under_curve_001.nii
% 
% Contrasts across these, using names stored in SPM.xCon.name:
% e.g., 
% con_htw_ampl_targetvsstandards5.nii
% con_htw_time_targetvsstandards5.nii	
% con_htw_widt_targetvsstandards5.nii
% con_htw_area_targetvsstandards5.nii	
%
% :Examples:
% ::
%
%    % RUN THIS IN COMMAND WINDOW TO BATCH
%    subj = dir('06*')
%    for i = 1:length(subj), cd(subj(i).name), spm_htw_from_fit, cd('..'); end
%
%    % ANOTHER BATCH EXAMPLE:
%    d = dir('remi*'); d = d(cat(2, d.isdir)); [mydirs{1:length(d)}] = deal(d.name)
%    for i = 1:length(mydirs), cd(mydirs{i}), spm_htw_from_fit('all', 'nodb', 'startend', [4 15]), cd('..'); end
%
%    %An example for an event-related design, specifying condition numbers to get HTW from:
%    spm_htw_from_fit('all', 'contrasts', 'condition_numbers', 1:14, 'startend', [4 10]);
%
%    % CALCULATE CONTRASTS ONLY ON ALREADY-ESTIMATED HTW IMAGES
%    spm_htw_from_fit('contrasts','condition_numbers',1:14);
%
% :References:
% ::
% Lindquist, M. A. & Wager, T. D. (2007). Validity and power in hemodynamic response modeling: 
% a comparison study and a new approach. Human Brain Mapping. 8:764-84.
%
% Lindquist, M. A., Meh Loh, J., Atlas, L. Y. & Wager T. D. (2009). Modeling the hemodynamic 
% response function in fMRI: efficiency, bias and mis-modeling. Neuroimage. 45:187-98.

% Notes: 
% This function is adapted from apply_derivative_boost.m, which was
% originally intended to implement Calhoun derivative boost, but this is
% now deprecated.
% Original documentation:
%        In addition, 'amplitudes' now has two separate parts:
%          - The first uses Vince Calhoun's derivative boost (Calhoun, 2004) to
%             estimate amplitudes.  NOTE: *We have not worked out the scaling yet, so
%             I'm not sure this is working right*
%             To turn this OFF, enter 'nodb' as an optional argument
%

    spmname = fullfile(pwd, 'SPM.mat');     % SETUP INPUTS
    if ~exist(spmname, 'file'), error('You must be in an SPM 1st-level results directory with SPM5 SPM.mat file.'); end
    
    load(spmname);
    
    nbf = size(SPM.xBF.bf, 2);
    
    docontrasts = 1;
    doamps = 1;
    nodb = 1;
    nohtw = 0;
    condition_numbers = [];
    do_downsample = [];          % downsample; default = 1 sec if units are seconds (1/dt)
    
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % reserved keywords
                case 'all', docontrasts = 1; doamps = 1;

                case 'amplitudes', doamps = 1; docontrasts = 0;

                case 'contrasts', docontrasts = 1; doamps = 0;

                case 'noamplitudes', doamps = 0;
                    
                case 'nocontrasts', docontrasts = 0;
                    
                case 'db', nodb = 0; % original DB estimation (legacy, deprecated)

                case 'nohtw', nohtw = 1; % skip HTW estimation

                case 'startend', startend = varargin{i + 1};  % starting and ending values in seconds for amplitude estimation window (for HTW estimation only).
                  
                case 'condition_numbers',  condition_numbers = varargin{i + 1}; 
                    
                case 'nodownsample', do_downsample = 0;
                    
                case 'downsample', do_downsample = varargin{i + 1};
                    
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end

    if isempty(do_downsample)
        % default downsampling
        do_downsample = round(1 ./ SPM.xBF.dt);
    end
    
    if ~(docontrasts || doamps)
        disp('Nothing to do! Enter ''contrasts'' ''amplitudes'' or ''all'' as input argument.');
        return
    end
    
    
    % ---------------------------------------------
    % FILE NAMES and imgtype
    % ---------------------------------------------
    
    imgs = dir(sprintf(['beta*img'])); imgs = char(imgs.name);
    
    if isempty(imgs)
        imgs = dir(sprintf(['beta*nii'])); imgs = char(imgs.name);
        imgtype = '.nii';
    else
        imgtype = '.img';
    end
    
    if isempty(imgs)
        error('Cannot find beta*img or beta*nii files in current folder')
    end
    
    n = size(imgs, 1);
    
        
    if doamps
        % ---------------------------------------------
        % ---------------------------------------------

        % ESTIMATE AMPLITUDES

        % ---------------------------------------------
        % ---------------------------------------------




        fprintf('Found %3.0f beta images', n); fprintf('\n');

        %load(spmname);
        nsess = length(SPM.Sess);
        fprintf('I think there are %3.0f sessions (runs)', nsess); fprintf('\n');

        if ~isempty(condition_numbers)
            wh_intercept = true(1, size(imgs, 1));  % exclude these
            wh_intercept(condition_numbers) = 0;
            wh_intercept = find(wh_intercept);

            fprintf('\nIncluding only these images: ');
            fprintf('%3.0f ', condition_numbers);
            fprintf('\n');

        else
            wh_intercept = SPM.xX.iB;
            fprintf('\nI think these images are intercepts, and am not using them: ');
            fprintf('%3.0f ', wh_intercept);
            fprintf('\n');
        end
        
        imgs(wh_intercept, :) = [];

        %not used; only when using image_eval_function
        %mask_img = './mask.img';

        n = size(imgs, 1);


        if n / nbf ~= round(n / nbf), error('Error!  Wrong number of images for the specified number of basis functions.'); end



        if nodb
            % skip Derivative Boost estimation and go straight to HTW
        else
            % DO DB
            
            switch nbf
                case 2
                    derivative_case = 'timeonly';
                case 3
                    derivative_case = 'timedispersion';
                    
                otherwise
                    warning('Deriv. Boost only works for SPM canonical hrf with time or time + dispersion derivatives.  This SPM.mat doesn''t match those specs.');
            end
            
            % Not used; only for image eval function
            % boost = @(b) sign(b(1)) .* sqrt(sum(b .^ 2));
            
            
    
            % ---------------------------------------------
            % CALCULATE
            % ---------------------------------------------
            cond_indx = 1;

            for i = 1: nbf : (n - nbf + 1)

                imgs_cond = imgs(i : i+nbf - 1, :);

                disp('Working on :')
                disp(imgs_cond)

                out_name = sprintf(['db_amplitude_%03d',imgtype ], cond_indx);

                % This code uses SCN lab tools to create images
                % ---------------------------------------------
                % %     y = image_eval_function(imgs_cond, boost, 'mask', mask_img, ...
                % %         'outimagelabels', {out_name});
                % ---------------------------------------------

                % This code uses SPM instead
                % ---------------------------------------------
                switch derivative_case
                    case 'timeonly'
                        spm_imcalc(imgs_cond, out_name, 'sign(i1) .* sqrt(i1.^2 + i2.^2)');
                    case 'timedispersion'
                        spm_imcalc(imgs_cond, out_name, 'sign(i1) .* sqrt(i1.^2 + i2.^2 + i3.^2)');
                    otherwise
                        error('Basis set is incompatible with DB estimation!');
                end

                fprintf('Created %s\n', out_name);
                % ---------------------------------------------

                cond_indx = cond_indx + 1;
            end

            % Get and save names
            db_amp_names = [];
            for i = 1:nsess
                sessnames = char(SPM.Sess(i).Fc.name);
                sessnames = [repmat(sprintf('Sess%02d_', i), size(sessnames, 1), 1) sessnames];
                db_amp_names = strvcat(db_amp_names, sessnames);
            end

            save db_amplitude_names db_amp_names
            disp(db_amp_names);
            disp(' ')
            disp('Saved DB amplitude condition names for each image in db_amplitude_names.mat');

            fprintf('\n*-----------------------------*\nApplied DB successfully\n*-----------------------------*\n')

        end

        if nohtw
            % skip this

        else

            % ---------------------------------------------
            % Estimated amplitude from fit: HTW
            % ---------------------------------------------
            % This code uses SCN lab tools to create images

            disp(' ')
            disp('Next: Estimating amplitude, time to peak, width, and area-under-curve images from fitted response using SCN lab code.')
            disp(' ')

            % downsample bf, if requested
            if do_downsample
                mytimeres = SPM.xBF.dt * do_downsample;
                SPM.xBF.bf = downsample(SPM.xBF.bf, do_downsample);
                
            else
                mytimeres = SPM.xBF.dt;
            end
           
            if exist('startend', 'var')
                % just check, and use input values
                if length(startend) ~= 2, error('Startend input must have two values, a starting and ending value in seconds for the amp. estimate window'), end
            
            else
                % Set range in sec
                htw_from_fit(SPM.xBF.bf, ones(size(SPM.xBF.bf, 2), 1), mytimeres, 'plot', 'verbose');

                disp(' ')
                disp('Enter the range in seconds within which to estimate peak amplitude.')
                disp('Example: type [3 12] and press return for a typical event-related setup.');
                disp('More sustained responses, like pain responses, may require a longer window.');
                disp('This estimates the amplitude of the IMPULSE RESPONSE, before convolution with the stimulus function')
                disp('so if you have an epoch design, a typical window of [3 12] sec is still appropriate.');
                disp('Also note: AUC images are calculated as the area under the curve within the window you specify.')
                disp(' ')
                disp('In the future, you can input ''startend'', [3 12] (for example) to skip interactive selection')
                startend = input('Enter your choice in [ ] and press return: ');
            end

            % Test your choice by showing you a plot
            htwfunction = @(b) htw_from_fit(SPM.xBF.bf, b, mytimeres, 'startval', startend(1), 'endval', startend(2), 'plot', 'verbose');
            htwfunction(ones(size(SPM.xBF.bf, 2), 1))
            drawnow

            % Create without plot option for loop through brain.
            htwfunction = @(b) htw_from_fit(SPM.xBF.bf, b, mytimeres, 'startval', startend(1), 'endval', startend(2));

            disp('Check the screen for a plot of your choice of window.')
            disp(' ')


            % ---------------------------------------------
            % CALCULATE
            % ---------------------------------------------
            cond_indx = 1;

            for i = 1: nbf : (n - nbf + 1)

                imgs_cond = imgs(i : i+nbf - 1, :);

                disp('Working on :')
                disp(imgs_cond)

                clear out_name
                out_name{1} = sprintf(['htw_amplitude_%03d',imgtype], cond_indx);
                out_name{2} = sprintf(['htw_time_to_peak_%03d',imgtype], cond_indx);
                out_name{3} = sprintf(['htw_width_%03d',imgtype], cond_indx);
                out_name{4} = sprintf(['htw_area_under_curve_%03d',imgtype], cond_indx);

                % ---------------------------------------------
                [h, t, w, auc] = image_eval_function(imgs_cond, htwfunction, 'mask', fullfile(pwd, ['mask', imgtype]), ...
                    'outimagelabels', out_name);

                h;t;w;auc;  % we need the outputs above to tell it to write 4 images.
                % ---------------------------------------------

                cond_indx = cond_indx + 1;
            end

            % Get and save names
            htw_amp_names = [];
            for i = 1:nsess
                sessnames = char(SPM.Sess(i).Fc.name);
                sessnames = [repmat(sprintf('Sess%02d_', i), size(sessnames, 1), 1) sessnames];
                htw_amp_names = strvcat(htw_amp_names, sessnames);
            end

            if ~exist(fullfile(pwd, 'db_amplitude_names.mat'), 'file')
                save db_amplitude_names htw_amp_names
            else
                save db_amplitude_names -append htw_amp_names
            end
            disp(htw_amp_names);
            disp(' ')
            disp('Saved HTW amplitude condition names for each image in db_amplitude_names.mat');

            fprintf('\n*-----------------------------*\nApplied HTW estimation successfully\n*-----------------------------*\n')

        end

    end  % amplitudes



    % ---------------------------------------------
    % ---------------------------------------------

    % CREATE CONTRAST FOR THIS SUBJECT

    % ---------------------------------------------
    % ---------------------------------------------

    if docontrasts

        % Load contrast vectors

        % -----------------------------------------
        %load(spmname);
        % nsess = length(SPM.Sess);

        % Define which indices to exclude from contrasts
        if ~isempty(condition_numbers)
            nconvals = length(SPM.xCon(1).c(:, 1));  % test contrast
            wh_intercept = true(1, nconvals);  % exclude these
            wh_intercept(condition_numbers) = 0;
            wh_intercept = find(wh_intercept);

            fprintf('\nIncluding only these images: ');
            fprintf('%3.0f ', condition_numbers);
            fprintf('\n');

        else
            wh_intercept = SPM.xX.iB;
            fprintf('\nI think these images are intercepts, and am not using them: ');
            fprintf('%3.0f ', wh_intercept);
            fprintf('\n');
        end


        if ~isfield(SPM, 'xCon')
            error('Enter F-contrasts first, with the first contrast vector in each F-contrast the contrast across the CANONICAL basis function.');
        else
            disp('Using contrasts entered as F-contrasts. Assuming the first contrast vector in each F-contrast ')
            disp('is a contrast of interest across the CANONICAL basis function regressors.')
        end

        wh_F = strmatch('F', char(SPM.xCon.STAT), 'exact');

        wh_T = strmatch('T', char(SPM.xCon.STAT), 'exact');
        
        wh_F = [wh_F wh_T];
        
        if isempty(wh_F) && isempty(wh_T), error('No F-contrasts or T-contrasts entered yet.'); end

        % All sets of contrast images
        % ------------------------------------------
        if ~nodb
            
            ampimgs_name = sprintf(['db_amplitude_*',imgtype]);
            
            ampimgs = dir(ampimgs_name); ampimgs = char(ampimgs.name);
            
            if ~isempty(ampimgs)
                contrast_image_names_dbamp = calc_contrasts(ampimgs_name, ampimgs, wh_F, SPM, wh_intercept, nbf, imgtype);
            end
            
        end
        
        
        % HTW amplitude
        % ------------------------------------------
        ampimgs_name = sprintf(['htw_amplitude_*',imgtype]);

        ampimgs = dir(ampimgs_name); ampimgs = char(ampimgs.name);

        if ~isempty(ampimgs)
            contrast_image_names_htwamp = calc_contrasts(ampimgs_name, ampimgs, wh_F, SPM, wh_intercept, nbf, imgtype);
        else
            disp(['Checked for but did not find: ' ampimgs_name]);
        end


        % HTW time
        % ------------------------------------------
        ampimgs_name = sprintf(['htw_time_to_peak_*',imgtype]);

        ampimgs = dir(ampimgs_name); ampimgs = char(ampimgs.name);

        if ~isempty(ampimgs)
            contrast_image_names_htwtime = calc_contrasts(ampimgs_name, ampimgs, wh_F, SPM, wh_intercept, nbf, imgtype);
        else
            disp(['Checked for but did not find: ' ampimgs_name]);
        end


        % HTW width
        % ------------------------------------------
        ampimgs_name = sprintf(['htw_width_*',imgtype]);

        ampimgs = dir(ampimgs_name); ampimgs = char(ampimgs.name);

        if ~isempty(ampimgs)
            contrast_image_names_htwwid = calc_contrasts(ampimgs_name, ampimgs, wh_F, SPM, wh_intercept, nbf, imgtype);
        else
            disp(['Checked for but did not find: ' ampimgs_name]);
        end


        % HTW area
        % ------------------------------------------
        ampimgs_name = sprintf(['htw_area_under_curve_*',imgtype]);

        ampimgs = dir(ampimgs_name); ampimgs = char(ampimgs.name);

        if ~isempty(ampimgs)
            contrast_image_names_htwarea = calc_contrasts(ampimgs_name, ampimgs, wh_F, SPM, wh_intercept, nbf, imgtype);
        else
            disp(['Checked for but did not find: ' ampimgs_name]);
        end
        
        check_it = whos('contrast_image_names*');
        
        if isempty(check_it)
            disp('No valid images found to create contrasts on');
            
        else
            save('db_amplitude_names', '-append', 'contrast_image_names*');
            disp('Saved lists of contrast image names in db_amplitude_names.mat');
        end
        
    end % end contrasts

    


    %% INLINE



end  % main function







function contrast_image_names = calc_contrasts(ampimgs_name, ampimgs, wh_F, SPM, wh_intercept, nbf, imgtype)


    n = size(ampimgs, 1);

    fprintf('Found %3.0f images:\n', n);
    disp(ampimgs)

    disp('Reading image data.')
    V = spm_vol(ampimgs);
    vols = spm_read_vols(V);

    % spm_check_registration(ampimgs);
    % colormap jet

    contrast_image_names = [];

    for i = 1:length(wh_F)


        name = SPM.xCon(wh_F(i)).name;
        disp(['Calculating contrast on: ' name])
        original_name = name;

        name = deblank(name);
        wh_bad = (name == ' ' | name == ',' | name == '.' | name == '^' | name == '~' | name == '''' | name == ':' | name == '*' | name == '%' | name == ';' | name == '@' | name == '&');
        name(wh_bad) = [];

        name = ['con_', ampimgs_name(1:8), '_', name, imgtype];

        c = SPM.xCon(wh_F(i)).c(:, 1);
        c(wh_intercept) = [];
        c = c(1 : nbf : end);

        if length(c) ~= n, error('Contrast is wrong length for some reason!  Coding error in this function?  Or wrong number of db_amplitude images.'); end

        fprintf('Contrast values: ')
        fprintf('%01d ', c)
        fprintf('\n')
        
        % calculate and save

        contrast_image_calc(name, c, vols, V, original_name)

        disp(['Written: ' name]);

        contrast_image_names = strvcat(contrast_image_names, name);

        disp(' ');

    end

    fprintf('\n*-----------------------------*\nContrasts Done successfully!\n*-----------------------------*\n')
    spm_check_registration(contrast_image_names);

end





function contrast_image_calc(Q, myc, vols, V, original_name)

    % FROM:
    % function contrast_image(Q, myc)
    %
    % Tor Wager
    %
    % Creates a contrast image called Q (do not include path)
    % Given a list of img files P (spm format, with path)
    % and a contrast vector myc
    % In the directory of 1st image in P.


    if ~(length(myc) == size(vols,4))
        error('Contrast vector length is not equal to number of image files.')
    end

    myc2 = zeros(size(vols));

    for i = 1:length(myc)
        myc2(:,:,:,i) = myc(i);
    end

    cvol = vols .* myc2;
    cvol = sum(cvol,4);

    % -------------------------
    % write
    % -------------------------
    dd = fileparts(V(1).fname);
    Q = fullfile(dd, Q);
    Vo = V(1);
    Vo.fname = Q;
    Vo.descrip = ['Contrast ' original_name];

    spm_write_vol(Vo,cvol);

end
