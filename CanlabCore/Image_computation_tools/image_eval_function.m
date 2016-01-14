function varargout = image_eval_function(imageNames, fhandle, varargin)
% Evaluate any function, defined by the function handle fhandle,
% on each in-mask voxel for a set of images.
%
% :Usage:
% ::
%
%     varargout = image_eval_function(imgnames, fhandle, ['mask', mask], ['preprochandle', preprochandle], varargin)
%
% :Other Optional args: 'outimagelabels' , 'connames'
% ::
%
%     varargout = image_eval_function(Y, fhandle, varargin)
%
%
% evaluate fhandle on paired columns of X and Y
%
% :Note: You must call image_eval_function with outputs, one output for
% each output you're requesting from the voxel-level function.
% Eg:
% ::
%     [t, df, betaorcontrast, Phi, sigma, stebeta, F, pvals] = ...
%     image_eval_function(imgs, fhandle, 'mask', maskimg, 'outimagelabels' , names);
%
% :Inputs:
%
%   **fhandle:**
%        is a function handle:
%        ::
%            fhandle = @(variable inputs) fit_gls(variable_input, fixed_inputs);
%
%            fhandle = @(y) fit_gls(y, X, c, p, PX);
%
%
%   **'outimagelabels':**
%        should be followed by a set of image names, one name
%        per output of fhandle per element. e.g., outnames{j}{i} is the output
%        image for output j and element (input image) i.  elements may be images for each
%        subject, if a set of one image per subject is entered, or something
%        else depending on the nature of input imgs.
%
% Note: do not include suffixes: no .img
%
% :Examples: 
%
% Generalized least squares fitting on 100 Y-variables, same X
% ::
%
%    % Get image list
%    imgs = filenames('trial_height*img', 'char')
%    imgs = sort_image_filenames(imgs)
%
%    % Get pre-stored design matrix
%    X = eventdesign{3};
%
%    preprochandle = @(y) trimts(y, 3, []);
%
% Generate an image with the number of in-analysis (valid)
% subjects in each voxel
%
% EXPT.SNPM.P{2} is a list of subject-level contrast images.
% ::
%
%    fhan = @(y) sum(abs(y) > 0 & ~isnan(y));
%    y = image_eval_function(EXPT.SNPM.P{2}, fhan, 'mask', EXPT.mask, 'outimagelabels', {{'sum_valid_antic.img'}});
%
%    y = rand(100, 100); X = y + rand(100, 1); X(:,end+1) = 1; c = [1 0]'; p = 2; PX = pinv(X);
%    fhandle = @(y) fit_gls(y, X, c, p, PX);
%    [t, df, beta, Phi, sigma, stebeta, F] = fhandle(y);
%
% ..
%    tor wager, jan 31, 2007
% ..

    % setup inputs
    nout = nargout;
    setup_inputs(nout);
    connames;
    outimagelabels;

    % test data: Use to check and define outputs before running
    % ---------------------------------------------------------
    disp('Evaluating test data to get output structure.');
    ytest = randn(size(imageNames, 1), 1);
    if ~isempty(preprochandle), ytest = preprochandle(ytest);  end

    % must hard-code outputs; not ideal...
    out1 = []; out2 = []; out3 = []; out4 = []; out5 = []; out6 = []; out7 = []; out8 = []; out9 = []; out10 = [];
    out11 = []; out12 = []; out13 = []; out14 = []; out15 = []; out16 = []; out17 = []; out18 = []; out19 = []; out20 = [];

    % outputs are called out1, out2, ... etc.
    build_test_eval_string(nargout);
    eval(fstr);

    % put all outputs together in cell array
    outputs = [];
    for i = 1:nout
        eval(['outputs{i} = out' num2str(i) ';']);
    end

    % get names of images

    define_output_names();
    outimagenames;

    % Create or get memory-mapped volumes for output images
    % Vout is the spm_vol structure for each ouput image
    % ---------------------------------------------------------------------
    maskInfo = iimg_read_img(imageNames(1,:), 2);
    
    
    
    %This doesn't seem to be necessary.  iimg_reconstruct_vols does it.
    %Vout = create_output_images(maskInfo, outimagenames, nimgs_this_output);


    fprintf('Outputs: %3.0f images\n', nout);
    for ii = 1:nout
        if ~isempty(outimagenames{ii})
        fprintf('  %s \t- %3.0f volume(s)\n', outimagenames{ii}(1,:), nimgs_this_output(ii));
        end
    end
    fprintf('\n');

    % Prepare data matrix
    % ---------------------------------------------------------------------
    fprintf('Loading data...  ')
    tic
    [Y, maskInfo] = iimg_get_data(mask, imageNames);
    fprintf('  %3.0f sec\n', toc);



    % perform preprocessing function at each voxel
    % ---------------------------------------------------------------------
    [n, v] = size(Y);
    if ~isempty(preprochandle)
        fprintf('Preprocessing:\n')
        disp(preprochandle);
        updateiterations = 1:round(v ./ 100):v;
        updateperc = round(linspace(0, 100, length(updateiterations)));
        fprintf('Working ...')

        for i = 1:v
            Y(:,i) = preprochandle(Y(:,i));

            update = ( updateiterations == i );
            if any(update), fprintf('\b\b\b\b%3d%%', updateperc(update)); end
        end

        fprintf(' Done. \n')
    end

    fprintf('Running analysis on function: ')
    disp(fhandle);


    % fit function at each voxel
    % ---------------------------------------------------------------------

    % build string (fstr) that defines outputs
    build_eval_string(nargout);

    % must hard-code outputs; not ideal...
    out1 = []; out2 = []; out3 = []; out4 = []; out5 = []; out6 = []; out7 = []; out8 = []; out9 = []; out10 = [];

    % outputs are called out1, out2, ... etc.
    eval(fstr);

    % put all outputs together in cell array
    outputs = [];
    for i = 1:nargout
        eval(['outputs{i} = out' num2str(i) ';']);
    end

    % clear outputs to save memory
    clear out1 out2 out3 out4 out5 out6 out7 out8 out9 out10

    % Example of what this might look like for a specific application:
    %[t, df, beta, Phi, sigma, stebeta, F] = matrix_eval_function(Y, fhandle);
    %
    % Evaluated at each voxel: fhandle = @(y) fit_gls(y, X, c, p, PX);



    % reconstruct images and write
    % ---------------------------------------------------------------------



    % write them
    fprintf('Writing output images... ');
    tic

    write_output_images;

    fprintf('  %3.0f sec\n', toc);

    varargout = outputs;

    %END OF MAIN FUNCTION CODE


    % -------------------------------------------------------------------
    %
    %
    % INLINE FUNCTIONS
    %
    %
    %
    % -------------------------------------------------------------------

    function setup_inputs(nout)
        mask = [];
        preprochandle = [];
        outimagelabels = cell(1, nout);
        connames = {};
        startslice = 1;

        for i = 1:length(varargin)
            if ischar(varargin{i})
                switch varargin{i}

                    % functional commands
                    case 'mask', mask = varargin{i+1}; varargin{i+1} = [];
                    case 'data', dat = varargin{i+1}; varargin{i+1} = [];
                    case 'preprochandle', preprochandle = varargin{i+1}; varargin{i+1} = [];

                    case {'outnames', 'outimagenames', 'outimagelabels'}, outimagelabels = varargin{i+1}; varargin{i+1} = [];
                    case 'connames', connames = varargin{i+1}; varargin{i+1} = [];

                    case {'start', 'startslice', 'slice'}, startslice = varargin{i+1};

                    otherwise, warning(['Unknown input string option:' varargin{i}]);
                end
            end
        end

        if isempty(mask), error('mask is empty. no-mask opt. not implemented yet.');
        else  fprintf('Found mask: %s\n', mask);
        end

        if nout ~= length(outimagelabels)
            disp('Warning! Length of outimagelabels does not match number of outputs requested.');
        end
        
        imageNames = char(imageNames);
    end

    % % -------------------------------------------------------------------
    % % build string to evaluate that defines outputs and inputs
    % % -------------------------------------------------------------------
    function build_eval_string(numargs)
        fstr = '[';
        for arg = 1:numargs
            fstr = [fstr 'out' num2str(arg)];
            if arg ~= numargs
                fstr = [fstr ', '];
            else
                fstr = [fstr '] = matrix_eval_function(Y, fhandle);'];
            end
        end
    end

    function build_test_eval_string(numargs)
        fstr = '[';
        for arg = 1:numargs
            fstr = [fstr 'out' num2str(arg)];
            if arg ~= numargs
                fstr = [fstr ', '];
            else
                fstr = [fstr '] = fhandle(ytest);'];
            end
        end
    end

    % % -------------------------------------------------------------------
    % % print banner and define update points for text output
    % % -------------------------------------------------------------------
    function print_banner()
        fprintf('matrix_eval_function.m')
        fprintf('\n_______________________________\n')
        fprintf('Evaluating this function on %3.0f variables:\n', v);
        disp(fhandle);
        fprintf('...using this command:\n%s\n', fstr);
        fprintf('_______________________________\n')
        str = sprintf('Running ... Done  %03d%%', 0); fprintf(str);
        updateiterations = 1:round(v ./ 100):v;
        updateperc = round(linspace(0, 100, length(updateiterations)));
    end


    % % -------------------------------------------------------------------
    % % Define image output names
    % % -------------------------------------------------------------------
    % % %     function define_output_names
    % % %
    % % %         isok = [];
    % % %         % enter full set of output image names in outputimagelabels{i}{j}
    % % %         % check and see if conditions are met:
    % % %         for i = 1:nout
    % % %             if length(outimagelabels) >= i
    % % %                 if length(outimagelabels{i}) == size(outputs{i}, 2)
    % % %
    % % %                     for j = 1 : length(outimagelabels{i})
    % % %                         if ischar(outimagelabels{i}{j})
    % % %                             isok(i, j) = 1;
    % % %                         end
    % % %                     end
    % % %                 end
    % % %             end
    % % %         end
    % % %
    % % %         % if all is good, do nothing further; else, create
    % % %         if ~isempty(isok) && all(isok(:))
    % % %             disp('Using input image labels as-is');
    % % %             outimagenames = outimagelabels;
    % % %             return
    % % %
    % % %         else
    % % %             disp('Attempting to create output image names.');
    % % %         end
    % % %
    % % %         % enter stem for output names for each output in
    % % %         % outputimagelabels{x}, and names for each element
    % % %         % (subject/contrast/etc) image within each output in connames
    % % %         % or no connames for numbers
    % % %
    % % %         if ~exist('connames', 'var'), connames = []; end
    % % %         nout = length(outputs);
    % % %
    % % %         nc = length(connames);
    % % %
    % % %         outimagenames = cell(1, nout);
    % % %         for i = 1:length(outputs)
    % % %             % same number of output params as connnames?
    % % %             if nc == size(outputs{i}, 2)
    % % %                 for j = 1:nc
    % % %                     outimagenames{i}{j} = [outimagelabels{i} '_' connames{j} '.img'];
    % % %                 end
    % % %             elseif size(outputs{i}, 2) == 1
    % % %                 % only one, don't number
    % % %                 outimagenames{i}{1} = [outimagelabels{i} '.img'];
    % % %             else
    % % %                 for j = 1:size(outputs{i}, 2)
    % % %                     outimagenames{i}{j} = [outimagelabels{i} num2str(j) '.img'];
    % % %                 end
    % % %             end
    % % %         end
    % % %
    % % %         nimgs_this_output = ones(1, nout);
    % % %
    % % %         for ii = 1:nout
    % % %             nimgs_this_output(ii) = numel(outputs{ii});
    % % %         end
    % % %
    % % %     end


    function define_output_names()
        for i = 1:nout
            if length(outimagelabels) >= i  % if we have an output image name for this output
                outnm = outimagelabels{i};

                % Make sure we have .img extensions
                for j = 1:size(outnm, 1)  % for each name in this set, if there is more than one
                    this_img = deblank(outnm(j, :));

                    if ~strcmp(this_img(end-3:end), '.img')
                        % we need to add img extenstion
                        error('You must add .img to the end of each name.');

                        %outnm(j, :) = [outnm(j, :) '.img'];
                    end
                end

                %outimagenames{i} = outnm;
            end
        end

        nimgs_this_output = ones(1, nout);

        for ii = 1:nout
            nimgs_this_output(ii) = numel(outputs{ii});
        end

        for ii = 1:nout
            if nimgs_this_output(ii) > 1
                % if long enough, leave it alone
                if size(outimagelabels{ii}, 1) == nimgs_this_output(ii)
                    % we're OK, leave the names alone
                    outimagenames{ii} = outimagelabels{ii};

                elseif size(outimagelabels{ii}, 1) == 1

                    % we  have single name, but multiple output vols requested; expand into multi-volume set
                    spm_imgs = expand_4d_filenames(outimagelabels{ii}, nimgs_this_output(ii));

                    outimagenames{ii} = spm_imgs;

                else
                    % something is wrong.
                    disp('Image names and requested output lengths do not match up.')
                    disp('Either enter a string matrix defining image names for each element of each output,')
                    disp('or a single image name for the 4-D image for multiple-element outputs.');
                    error('Quitting.');

                end
            else
                outimagenames{ii} = outimagelabels{ii};
            end
        end
    end


    % % -------------------------------------------------------------------
    % % Write to disk
    % % -------------------------------------------------------------------
    function write_output_images
        nout = length(outputs);
        nc = length(connames);
        for i = 1:nout
            nimgs = size(outputs{i}, 2);
            if nimgs > 1, fprintf(' %03d', 0); end

            for j = 1:nimgs    % changed output file names from previous cell array to list

                if nimgs > 1, fprintf('\b\b\b\b %03d', j); end

                % Note: 12/5/2007; with SPM2 at least, the range of a
                % multi-volume image seems to be clipped for each image to
                % the range of the first volume when the image is created.
                % Therefore, it is important to create the images with the
                % appropriate range!
                % has something to do with spm scaling factor:
                % pinfo scale and offset seem to be set based on first
                % volume only.  creating first using spm_create_vol may not
                % help...?
                % if multiple images, re-create first image with range of
                % data across images...
                % kludgy fix...
                % Note: 4/1/08: This was an SPM scaling factor issue.
                % Images are now created by iimg_reconstruct_vols with
                % scaling factors of intercept = 0, slope = 1

                if nimgs > 1 && j == 1
                    % % %                     %mymax = max(abs(outputs{i}(:)));  % do the below to
                    % % %                     %avoid out of memory errors
                    % % %                     mymax = max(  max(max(outputs{i})), abs(min(min(outputs{i})))  );
                    tmp = outputs{i}(:,j)';
                    % % %                     tmp(1) = mymax; tmp(2) = -mymax;
                    iimg_reconstruct_vols(tmp', maskInfo, 'outname', deblank(outimagenames{i}(j,:)));

                    % code for checking stuff:
                    %i, j, tmp = outputs{i}(:,j)'; create_figure('plot'); plot(tmp), spm_image('init', deblank(outimagenames{i}(j,:))), global VV, VV = spm_vol(outimagenames{i}(j,:)); spm_type(VV.dim(4)) , spm_type(VV.private.hdr.dime.datatype), global vv, vv = spm_read_vols(VV); min(vv(:)), max(vv(:))
                else
                    iimg_reconstruct_vols(outputs{i}(:,j), maskInfo, 'outname', deblank(outimagenames{i}(j,:)));
                end


                % re-write first image, if needed
                if nimgs > 1 && j == nimgs
                    iimg_reconstruct_vols(outputs{i}(:,1), maskInfo, 'outname', deblank(outimagenames{i}(1,:)));
                end
            end
        end
    end
end





% % -------------------------------------------------------------------
% % Create new output images or check for existing ones
% % -------------------------------------------------------------------

function Vout = create_output_images(maskInfo, outimagenames, nimgs_this_output)

    Vout = {};

    for i = 1:length(outimagenames)

        if size(outimagenames{i}, 1) == nimgs_this_output(i)
            % One named image per output volume
        elseif size(outimagenames{i}, 1) == 1

        else
            % Mismatch!
            error('There seems to be an image name / output number mismatch.')
        end

        fprintf('Output %3.0f : %3.0f output volumes ', i, nimgs_this_output(i))

        for j = 1:size(outimagenames{i}, 1)
            % Named images or 4-D volumes (in which case, create first frame only)

            % % %         % enforce no filename extension; ext of type .img will be added
            % % %         % below
            % % %
            % % %         [dummy, outimagenames{i}, e] = fileparts(outimagenames{i}(1,:));
            % % %         namestr = ['.' filesep outimagenames{i} '.img'];

            namestr = deblank(outimagenames{i}(j, :));

            if exist(namestr, 'file')
                % image already exists; add to current
                fprintf('   ...Found existing: %s\n', namestr);

            else
                % don't need to create ALL images, cause will write whole-images
                % later...
                fprintf('   ...Creating: %s\n', namestr);

                Vout{i}{j} = make_output_image(maskInfo, namestr, ' ', 1); % single-volume only
            end
        end
    end

end

% % Create one output image
% % -------------------------------------------------------------------
function V = make_output_image(maskInfo, fname, descrip, n)
    V = struct('fname', '', 'dim', maskInfo.dim, 'mat', maskInfo.mat, 'pinfo', maskInfo.pinfo);

    %V.dim(4)  = spm_type(Type);

    % set data type to float
    switch(spm('Ver'))
        case 'SPM2'
            Type = 'float'; %'double';
            V.dim(4) = spm_type(Type);
        case 'SPM5'
            Type = 'float32';
            V.dt(1) = spm_type(Type); % NOT TESTED.  %'float32');
            V.dt(2) = maskInfo.dt(2);
        otherwise
            error('Unknown SPM version "%s": neuroscientists of the future, fix me!', spm('Ver'));
    end

    V.fname   = fname;
    V.descrip = descrip;
    V.n       = n;
    spm_create_vol(V);

    V.pinfo(1) = 1; % set image scaling factor to 1
    V.pinfo(2) = 0; % set image offset to 0
    V.pinfo(3) = prod(V.dim(1:3)) * spm_type(Type, 'bits') / 8 * (n-1);

    dat = NaN .* zeros(V.dim(1:3));
    spm_write_vol(V, dat);

end

% % % .img files can contain multiple volumes (indexed by .n in spm_vol structure).
% % % how many are in this volume?
% % % see also: scn_num_volumes.m
% % % % -------------------------------------------------------------------
% % function n = Nvol(V)
% %     % number of images n stored in this file
% %
% %     if ~isstruct(V)
% %         V = spm_vol(V);
% %         spm_close_vol(V);
% %     end
% %
% %     fp   = fopen(V.fname);
% %     fseek(fp, 0, 'eof');
% %     Len  = ftell(fp);
% %     fclose(fp);
% %
% %     switch(spm('Ver'))
% %         case 'SPM2'
% %             mydt = V.dim(4);
% %         case 'SPM5'
% %             mydt = V.dt(1);
% %         otherwise
% %             error('Unknown SPM version "%s": neuroscientists of the future, fix me!', spm('Ver'));
% %     end
% %     n    = Len/(prod(V.dim(1:3))*spm_type(mydt, 'bits')/8);
% % end

% .img files can contain multiple volumes (indexed by .n in spm_vol structure).
% how many are in this volume?
% % -------------------------------------------------------------------
function n = Nvol(V)
    % number of images n stored in this file

    if ~isstruct(V)
        V = spm_vol(V);
        %spm_close_vol(V); % spm2
    end

    switch(spm('Ver'))
        case 'SPM2'
            fp   = fopen(V.fname);
            fseek(fp,0,'eof');
            Len  = ftell(fp);
            fclose(fp);
            n    = Len/(prod(V.dim(1:3))*spm_type(V.dim(4),'bits')/8);

        case 'SPM5'
            n = length(V);

        otherwise
            error('Unknown SPM version "%s": neuroscientists of the future, fix me!', spm('Ver'));
    end

end

