% out_names = write_residual_images(imgs, X, maskimg, varargin)
%
% Purpose: To create a set of residual images using a GLM fit.
% 
% INPUTS:
% imgs         image names, in string matrix
% X            the design matrix for regressors to remove, including
%               intercept
% maskimg      String; name of mask image.  Everything outside mask is
%               written as zero in output images. Does not have to be in
%               the same 'space' as the input; mask is resampled.
% Optional:
%   weights    followed by weights on observations
%   outputbase  followed by base name for output residual 4-D image
%               default: 'resid'
%
% ---------------------------------------------------------
% If you are creating your own whole-brain shell function 
% for your own application, here is the generic format:
%
% function my_function_name(imgs, fixed_inputs, varargin)
%
% imgs         image names, in string matrix
% fixed_inputs one or more inputs that should always be entered
% varargin     optional inputs. Entered as 'name', value pairs of input
%               arguments
%
% See comments in the code for more details.
%
% Will write the following files into the current directory:
% 
% Example:
% --------------------------------------------------------
% imgs = filenames('swtrial_AUC*.img', 'char', 'absolute')
% X = rand(size(imgs, 1), 3); X(:, 4) = 1;  % fake design matrix with 4 columns
% maskimg = which('brainmask.nii')
% mkdir test_residuals
% cd test_residuals
% out_names = write_residual_images(imgs, X, maskimg);

function out_names = write_residual_images(imgs, X, maskimg, varargin)

    % -------------------------------------------------------------
    % Define Optional Inputs Here
    % -------------------------------------------------------------
    
    % Default values for all optional inputs
    % ---------------------------------------
    nobs = size(imgs, 1);
    weights = ones(nobs, 1);        % default weights
    outputbase = 'resid';
    
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % functional commands
                case 'weights', weights = varargin{i+1};

                case {'outputbase', 'outname'}, outputbase = varargin{i+1};
                        
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end
    
    % Map mask image into space of imgs, if not already, and write
    % 'mask.img'
    % -------------------------------------------------------------
    scn_map_image(maskimg, imgs(1, :), 'write', 'mask.img');
    maskimg = 'mask.img';
    
    
    
    % =============================================================    
    % This next part processes the fixed inputs (same for each voxel)
    % and creates a function handle that can be run with any arbitrary
    % set of data (i.e., data from one voxel), and will produce 
    % whatever outputs you'd like to save in output images files.
    % Everything below is specific to your specific application.
    % =============================================================
    
    % Check weights
    % -------------------------------------------------------------
    wtol = max(sqrt(weights))*eps^(1/3);
    if any(weights == 0) || any(weights < wtol)
        fprintf('Some weights are too low or badly scaled. Setting zero weights to %3.4f\n.', wtol);
        weights(weights < wtol) = wtol;
    end

    W = diag(weights);


    % Setup inputs and function handle
    % -------------------------------------------------------------
    wh = intercept(X, 'which');
    if isempty(wh), disp('Warning! No intercept found in model.'); end

    %PX = pinv(X);
    PX = inv(X' * W * X) * X' * W;  % General weighted
    
    % Check if X is OK
    condx = cond(X);
    fprintf('Condition number of design matrix: %3.2f\n', condx);
    if condx > 10000, error('Design matrix is not estimable!'); end
    
    xpx = X * PX;
    
    % This is the function that will be run over the whole brain,
    % voxel by voxel
    % You could also define your own named function here, and call that
    % function.
    fhandle = inline('y - xpx * y', 'y', 'xpx');
    fhan = @(y) fhandle(y, xpx); 
    
    % In this example, y is the variable input (data), and everything else
    % is fixed (same for each voxel).
    %fhandle = @(y) fit_gls(y, X, contrasts, arorder, PX, weights);
    
    % =============================================================    
    % This next part generates output image names
    % In general, there is one set of names for each named output
    % in the function called by your function handle.
    % If only one output is returned (e.g., y), then one set of image names
    % need be specified. If the output is a scalar, one image will be created.
    % If the output is a vector, then you need to specify a list of names
    % for each 3-D volume. These can simply refer to the names of different
    % volumes in a single 4-D image. You can also specify a single file
    % name for a 4-D image, and multiple volumes will be created.
    %
    % If multiple outputs are returned, you need to specify one set of 
    % output names for each separate output. Each output gets its own cell.
    % 
    % For example, if your function returned 4 outputs, names might look
    % like this:
    %     names = cell(1, 4);
    %     names{1} = 'sigma.img';
    %     names{2} = 'phi.img';
    %     names{3} = 'df.img';
    %     names{4} = 'omnibus_F.img';
    % Some cells could have single output image names, and others could have
    % and lists of image names for vector outputs.
    %
    % =============================================================    

    
    % Setup output image names
    % This can be a single name for a 4-D file
    % -------------------------------------------------------------
    out_name = [outputbase '.img'];

    % =============================================================
    % it's often a good idea to save the input variables in a SETUP file so you
    % can refer to them later.
    % =============================================================

    
    save write_residual_images_SETUP imgs X weights outputbase maskimg

    % =============================================================
    % Last step: Run
    % =============================================================
    % This command actually has to have one named output for each output
    % you want to save in an image file.  SO this command is
    % application-specific.  For example, if you want to save a whole
    % series of image names (12 outputs with 12 images or sets of images), 
    % you might run something like this:
%     [beta, t, pvals, convals, con_t, con_pvals, sigma, Phi, df, stebeta, conste, F] = ...
%         image_eval_function(imgs, fhandle, 'mask', maskimg, 'outimagelabels' , names);
    
    y = image_eval_function(imgs, fhan, 'mask', [pwd filesep 'mask.img'], 'outimagelabels', {out_name});

    % This also returns names of output images.
    % it is application-specific for this function (not necessary):
    out_names = expand_4d_filenames(fullfile(pwd, out_name), nobs);
    
end
