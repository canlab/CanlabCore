function image_eval_function_multisubj(imageNames,fhandle,varargin)
% Evaluate any function, defined by the function handle fhandle,
% on each in-mask voxel for a set of images.
% imageNames is a cell array of N cells, each containing images for one
% replication (i.e., subject)
%
% :Usage:
% ::
%
%     image_eval_function_multisubj(imgnames,fhandle,['mask',mask],['preprochandle',preprochandle],['outnames',outimagelabels],varargin)
%
% other optional args: 'outimagelabels' , 'connames'
%
% At each voxel, a cell array is formed, one cell per subject.
% This would correspond to a matrix is formed of t x N, where t is time and N is
% replication (subject) but cells can deal with unequal data vector lengths for each subject.
% The anonymous function in fhandle should operate on data from each cell (subject).
% ::
%
%    varargout = image_eval_function(Y,fhandle,varargin)
%
% evaluate fhandle on paired columns of X and Y
%
% fhandle is a function handle:
% ::
%
%    fhandle = @(variable inputs) fit_gls(variable_input,fixed_inputs);
%    fhandle = @(y) fit_gls(y,X,c,p,PX);
%
% specify the outputs by adding them as output image names.
% The number of outputs returned is determined by the number of named
% images in the list entered following the 'outnames' keyword.
%
% :Note on output images: You specify the names of the output images
% One image will be written per output of fhandle.
% The images will have one volume per element of the output variable.
% If you are returning an output with one value per subject, for
% example, then a single image will be written with one volume in it
% per subject.
%
% preprochandle is a function handle.
% it encapsulates the preprocessing function.
% the function should work on each cell (subject) of a t x N cell array of time courses for each
% subject, where each  cell contains a t x v matrix of data from a
% slice. The preproc function should thus be able to handle a whole slice as
% input.
% the function can itself be a cell array with multiple handles in
% different cells
%
% :Examples: Generalized least squares fitting on 100 Y-variables, same X
% ::
%
%
%    % Get image list
%    imgs = filenames('trial_height*img','char')
%    imgs = sort_image_filenames(imgs)
%
%    % Get pre-stored design matrix
%    X = eventdesign{3};
%
%    preprochandle = @(y) trimts(y,3,[]);
%
%    y = rand(100,100); X = y + rand(100,1); X(:,end+1) = 1; c = [1 0]'; p = 2; PX = pinv(X);
%    fhandle = @(y) fit_gls(y,X,c,p,PX);
%    [t, df, beta, Phi, sigma,stebeta, F] = fhandle(y);
%
% ..
%    tor wager, jan 31, 2007
% ..


    setup_inputs; % setup inputs
    connames;
    outimagenames;
    
    switch spm('Ver')
        case 'SPM2'
            % spm_defaults is a script
            disp('WARNING: spm defaults not set for spm2. Make sure your defaults are set correctly');

        case {'SPM5', 'SPM8', 'SPM12'}
            % spm_defaults is a function
            spm_defaults()

        otherwise
            % unknown SPM
            disp('Unknown version of SPM!');
            spm_defaults()
    end

    doload = 0;  % load key vars from file -- TEMPORARY

    % Prepare data images and map into volumes V
    % ---------------------------------------------------------------------
    disp('Preparing data.')
    
    % Make sure mask is in same dims as input images, and write mask.img to
    % current directory
    
    scn_map_image(mask, imageNames{1}(1,:), 'write', 'mask.img');
    maskInfo = iimg_read_img('mask.img', 2);
    

    fprintf('Mask info:\nDimensions: %3.0f %3.0f %3.0f\n', maskInfo.dim(1:3))
    fprintf('Voxels in mask: %3.0f\n', maskInfo.n_inmask)

    v =  maskInfo.n_inmask;

    N = length(imageNames);
    n = zeros(1, N);
    V = cell(1, N);

    if doload, load V_tmp, end

    fprintf('Mapping all input data volumes to memory\n')

    for ii = 1:N

        n(ii) = size(imageNames{ii}, 1);

        fprintf('Subject %3.0f, %3.0f images.\n', ii, n(ii));

        if ~doload,
            V{ii} = spm_vol(imageNames{ii});

            switch lower(spm('Ver'))

                case 'spm2'
                    % OK
                case 'spm5'
                    % update image number
                    n(ii) = length(V{ii});

                case 'spm8'
                    % update image number
                    n(ii) = length(V{ii});
                    
                otherwise
                    error('I don''t recognize your version of SPM!');

            end


        end


    end


    % set up preprocessing
    % ---------------------------------------------------------------------
    ynstr = {'Off' 'On'};
    fprintf(1,'Preprocessing is %s\n', ynstr{~isempty(preprochandle) + 1})
    if ~isempty(preprochandle)
        fprintf(1,'Preprocessing handle for subject 1:\n');
        disp(preprochandle{1})
    end


    % set up function evaluation
    % ---------------------------------------------------------------------

    % build string (fstr) that defines outputs
    fstr = build_eval_string(nout);

    % must hard-code outputs; not ideal...
    out1 = []; out2 = []; out3 = []; out4 = []; out5 = []; out6 = []; out7 = []; out8 = []; out9 = []; out10 = [];
    out11 = []; out12 = []; out13 = []; out14 = []; out15 = []; out16 = []; out17 = []; out18 = []; out19 = []; out20 = [];
    out21 = []; out22 = []; out23 = []; out24 = []; out25 = []; out26 = []; out27 = []; out28 = []; out29 = []; out30 = [];
    out31 = []; out32 = []; out33 = []; out34 = []; out35 = []; out36 = []; out37 = []; out38 = []; out39 = []; out40 = [];
    
    if nout > 40, error('Can only have up to 40 outputs.'); end

    % outputs are called out1, out2, ... etc.

    % Example of what this might look like for a specific application:
    %[t, df, beta, Phi, sigma,stebeta, F] = matrix_eval_function(Y,fhandle);
    %
    % Evaluated at each voxel: fhandle = @(y) fit_gls(y,X,c,p,PX);


    % do dummy test to get sizes of outputs
    % which are determined flexbily based on what function handle does.
    %
    % each element in output is stored in an image volume
    % elements in the same output (out1, out2, etc.) are stored as volumes
    % in the same .img file, indexed by .n fiels in the spm_vol structure
    % (these are stored during processing on 3rd dim of sliceoutput)
    %
    % to access a volume, try spm_vol('imagename.img, 3') for 3rd volume in
    % .img file.
    
    Y = cell(1, N);
    for ii = 1:N
        Y{ii} = rand(n(ii), 1);
    end

    if ~isempty(preprochandle)

        for ii = 1:N
            str = [num2str(ii) ' Preprocess Test'];
            fprintf('%s',str);

            if iscell(preprochandle)
                Y{ii} = preprochandle{ii}(Y{ii});
            else
                Y{ii} = preprochandle(Y{ii});
            end

            erase_string(str)
        end
    end

    eval(fstr);  % this evaluates the function

    for ii = 1:nout
        str = ['outputs{ii} = out' num2str(ii) ';'];
        eval(str)
    end

    % not necessary if we're writing multiple vols to same image
    %outimagenames = define_output_names;

    nimgs_this_output = ones(1, nout);

    for ii = 1:nout
        nimgs_this_output(ii) = numel(outputs{ii});
    end


    % Create or get memory-mapped volumes for output images
    % Vout is the spm_vol structure for each ouput image
    % ---------------------------------------------------------------------
    Vout = create_output_images(maskInfo, outimagenames, nimgs_this_output);

    fprintf('Outputs: %3.0f images\n', nout);
    for ii = 1:nout
        fprintf('  %s \t- %3.0f volume(s)\n', outimagenames{ii}, nimgs_this_output(ii));
    end
    fprintf('\n');



    disp('This will be evaluated at each voxel:')
    disp(fstr)
    disp(' ');



    % for each slice
    % ---------------------------------------------------------------------

    nZ = maskInfo.dim(3);

    %if doload
    %    which_slices = 3;
    %else
    which_slices = startslice : nZ;
    %end

    for whslice = which_slices

        % Load slice
        % ---------------------------------------------------------------------
        t1 = clock;

        clear Yslice
        [Yslice, wh_in_slice] = load_slice_data(whslice, V, maskInfo, preprochandle);

        fprintf(' Elapsed: %3.0f s\n', etime(clock, t1));
        vox_in_slice = sum(wh_in_slice);


        % get analyzed data for this slice (or empty values if no eligible voxels)
        % ---------------------------------------------------------------------
        sliceoutputs = analyze_slice_data;


        % write data from slice
        % ---------------------------------------------------------------------
        write_output_slice(Vout, sliceoutputs, whslice);

    end  % end slices


    %END OF MAIN FUNCTION CODE







    % -------------------------------------------------------------------
    %
    %
    % INLINE FUNCTIONS
    %
    %
    %
    % -------------------------------------------------------------------

    function setup_inputs
        mask = [];
        preprochandle = [];
        outimagenames = cell(1);
        connames = {};

        startslice = 1;

        for i = 1:length(varargin)
            if ischar(varargin{i})
                switch varargin{i}

                    % functional commands
                    case 'mask', mask = varargin{i+1}; varargin{i+1} = [];

                    case 'preprochandle', preprochandle = varargin{i+1}; varargin{i+1} = [];

                    case {'outnames' 'outimagenames'}, outimagenames = varargin{i+1}; varargin{i+1} = [];
                        nout = length(outimagenames);

                    case {'start','startslice','slice'}, startslice = varargin{i+1};

                    case 'connames', connames = varargin{i+1}; varargin{i+1} = []; % Not used anymore

                    otherwise, warning(['Unknown input string option:' varargin{i}]);
                end
            end
        end

        if isempty(mask), error('mask is empty. no-mask opt. not implemented yet.');
        else fprintf(1,'Found mask: %s\n',mask);
        end

    end





    % % -------------------------------------------------------------------
    % % Analyze data for this slice (inline)
    % % -------------------------------------------------------------------

    function sliceoutputs = analyze_slice_data

        updateiterations = 1:round(vox_in_slice ./ 100):vox_in_slice;
        updateperc = round(linspace(0,100,length(updateiterations)));


        fprintf(1,'Analysis ...')

        % initialize output images
        sliceoutputs = cell(1, nout);
        for ii = 1:nout
            sliceoutputs{ii} = NaN .* zeros([maskInfo.dim(1:2) nimgs_this_output(ii)]);
        end

        x = maskInfo.xyzlist(wh_in_slice, 1);
        y = maskInfo.xyzlist(wh_in_slice, 2);

        t1 = clock;

        % for each voxel in the slice
        for ii = 1:vox_in_slice

            % get data

            for subj = 1:N, Y{subj} = Yslice{subj}(:, ii);  end

            % fit function at this voxel
            % ---------------------------------------------------------------------

            eval(fstr)

            % save output in slice matrices
            % ---------------------------------------------------------------------
            for jj = 1:nout
                str = ['sliceoutputs{jj}(x(ii), y(ii), :) = out' num2str(jj) '(:);'];
                eval(str)
            end

            update = ( updateiterations == ii );
            if any(update), fprintf(1,'\b\b\b\b%3d%%',updateperc(update)); end

        end

        fprintf(' Elapsed: %3.0f s\n', etime(clock, t1));

    end




end     % END MAIN FUNCTION





% -------------------------------------------------------------------
%
%
% SUB-FUNCTIONS
%
%
%
% -------------------------------------------------------------------


% % -------------------------------------------------------------------
% % build string to evaluate that defines outputs and inputs
% % -------------------------------------------------------------------

function fstr = build_eval_string(numargs)
    fstr = '[';
    for arg = 1:numargs
        fstr = [fstr 'out' num2str(arg)];
        if arg ~= numargs
            fstr = [fstr ','];
        else
            fstr = [fstr '] = fhandle(Y);'];
        end
    end


end



% % -------------------------------------------------------------------
% % Load and preprocess (if specified) data for this slice
% % -------------------------------------------------------------------

function [Y, wh_in_slice] = load_slice_data(whslice, V, maskInfo, preprochandle)

    wh_in_slice = maskInfo.xyzlist(:,3) == whslice;
    xyzslice = maskInfo.xyzlist(wh_in_slice, :);
    xyzslice(:, end+1) = 1;
    xyzslice = xyzslice';

    fprintf('Slice %3.0f, %3.0f in-mask voxels\n', whslice, sum(wh_in_slice));
    fprintf('Loading and preprocessing for this slice for subject: ');

    N = length(V);

    Y = cell(1, N);

    if isempty(xyzslice)
        fprintf(' No in-mask data ');
        return
    end

    for i = 1:N
        % load
        str = [num2str(i) ' Load '];
        fprintf('%s',str);
        Y{i} = spm_get_data(V{i}, xyzslice);

        % preprocess
        erase_string(str)

        if ~isempty(preprochandle)

            str = [num2str(i) ' Preprocess '];
            fprintf('%s',str);

            if iscell(preprochandle)
                Y{i} = preprochandle{i}(Y{i});
            else
                Y{i} = preprochandle(Y{i});
            end

            erase_string(str)
        end


    end
    
    
    % Check for missing data in each set
    wh_bad_vox = cell(1, N);
    wh_missing_data = cell(1, N);
    
    fprintf('Summary of missing or bad voxels:\n');
    
    for i = 1:N
        
        wh_bad_vox{i} = all(Y{i} == 0) | all(isnan(Y{i}));
        wh_missing_data{i} = any(Y{i} == 0) | any(isnan(Y{i}));
        
        fprintf('Dataset %3.0f: %3.0f voxels have no data, and %3.0f voxels have missing data.\n', i, sum(wh_bad_vox{i}), sum(wh_missing_data{i}));
        
    end
    
        
    
end


% % -------------------------------------------------------------------
% % Create new output images or check for existing ones
% % -------------------------------------------------------------------

function Vout = create_output_images(maskInfo, outimagenames, nimgs_this_output)

    for i = 1:length(outimagenames)

        % enforce no filename extension; ext of type .img will be added
        % later.
        [dummy, outimagenames{i}, e] = fileparts(outimagenames{i});
        
        namestr = ['.' filesep outimagenames{i} '.img'];

        if exist(namestr, 'file')
            % image already exists; add to current
            fprintf(1,'   ...Found existing: %s\n', outimagenames{i});

            % check how many images are in it
            n = Nvol(namestr);

            for j = 1:nimgs_this_output(i)

                if j > n
                    % for any un-made volumes
                    fprintf('   ...adding volume %3.0f to %s\n', j, outimagenames{i});
                    Vout{i}{j} = make_output_image(maskInfo, outimagenames{i}, ' ', j);
                else
                    Vout{i}{j} = spm_vol(sprintf('%s, %d',namestr, j));
                end

            end


        else

            
            fprintf(1,'   ...Creating: %s\n', outimagenames{i});
            for j = 1:nimgs_this_output(i)
                % images stored in this file.  Each is stored using 'n' field
                % in same 4-D .img file.
                Vout{i}{j} = make_output_image(maskInfo, outimagenames{i}, ' ',j);
            end
        end
    end

end

% % Create one output image
% % -------------------------------------------------------------------
function V = make_output_image(maskInfo, fname, descrip,n)
    

    V = struct('fname','', 'dim', maskInfo.dim, 'mat',maskInfo.mat, 'pinfo', [1 0 maskInfo.pinfo(3, 1)]'); % scaling factors 1 0, keep data type from mask 

    % set data type to float
    switch(spm('Ver'))
        case 'SPM2'
            Type = 'double';
            V.dim(4) = spm_type(Type);
        case {'SPM5', 'SPM8'}
            Type = 'float32';
            V.dt(1) = spm_type(Type); % NOT TESTED.  %'float32');
            V.dt(2) = 1;
        otherwise
            error('Unknown SPM version "%s": neuroscientists of the future, fix me!', spm('Ver'));
    end
    
    V.fname   = [fname '.img'];
    V.descrip = descrip;
    V.n       = n;
    spm_create_vol(V);

    V.pinfo(3)=prod(V.dim(1:3))*spm_type(Type,'bits')/8*(n-1);

    dat = NaN .* zeros(V.dim(1:3));
    spm_write_vol(V, dat);
    
    % re-load to get full volinfo structure with private (Nifti), for this
    % 3-D volume only
    V = spm_vol([V.fname ', ' num2str(n)]);

end

% % % .img files can contain multiple volumes (indexed by .n in spm_vol structure).
% % % how many are in this volume?
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
% %     fseek(fp,0,'eof');
% %     Len  = ftell(fp);
% %     fclose(fp);
% %     n    = Len/(prod(V.dim(1:3))*spm_type(V.dim(4),'bits')/8);
% % 
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

        case {'SPM5', 'SPM8'}
            n = length(V);

        otherwise
            error('Unknown SPM version "%s": neuroscientists of the future, fix me!', spm('Ver'));
    end

end



% % -------------------------------------------------------------------
% % Write to disk
% % -------------------------------------------------------------------
function write_output_slice(Vout, sliceoutputs, whslice)

    fprintf(1,'Writing slice %3.0f in output images.\n', whslice);
    
    for i = 1:length(Vout)  % length(sliceoutputs)

        scn_write_plane(cat(2, Vout{i}{:}), sliceoutputs{i}, whslice);

        % % %         for j = 1:length(Vout{i})  % size(sliceoutputs{i},3)    % for each image in 3rd dim
        % % %
        % % %             myslice = sliceoutputs{i}(:,:,j);
        % % %
        % % %             spm_write_plane(Vout{i}{j}, myslice, whslice);
        % % %
        % % %             %iimg_reconstruct_3dvol(myslice(:), maskInfo, 'outname', outimagenames{i}{j}, 'slice', whslice);
        % % %
        % % %         end
    end
end





% OLD inline functions no longer used

% % -------------------------------------------------------------------
% % print banner and define update points for text output (inline)
% % -------------------------------------------------------------------
%     function print_banner
%
%         fprintf(1,'matrix_eval_function.m')
%         fprintf(1,'\n_______________________________\n')
%         fprintf(1,'Evaluating this function on %3.0f variables:\n',v);
%         disp(fhandle);
%         fprintf(1,'...using this command:\n%s\n',fstr);
%         fprintf(1,'_______________________________\n')
%         str = sprintf('Running ... Done  %03d%%',0); fprintf(1,str);
%         updateiterations = 1:round(v ./ 100):v;
%         updateperc = round(linspace(0,100,length(updateiterations)));
%     end
%
%
%     % % -------------------------------------------------------------------
%     % % Define image output names (inline)
%     % % -------------------------------------------------------------------
%     function outimagenames = define_output_names
%
%         %%% OLD: not used now
%         %%% outimagelabels should be one per output
%         % names would be created for each image
%         % this has been replaced by a scheme which writes multiple vols to
%         % same image file
%
%         % outputs for one voxel; for formatting
%         outputs = cell(1, nout);
%         for i = 1:nout, str = ['outputs{i} = out' num2str(i) ';']; eval(str); end
%
%         if ~exist('connames','var'), connames = []; end
%
%         nc = length(connames);
%
%         outimagenames = cell(1,nout);
%         for i = 1:nout
%             % same number of output params as connnames?
%             if nc == size(outputs{i},2)
%                 for j = 1:nc
%                     outimagenames{i}{j} = [outimagelabels{i} '_' connames{j} '.img'];
%                 end
%             elseif size(outputs{i},2) == 1
%                 % only one, don't number
%                 outimagenames{i}{1} = [outimagelabels{i} '.img'];
%             else
%                 for j = 1:size(outputs{i},2)
%                     outimagenames{i}{j} = [outimagelabels{i} num2str(j) '.img'];
%                 end
%             end
%         end
%
%     end
