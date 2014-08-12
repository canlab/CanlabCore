function [clusters,SPM,xX,xCon] = tor_extract_rois(imnames,varargin)
    % function [clusters, SPM, xX, xCon] = tor_extract_rois(imnames [can be empty],[opt] SPM, [opt] VOL, [opt] xX)
    %
    % this function gets timeseries data from all clusters in an SPM results output.
    % input:
    %	imnames: a matrix of image names, in spm_list_files output format
    %		 if empty, no timeseries data will be extracted.
    %
    %   clusters: if only 2 arguments, clusters structure is 2nd arg, and we
    %   extract data using existing clusters structure
    %
    %   If 3 arguments, enter SPM and VOL to extract data from VOXEL
    %   coordinates in these structures
    %	SPM:	 SPM variable from loaded results
    %	VOL:	 VOL variable from loaded results
    %
    %   Optional 4th argument fits a design matrix and returns betas
    %   xX:      xX design matrix to fit to timeseries
    %        OR 4th argument can be 0 (or any non-structure arg), suppressing verbose output.
    %
    %	[Last 2 arguments are optional.  Use if results are already loaded into workspace]
    %
    % Automatic fitting of model to cluster timeseries average using analyze_cluster_rois
    % with High-Pass filter length of your choice.
    % This only works if you input only the file names or input all optional arguments, including xX
    %
    % 10/17/01 by Tor Wager
    % Last modified 3/19/04 by Tor, to get clusters structure as input and use
    % timeseries3 instead of (bad on mac osx) timeseries2
    %
    % NOTE (WARNING): WORKS ON XYZ VOXEL COORDINATES - TRANSFORMATION TO
    % DIFFERENT SPACES ONLY IF ENTERING 2 ARGS, 1st one names, 2nd one clusters
    %
    % see transform_coordinates.m for transformation, or check_spm_mat, or cluster_interp.
    %
    % Functions called
    %   C:\matlabR12\toolbox\matlab\datatypes\squeeze.m
    %   c:\tor_scripts\voistatutility\nanmean.m
    %   (calls other spm functions)
    %   center_of_mass.m
    %   analyze_cluster_rois.m

    persistent prev_imnames
    persistent Vimnames
    
    verbose = 0;
    clusters = [];
    SPM = [];
    xX = [];
    xCon = [];
    clustersin = [];

    if nargin == 1
        % ----------------------------------------------------------------------------------
        % get SPM info from SPM.mat
        % ----------------------------------------------------------------------------------
        [SPM,VOL,xX,xCon,xSDM] = spm_getSPM;
    elseif nargin == 2
        clustersin = varargin{1};

        if isempty(clustersin), return, end
        
        if ~isempty(imnames)
            % ----------------------------------------------------------------------------------
            % check mat files for compatibility
            % -----------------------------------------------------------------
            if isfield(clustersin,'M')
                V = spm_vol(imnames(1,:));
                [chk,clustersin] = check_spm_mat(clustersin(1).M,V(1).mat,clustersin);
                if chk, disp('Warning!  Mat files for clusters and images do not match; using XYZmm to adjust cluster voxel coordinates'),end
            end
        end

        allxyz = cat(2,clustersin.XYZ);
    elseif nargin == 3
        SPM = varargin{1};
        VOL = varargin{2};
        if isempty(SPM), return, end
        if isempty(SPM.XYZ), disp('tor_extract_rois: no voxels to extract'), return, end
        allxyz = SPM.XYZ;
    elseif nargin == 4
        SPM = varargin{1};
        VOL = varargin{2};
        if isempty(SPM), return, end
        if isempty(SPM.XYZ), disp('tor_extract_rois: no voxels to extract'), return, end
        allxyz = SPM.XYZ;

        a = varargin{3};
        if isstruct(a), xX = a;
        else verbose = 0;
        end
    elseif nargin ~= 3
        error('Wrong number of arguments.  Use 1 if loading from SPM.mat, 3 args for no analysis, 4 with analysis.')
    end

    % ----------------------------------------------------------------------------------
    % get cluster index number from each voxel
    % ----------------------------------------------------------------------------------
    if isempty(clustersin)
        try
            %  WORKS IN SPM5 -- CHANGING from old version which did not
            %  cluster for large voxel sets
            cl_index = spm_clusters(SPM.XYZ);
            
%             if size(SPM.XYZ,2) > 60000
%                 warning('Too many voxels to cluster! SPM bug.');
%                 cl_index = ones(1,size(SPM.XYZ,2));
%             else
%                 cl_index = spm_clusters(SPM.XYZ);
%             end
        catch
            disp('Error evaluating: cl_index = spm_clusters(SPM.XYZ);')
            disp('No significant XYZ coordinates in SPM.XYZ!?')
            clusters = [];
            return
        end
    else
        cl_index = 1:length(clustersin);
    end


    % ----------------------------------------------------------------------------------
    % get the data from ALL voxels; all_data stores raw data for all
    % clusters
    % ----------------------------------------------------------------------------------
    O.coords = allxyz'; 
    clear allxyz
    if ~isempty(imnames)
        if ( any(size(prev_imnames) ~= size(imnames)) || any(prev_imnames(:) ~= imnames(:)) || length(Vimnames) ~= size(imnames, 1))
            prev_imnames = imnames;
            Vimnames = spm_vol(imnames);
        else
            disp('Using cached image names.');
        end
        all_data = spm_get_data(Vimnames, O.coords');
    end


    % ----------------------------------------------------------------------------------
    % define each cluster as cell in array.
    % ----------------------------------------------------------------------------------

    for i = 1:max(cl_index)
        if verbose && i == 1
            fprintf(1,'\n\t%3.0f clusters: Extracting %03d ',max(cl_index),i),
        elseif verbose
            fprintf(1,'\b\b\b%03d',i)
            if i == max(cl_index), fprintf(1,'\n'), end
        end

        % voxels in this cluster
        a = find(cl_index == i);


        % make cluster, if we don't have it
        if isempty(clustersin)

            if isfield(SPM,'title'), cl.title = SPM.title;  else  cl.title = '';  end
            if isfield(SPM,'u'),cl.threshold = SPM.u; else  cl.u = [];  end

            if isfield(SPM,'Z_descrip'), cl.Z_descrip = SPM.Z_descrip;  else  cl.title = '';  end

            if isfield(VOL,'VOX') && isfield(VOL,'M')
                cl.voxSize = VOL.VOX;
                cl.M = VOL.M;
            elseif isfield(VOL,'V') && isfield(VOL.V,'mat')
                cl.voxSize = diag(VOL.V.mat(1:3,1:3))';
                cl.M = VOL.V.mat;
            else
                error('Cannot find mat file info in VOL');
            end

            cl.name = [cl.title '_' num2str(i) '_' mat2str(size(a,2)) '_voxels'];
            cl.numVox = size(a,2);
            cl.Z = SPM.Z(a);
            cl.XYZmm = SPM.XYZmm(:,a);
            cl.XYZ = SPM.XYZ(:,a);
            if isfield(SPM,'df') && isfield(SPM,'STAT') && isfield(VOL,'R') && isfield(SPM,'n') && isfield(SPM,'u') && isfield(VOL,'FWHM')
                cl.pVoxelLev = spm_P(1,0,max(cl.Z),SPM.df,SPM.STAT,VOL.R,SPM.n);
                cl.pClustLev = spm_P(1,cl.numVox/prod(VOL.FWHM),SPM.u,SPM.df,SPM.STAT,VOL.R,SPM.n);
            end

            if ~isempty(cl.Z),
                % report number of sub-cluster peaks within cluster
                [N] = spm_max(cl.Z,cl.XYZ);
                cl.numpeaks = length(N);
            end

        else  % we already have clusters
            cl = clustersin(i);
        end

        % ----------------------------------------------------------------------------------
        % get the timeseries for the cluster
        % ----------------------------------------------------------------------------------
        if ~isempty(imnames)

            if isempty(clustersin)
                % no clusters input
                cl.all_data = all_data(:,a);

            else
                % clusters input
                nv = size(clustersin(i).XYZ,2);
                cl.all_data = all_data(:,1:nv);         % take first n data vectors
                all_data(:,1:nv) = [];                 % remove from all_data
            end



            cl.timeseries = nanmean(cl.all_data',1)';

            if size(imnames,1) < 200 && i == 1
                cl.imnames = imnames;
            elseif i == 1
                cl.imnames = imnames([1 end],:);
            else
                cl.imnames = [];
            end


            % ----------------------------------------------------------------------------------
            % get the SNR for the cluster if it looks like rfx data rather than individual data
            % ----------------------------------------------------------------------------------
            if size(imnames,1) < 60,
                try
                    cl.snr_avgts = get_snr(cl.timeseries);
                    cl.snr = get_snr(cl.all_data);
                    if verbose, fprintf(1,' : avg SNR = %3.2f, range %3.2f - %3.2f \n   ',cl.snr_avgts,min(cl.snr),max(cl.snr)), end

                    % calculate # subjects w/ estimates above zero
                    cl.numpos = sum(cl.timeseries > 0);

                    % calc 80% power level, no mult. comp corr
                    % from Jerry Dallal @ Tufts U.  http://www.tufts.edu/~gdallal/SIZE.HTM
                    % (16 * sd^2 / est^2) + 1
                    cl.power80 = (4 ./ cl.snr_avgts)^2 + 1;

                catch
                    warning('Problem with SNR calculation. Skipping.')
                end
            end

        else
            % disp('No timeseries data extracted - image names empty.')
            if isfield(SPM, 'all_data') && size(SPM.all_data, 2) == size(SPM.XYZmm, 2)
                cl.all_data = SPM.all_data(:, a);
                
                cl.timeseries = nanmean(cl.all_data, 2);
            end
            
        end

        cl.center = mean(cl.XYZ, 2)';

        % in case we've added a full matrix of Z-scores with
        % cluster_barplot, etc.
        if min(size(cl.Z)) > 1, cl.Z = ones(1,size(cl.XYZ,2));  end

        if size(cl.XYZmm,2) > 1, cl.mm_center = center_of_mass(cl.XYZmm,cl.Z);  % mean(cl.XYZmm');
        else cl.mm_center = cl.XYZmm';
        end
        clusters = [clusters, cl];
    end



    if ~isempty(imnames) && exist('xX') && ~isempty(xX)
        % ----------------------------------------------------------------------------------
        % adjust timeseries for each cluster and fit xX.X model to timeseries data
        % ----------------------------------------------------------------------------------
        try
            clusters = analyze_cluster_rois(clusters,xX);
        catch
            disp('Error analyzing timeseries clusters - skipping analysis.')
        end
    end
end

