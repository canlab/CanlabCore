function cluster_orthviews(varargin)
% :Usage:
% ::
%
%    cluster_orthviews(inputs in any order)
%
% This function uses spm_orthviews to display activation blobs from a
% clusters structure on a canonical structural brain image. Multiple
% clusters may be plotted in multiple colors, and blobs may be added to an
% existing orthviews figure.
%
% :Inputs:
%
%   clusters structures
%
%   colors cell array, e.g., {[0 0 1] [1 0 0] [0 1 0]}
%   - if no colors cell array, uses Z- or t-scores to color map
%   - if colors cell is same length as clusters structure, uses these
%   colors as fixed colors for each region
%
% :Optional Inputs:
%
%   **'add':**
%        to suppress making new orthviews
%
%   **'copy':**
%        to copy to smaller subfigure with empty axis beside it
%
%   **'unique':**
%        to display in unique colors for each cluster
%
%   **'overlay':**
%        followed by name of image, to use a custom anatomical overlay
%
%   **'bivalent':**
%        to plot increases in solid colors specified by colors cell
%        {1} and {2}
%
%        OR, if no colors entered, use hot/cool map in 
%        spm_orthviews_hotcool_colormap
%
% :Options if you specify colors (addColouredBlobs):
%
%   **'trans':**
%        to make blobs transparent
%
%   **'solid':**
%        to make them solid
%
% :Options for using spm's color map (addBlobs):
%
%   **'blue':**
%        display in blue split color map instead of default red-yellow
%
%
%   **'handle':**
%        followed by integer for which orthviews window to use (default = 1)
%
% :Options for matching colors across left/right hemispheres (region objects only):
% all_colors = match_colors_left_right(r);
% orthviews(r, all_colors);


% Programmers' notes
%
% ..
%    By Tor Wager
%    April 2007 by TW : add hot/cool colormap to bivalent option
%               Add handles option
%    April 2014 by TW : add custom colors option for unique colors
%    Jan 2020: changed to spm('Defaults','fmri')
% ..

    spm('Defaults','fmri')
    
    % overlay = which('scalped_single_subj_T1.img');
    % overlay = which('SPM8_colin27T1_seg.img');  % spm8 seg cleaned up
%     overlay = which('keuken_2014_enhanced_for_underlay.img');
                overlay = which('spm152.nii');

    donew = 1;  % new fig
    docopy = 0;  % copy to new axis
    douniquecolors = 0;
    dobiv = 0;
    dotrans = 0;
    dosolid = 0;
    doblue = 0;
    initwhitebg = 0; % automatic white background for spm_orthviews
    skipempty = 0;
    
    wh_handle = 1;
    
    cl = [];
    cols = [];

    for i = 1:length(varargin)
        if iscell(varargin{i})
            
            cols = varargin{i}; % colors
            
        elseif isstruct(varargin{i}) || strcmp(class(varargin{i}), 'region') 
            
            cl{end+1} = varargin{i};
            
        elseif ischar(varargin{i})
            
            switch(varargin{i})
                case 'add'
                    donew = 0;
                case 'copy'
                    docopy = 1;
                case 'unique'
                    douniquecolors = 1;
                case 'overlay'
                    overlay = deblank(varargin{i+1});
                case 'bivalent'
                    dobiv = 1;
                case 'trans'
                    dotrans = 1;
                case 'solid'
                    dosolid = 1;
                    
                case 'blue'
                    doblue = 1;
                    
                case 'handle'
                    wh_handle = varargin{i+1};
                
                case 'skipempty'
                    skipempty = 1;    
                    
                case {'whitebg', 'white', 'initwhitebg'}
                    initwhitebg = 1;

            end
        end
    end

    if ischar(overlay)
        [path, name, ext] = fileparts(overlay);
        if(isempty(path))
            overlay = which(overlay);
        end
    elseif isempty(overlay)
        overlay = which('spm2_single_subj_T1_scalped.img');  % spm2
    end

    if isempty(overlay)
        
        disp('No valid overlay image specified, and I cannot find the default one.');
        error('The default is SPM8_colin27T1_seg.img, and should be on your path.');
        
    end
    
    if isempty(cl)
        if skipempty
            % do nothing
        else
             spm_check_registration(overlay); 
        end
        
        if initwhitebg, doinitwhitebg; end
        return
    end

    for i = 1:length(cl)
        if dotrans, cl{i} = setup_trans(cl{i}); end
        if dosolid, cl{i} = setup_solid(cl{i}); end
    end

    for i = 1:length(cl)
        CLU{i} = clusters2CLU(cl{i});
        if ~isfield(CLU{i}, 'Z'), CLU{i}.Z = ones(1, size(CLU{i}.XYZ, 2)); end
        if min(size(CLU{i}.Z)) > 1, CLU{i}.Z = CLU{i}.Z(1, :); end
    end



    % -----------------------------------------------------------
    % Bivalent response: red = activation, blue = deactivation
    % -----------------------------------------------------------
    if dobiv

        if ~isempty(cols)
            % we have solid colors (addColoredBlobs, not the heat-mapping)
            [poscl, negcl] = cluster_separate_posneg(cl{1}, 'z=1');

            if length(cols) < 2
                cols = [cols {[0 0 1]}];
            end
            
            %view clusters
            cluster_orthviews(poscl, cols{1}, 'overlay', overlay);
            if ~isempty(poscl)
                cluster_orthviews(negcl, cols{2}, 'add', 'overlay', overlay);
            else
                cluster_orthviews(negcl, cols{2}, 'overlay', overlay);
            end

        else
            % we have color mapping mapping
            cluster_orthviews(cl{1}, 'overlay', overlay);
            
            try
                Zvals = cat(2,cl{1}.Z); 
            catch
                error('Invalid clusters.Z field'); 
            end
            spm_orthviews_hotcool_colormap(Zvals, min(abs(Zvals)) .* .95);

        end

        return
    end


    if donew, warning off, spm_check_registration(overlay), set(gcf, 'Resize', 'on'); warning on, end %spm_image('init', overlay), end

    if length(cols) == length(cl{1})
        douniquecolors = 1;
    end
    
    if douniquecolors
        % unique colors for each blob
        
        if length(cols) == length(cl{1}) % use input colors
            colors = cols; 
        else
            colors = scn_standard_colors(length(cl{1})); % generate colors
        end

        ind = 1;
        for i = 1:length(cl)
            while length(colors) < length(cl{i})
                colors = scn_standard_colors(length(cl{i}));
                %colors = [colors {rand(1, 3)}];
            end
            for j = 1:length(cl{i})
                CLUtmp = clusters2CLU(cl{i}(j));
                if ~isfield(CLUtmp, 'Z'), CLUtmp.Z = ones(1, size(CLUtmp.XYZ, 2)); end
                if min(size(CLUtmp.Z)) > 1, CLUtmp.Z = CLUtmp.Z(1, :); end
                spm_orthviews('AddColouredBlobs', wh_handle, CLUtmp.XYZ, CLUtmp.Z(1, :), CLUtmp.M, colors{ind})
                ind = ind + 1;
            end
        end
    else
        for i = 1:length(cl)
            if isempty(cols) || length(cols) < i
                spm_orthviews('AddBlobs', wh_handle, CLU{i}.XYZ, CLU{i}.Z, CLU{i}.M)
                
                if doblue
                    % seems to work when we copy into window, but
                    % not in this script...
                    myfig = findobj('Tag','Graphics');
                    figure(myfig);
                    cm = get(myfig,'Colormap');
                    cm(65:end,:) = cm(65:end,[3 2 1]);
                    colormap(cm);
                end

            else
                spm_orthviews('AddColouredBlobs', wh_handle, CLU{i}.XYZ, CLU{i}.Z(1, :), CLU{i}.M, cols{i})
            end
        end
    end


    if docopy
        hh = gcf;
        h1 = get(hh, 'Children');
        f1 = figure('Color', 'w');
        c = copyobj(h1, f1);

        set(gcf, 'Position', [210 322 929 407]);
        for i = 1:length(c)
            copypos = get(c(i), 'Position');
            copypos(3) = copypos(3) .* .5;
            set(c(i), 'Position', copypos);
        end
        keyboard
    end

    mypos = mean(cl{1}(1).XYZmm, 2);
    if length(mypos) < 3, mypos = cl{1}(1).XYZmm; end
    if isempty(mypos), mypos = [0 0 0]'; end
    spm_orthviews('Reposition', mypos);

    % try to set the window button up function to show x, y, z position
    % coordinates
    fh = findobj('Tag', 'Graphics');
    if isempty(get(fh, 'WindowButtonUpFcn'))
        set(fh, 'WindowButtonUpFcn', 'spm_orthviews_showposition;');
        spm_orthviews_showposition;
        
    elseif strcmp(get(fh, 'WindowButtonUpFcn'), 'spm_orthviews_showposition;')
        % we already have it; do nothing, and no warning
    else
        % warning
        disp('A windowbtnup function already exists; not showing position coordinates.')
    end

    
    if initwhitebg, doinitwhitebg; end

end % main function







% ---------------------------------------------------------------------
%
%
% Sub-functions
% (Needed only for bivalent option)
%
% ---------------------------------------------------------------------


function [poscl, negcl] = cluster_separate_posneg(cl, varargin)
    clu = clusters2CLU(cl);

    whpos = find(clu.Z > 0);
    whneg = find(clu.Z < 0);
    N = fieldnames(clu);


    if length(varargin) > 0
        switch varargin{1}
            case 'z=1'
                clu.Z = ones(size(clu.Z));
        end
    end

    poscl = [];
    negcl = [];
    for i = 1:length(N)
        if length(clu.(N{i})) == length(clu.Z)
            poscl.(N{i}) = clu.(N{i})(:, whpos);
            negcl.(N{i}) = clu.(N{i})(:, whneg);
        else
            poscl.(N{i}) = clu.(N{i});
            negcl.(N{i}) = clu.(N{i});
        end
    end

    poscl = tor_extract_rois([], poscl, poscl);
    negcl = tor_extract_rois([], negcl, negcl);
end



function cl = setup_trans(cl)
    z = cat(2, cl.Z); z = z(:);
    %sd = nanstd(z);
    %z = z./sd;
    %for i = 1:length(cl)
    %    cl(i).Z = cl(i).Z ./ sd;
    %end
    cl(1).Z(1) = max(z) * 2.5;
end


function cl = setup_solid(cl)
    for i = 1:length(cl)
        cl(i).Z = ones(1, size(cl(i).XYZmm, 2));
    end
end




% function CLU = clusters2CLU(clusters, [opt] M)
%
% Inputting an M matrix will transform the coordinates
% by that M, to convert between voxel sizes, etc.
%
% by Tor Wager
function CLU = clusters2CLU(clusters, varargin)

    if ~isfield(clusters(1), 'threshold'), clusters(1).threshold = 1; end

    CLU.XYZmm = cat(2, clusters.XYZmm);
    CLU.XYZ = cat(2, clusters.XYZ);
    try
        CLU.Z = cat(2, clusters.Z);
    catch
        %CLU.Z = cat(1, clusters.Z)';
        for i = 1:length(clusters)
            if size(clusters(i).Z, 1) > size(clusters(i).Z, 2)
                clusters(i).Z = clusters(i).Z';
            end
        end
        CLU.Z = cat(2, clusters.Z);
    end
    CLU.title = clusters(1).title;

    CLU.u = clusters(1).threshold;
    CLU.threshold = clusters(1).threshold;

    CLU.voxSize = clusters(1).voxSize;
    CLU.VOX = clusters(1).voxSize;

    if nargin > 1
        CLU.M = varargin{1};
        CLU = transform_coordinates(CLU, CLU.M);
    else
        switch class(clusters)
            case 'struct'
                has_M_field = isfield(clusters, 'M');
                
            case 'region'
                has_M_field = ~isempty(strmatch('M', fieldnames(clusters)));
                
            otherwise
                error('cluster_orthviews: cl input is wrong class...');
        end
        
        if has_M_field
            CLU.M = clusters(1).M;
        else
            error('invalid clusters variable: clusters must have defined M field/attribute.');
        end
    end
    
    CLU.numVox = size(CLU.XYZmm, 2);

    if isempty(CLU.Z), 
        disp('clusters Z field is empty.  Filling with ones as a placeholder.')
        CLU.Z = ones(1, size(CLU.XYZmm, 2));
    end

    if size(CLU.Z, 1) > 1, CLU.allZ = CLU.Z; CLU.Z = CLU.Z(1, :); end
end


% function [clusters, SPM, xX, xCon] = tor_extract_rois(imnames [can be empty], [opt] SPM, [opt] VOL, [opt] xX)
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
function [clusters, SPM, xX, xCon] = tor_extract_rois(imnames, varargin)

    verbose = 0;
    clusters = []; clustersin = [];

    if nargin == 1
        % ----------------------------------------------------------------------------------
        % get SPM info from SPM.mat
        % ----------------------------------------------------------------------------------
        [SPM, VOL, xX, xCon, xSDM] = spm_getSPM;
    elseif nargin == 2
        clustersin = varargin{1};

        if ~isempty(imnames), 
            % ----------------------------------------------------------------------------------
            % check mat files for compatibility
            % -----------------------------------------------------------------
            if isfield(clustersin, 'M')
                V = spm_vol(imnames(1, :));
                [chk, clustersin] = check_spm_mat(clustersin(1).M, V(1).mat, clustersin);
                if chk, disp('Warning!  Mat files for clusters and images do not match; using XYZmm to adjust cluster voxel coordinates'), end
            end
        end

        allxyz = cat(2, clustersin.XYZ);
    elseif nargin == 3
        SPM = varargin{1};
        VOL = varargin{2};
        if isempty(SPM.XYZ), disp('tor_extract_rois: no voxels to extract'), return, end
        allxyz = SPM.XYZ;
    elseif nargin == 4
        SPM = varargin{1};
        VOL = varargin{2};
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
            cl_index = spm_clusters(SPM.XYZ);
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
    % load image files, if possible
    % ----------------------------------------------------------------------------------
    if ~isempty(imnames)
        if size(imnames, 1) < 100
            if verbose, fprintf(1, '\nReading %3.0f images...', size(imnames, 1)); end
            % Image loading a la SPM, and manual extraction.
            V = spm_vol(imnames);
            vols = spm_read_vols(V);
        end
    end



    % ----------------------------------------------------------------------------------
    % get the data from ALL voxels; all_data stores raw data for all
    % clusters
    % ----------------------------------------------------------------------------------
    O.coords = allxyz'; clear allxyz
    if ~isempty(imnames), 

        if size(imnames, 1) < 100
            % we have already loaded images
            ts = timeseries4(O.coords, vols);
            all_data = ts.indiv;
        else    % load the images and extract data
            ts = timeseries4(O.coords, imnames);
            all_data = ts.indiv;

            % timeseries2 does not work on Mac OSX
            % timeseries2 - maybe slower?, but more memory efficient for large n of images
            % does not support multiple datatypes w.i timeseries - e.g., for masked 1st subject
            %try
            %    ts = timeseries2('multi', imnames, O);
            %    cl.timeseries = ts.avg;
        end
    end


    % ----------------------------------------------------------------------------------
    % define each cluster as cell in array.
    % ----------------------------------------------------------------------------------

    for i = 1:max(cl_index)
        if verbose && i == 1
            fprintf(1, '\n\t%3.0f clusters: Extracting %03d ', max(cl_index), i), 
        elseif verbose
            fprintf(1, '\b\b\b%03d', i)
            if i == max(cl_index)
                fprintf(1, '\n')
            end
        end

        % voxels in this cluster
        a = find(cl_index == i);


        % make cluster, if we don't have it
        if isempty(clustersin)
            cl.title = SPM.title;
            cl.threshold = SPM.u;
            cl.voxSize = VOL.VOX;
            cl.M = VOL.M;
            cl.name = [cl.title '_' num2str(i) '_' mat2str(size(a, 2)) '_voxels'];
            cl.numVox = size(a, 2);
            cl.Z = SPM.Z(a);
            cl.XYZmm = SPM.XYZmm(:, a);
            cl.XYZ = SPM.XYZ(:, a);

            try
                cl.pVoxelLev = spm_P(1, 0, max(cl.Z), SPM.df, SPM.STAT, VOL.R, SPM.n);
                cl.pClustLev = spm_P(1, cl.numVox/prod(VOL.FWHM), SPM.u, SPM.df, SPM.STAT, VOL.R, SPM.n);
            catch
                % warning('Can''t get SPM p voxel and cluster values.  Skipping spm_P.')
            end

            if ~isempty(cl.Z), 
                % report number of sub-cluster peaks within cluster
                [N] = spm_max(cl.Z, cl.XYZ);
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
                cl.all_data = all_data(:, a);
            else
                % clusters input
                nv = size(clustersin(i).XYZ, 2);
                cl.all_data = all_data(:, 1:nv);         % take first n data vectors
                all_data(:, 1:nv) = [];                 % remove from all_data
            end


            cl.timeseries = nanmean(cl.all_data', 1)';

            if size(imnames, 1) < 200 && i == 1
                cl.imnames = imnames;
            elseif i == 1
                cl.imnames = imnames([1 end], :);
            else
                cl.imnames = [];
            end


            % ----------------------------------------------------------------------------------
            % get the SNR for the cluster if it looks like rfx data rather than individual data
            % ----------------------------------------------------------------------------------
            if size(imnames, 1) < 60
                try
                    cl.snr_avgts = get_snr(cl.timeseries);
                    cl.snr = get_snr(cl.all_data);
                    if verbose, fprintf(1, ' : avg SNR = %3.2f, range %3.2f - %3.2f ', cl.snr_avgts, min(cl.snr), max(cl.snr)), end

                    % calculate # subjects w/ estimates above zero
                    cl.numpos = sum(cl.timeseries > 0);

                    % calc 80% power level, no mult. comp corr
                    % from Jerry Dallal @ Tufts U.  http://www.tufts.edu/~gdallal/SIZE.HTM
                    % (16 * sd^2 / est^2) + 1
                    cl.power80 = (4 ./ cl.snr_avgts)^2 + 1;
                catch
                    warning('Problem with SNR calculation.')
                    keyboard
                end
            end

        else
            % disp('No timeseries data extracted - image names empty.')
        end

        cl.center = mean(cl.XYZ', 1);

        % in case we've added a full matrix of Z-scores with
        % cluster_barplot, etc.
        if min(size(cl.Z)) > 1, cl.Z = ones(1, size(cl.XYZ, 2));  end

        if size(cl.XYZmm, 2) > 1
            cl.mm_center = center_of_mass(cl.XYZmm, cl.Z);
        else
            cl.mm_center = cl.XYZmm';
        end
        clusters = [clusters, cl];
    end




    if ~isempty(imnames) && exist('xX') == 1

        % ----------------------------------------------------------------------------------
        % adjust timeseries for each cluster and fit xX.X model to timeseries data
        % ----------------------------------------------------------------------------------
        try
            clusters = analyze_cluster_rois(clusters, xX);
        catch
            disp('Error analyzing timeseries clusters - skipping analysis.')
        end
    end

    % ----------------------------------------------------------------------------------
    % save this data in current directory
    % ----------------------------------------------------------------------------------
    %matname = ['clusters_' deblank(SPM.title(~(SPM.title==' ')))];
    %str = ['save ' matname ' clusters'];
    %eval(str)
end



% function [ts, vols, chunksize] = timeseries4(coords, P, [chunksize], [nochk])
%
% Simple extraction from images named in str mtx P or vols
% from voxel coordinates (not mm!) listed in coords
%
% P can be filenames in str matrix (char array)
% or 4-D array of all volume info (vols)
% (i.e., put vols in output back in as P)
%
% ts is timeseries, with fields avg and indiv for average cluster
% and individual voxels
%
% Loads images 'chunksize' at a time; default is based on memory
% size of 2^29 bytes; empty uses default
% Optional 4th argument suppresses data validity checking
%
% vols is 4-D array of all data, [x y z time(image)]
%
% Uses spm_vol and spm_read_vols and spm_slice_vol
%
% Tor Wager, 2/9/05     change from timeseries3: uses slice-by-slice
%                       method, faster.
% Tor Wager, 12/11/05   speedup for getdata with large n. voxels; cosmetic
%                       changes to output of volume method.

function [ts, vols, chunksize] = timeseries4(coords, P, varargin)

    % -------------------------------------------------------------------
    % * set up input arguments
    % -------------------------------------------------------------------
    maxmem = 2^28;      % note:tested the G5s up to 200 imgs, no slowdown, so no need to chunk...
    chunksize = [];

    global defaults
    if isempty(defaults), spm_defaults, end

    if ischar(P), 
        fprintf(1, 'Map vols: '); t1= clock;
        V = spm_vol(P);
        nimages = length(V);   % number of images in data
        fprintf(1, '%3.0f s.\n', etime(clock, t1));
    elseif ismatrix(P), 
        chunksize = NaN;
        nimages = size(P, 4);    % number of images in data
    else
        error('P input must be string or data matrix.');
    end

    if length(varargin) > 0, 
        chunksize = varargin{1};
    end


    if size(coords, 2) ~= 3, coords = coords';  end


    % -------------------------------------------------------------------
    % * get extraction method
    % -------------------------------------------------------------------
    whslices = unique(coords(:, 3));     % which slices to extract
    nslices = length(whslices);         % number of z slices to extract from

    extract_type = 'slice';
    if ~ischar(P), extract_type = 'volume'; end
    if nslices > 15, extract_type = 'volume'; end

    % -------------------------------------------------------------------
    % * read images
    % -------------------------------------------------------------------

    switch extract_type
        % * slice loading method
        % -------------------------------------------------------------------
        case 'slice'

            ts.indiv = NaN .* zeros(nimages, size(coords, 1)); % placeholder for data extracted

            for i = 1:nslices

                sliceno = whslices(i);               % slice number


                whcoords = find(coords(:, 3) == sliceno);    % indices of in-slice coordinates
                slcoords = coords(whcoords, :);              % coordinates to extract data from for this slice

                sl = timeseries_extract_slice(V, sliceno);

                for c = 1:size(slcoords, 1)

                    ts.indiv(:, whcoords(c)) = sl(slcoords(c, 1), slcoords(c, 2), :);

                end

                ts.avg = nanmean(ts.indiv')';
                vols = sl;
            end


            % * whole-brain loading method
            % -------------------------------------------------------------------
        case 'volume'
            if isempty(chunksize)
                v = spm_read_vols(V(1));

                tmp = whos('v'); imgsize = tmp.bytes;
                chunksize =   floor(maxmem ./ imgsize);  %round(maxmem ./ (size(coords, 1)^3 .* 16));   % images to load at a time
            end

            if ischar(P)
                if length(V) < chunksize
                    fprintf(1, '\tChunking %3.0f into %3.0f imgs: ', length(V), chunksize);
                    t1 = clock;
                    fprintf(1, 'Load: ');
                    vols = spm_read_vols(V);
                    fprintf(1, '%3.0f s. Cluster. ', etime(clock, t1));
                    t1 = clock;
                    [ts.indiv, ts.avg] = getdata(vols, coords, maxmem);
                    fprintf(1, '%3.0f s. ', etime(clock, t1));
                else
                    % chunk it!  load chunksize images as a whole, extract, and
                    % concatenate
                    ts.indiv = []; ts.avg = [];
                    ind = 1;
                    for i = 1:chunksize:length(V)
                        t1 = clock;
                        fprintf(1, 'Load: %3.0f ', i);
                        e = min(i+chunksize-1, length(V));
                        wh = i:e;
                        vols = spm_read_vols(V(wh));
                        fprintf(1, '%3.0f s. Cluster. ', etime(clock, t1));
                        t1 = clock;
                        [indiv{ind}, avg{ind}] = getdata(vols, coords, maxmem);
                        ind = ind + 1;
                        fprintf(1, '%3.0f s. ', etime(clock, t1));
                    end
                    fprintf(1, 'Cat. ')
                    ts.indiv = cat(1, indiv{:});
                    ts.avg = cat(1, avg{:});
                    clear indiv; clear avg;
                end
            else
                % volumes already loaded
                vols = P; P = 1;
                [ts.indiv, ts.avg] = getdata(vols, coords, maxmem);
            end
    end         %   end switch extraction type


    % -------------------------------------------------------------------
    % * check for proper extraction against spm_read_vols
    % -------------------------------------------------------------------

    if length(varargin) > 1 || ~ischar(P)
        % already loaded or suppress checking.
    else
        fprintf(1, 'Chk.\n')
        chk = check_timeseries_vals(V, ts.indiv, coords);
        if chk, keyboard, end
    end
end





function [ind, avg] = getdata(vols, coords, maxmem)
    if size(coords, 1) == 1  % only one voxel
        co = 1;
        ind(:, co) = squeeze(vols(coords(co, 1), coords(co, 2), coords(co, 3), :));
    elseif size(coords, 1)^3 < inf  % always do this.
        % time increases linearly with size of matrix; so do it in chunks.
        csz = round(sqrt(size(coords, 1)));  % optimal chunk size to keep arrays as small as possible.
        indx = 1;
        for i = 1:csz:size(coords, 1)
            tmp = [];
            for co = i:min(i+csz-1, size(coords, 1))
                %t1 = clock;
                tmp = [tmp squeeze(vols(coords(co, 1), coords(co, 2), coords(co, 3), :))];
                %et(co) = etime(clock, t1);
            end
            ind{indx} = tmp;
            indx = indx + 1;
            %ind = [ind squeeze(vols(coords(co, 1), coords(co, 2), coords(co, 3), :))];
        end
        ind = cat(2, ind{:});
    else    % not a good idea speed-wise, apparently.
        tmp =  vols(coords(:, 1), coords(:, 2), coords(:, 3), :); % the values of interest are on the 3-D diagonals of tmp
        s = size(coords, 1);      % we can get the index values for diagonals by skipping elements of sv, below
        i = 1:s.^2 + s + 1:s^3;  % same as t1 = [1:size(coords, 1)]'; i = sub2ind(size(tmp), t1, t1, t1)

        sv = (1:s^3:prod(size(tmp))) - 1; % starting values for each volume (minus one, so we add this to i)
        sv = repmat(sv, s, 1) + repmat(i', 1, size(sv, 2));  % get the matrix of index values for each voxel at each time
        ind = tmp(sv)';
    end

    if size(ind, 2) > 1, 
        avg = nanmean(ind')';
    else
        avg = ind;
    end
end




function sl = timeseries_extract_slice(V, sliceno)
    % function sl = timeseries_extract_slice(V, sliceno)
    %
    % For a given set of image names or memory mapped volumes (V)
    % extracts data from slice # sliceno and returns an X x Y x time
    % matrix of data.
    % uses spm_slice_vol.m

    if ischar(V), V = spm_vol(V); end

    mat = spm_matrix([0 0 sliceno]);     % matrix for spm_slice_vol

    for i = 1:length(V)
        sl(:, :, i) = spm_slice_vol(V(i), mat, V(i).dim(1:2), 0);
    end
end




function chk = check_timeseries_vals(V, dat, coords)
    % function chk = check_timeseries_vals(V, dat, coords)
    %
    % Checks a random subset of up to 5 images against extracted data values
    % to make sure the right data is being extracted.
    %
    % V is memory-mapped volumes (see spm_vol)
    % dat is extracted data
    % coords is voxel coordinates, n rows x 3 columns
    %
    % tor wager

    n = min(5, length(V));

    % get random set of n images, up to 5
    wh = randperm(length(V));
    wh = wh(1:n);

    % select these rows in V and dat
    % dat is [images, coordinates]
    V = V(wh);
    dat = dat(wh, :);

    % get random set of nn coordinates, up to 5

    nc = size(coords, 1);
    nn = min(5, nc);
    whc = randperm(nc);
    whc = whc(1:nn);
    coords = coords(whc, :);

    % select these columns in dat
    dat = dat(:, whc);

    v = spm_read_vols(V);
    for i = 1:nn    % for each coordinate
        dat2(:, i) = squeeze(v(coords(i, 1), coords(i, 2), coords(i, 3), :));
    end

    chk = dat - dat2;
    chk = any(chk(:));

    if chk, warning('Problem with timeseries3!! Extracted data do not match expected values.  Quitting at error'); end
end



function snr = get_snr(data)
    % snr = get_snr(data)
    %
    % data is a matrix whos columns index voxels, and rows index subjects (or trials, etc.)
    %
    % Tor Wager

    mystd = nanstd(data);
    mystd(mystd == 0) = NaN;
    snr = nanmean(data) ./ mystd;
end


function doinitwhitebg
    try
        global st
        st.vols{1}.black2white = 1;
        bwexist = strfind(st.plugins, 'black2white');
        bwexist = any(cat(2, bwexist{:}));
        if ~bwexist
            st.plugins{end+1} = 'black2white';
        end
        fh = findobj('Type', 'Figure', 'Tag', 'Graphics'); % spm fig
        if exist('fh') && ishandle(fh), set(fh, 'Color', 'w'); end
        
    catch
        disp('Error in development version of initwhitebg!!')
    end
end


% ISMATRIX: Returns 1 if the input matrix is 2+ dimensional, 0 if it is a scalar 
%           or vector.
%
%     Usage ismat = ismatrix(X)
%
% RE Strauss, 5/19/00

function ismat = ismatrix(X)
  [r,c] = size(X);
  if (r>1 && c>1)
    ismat = 1;
  else
    ismat = 0;
  end

end
