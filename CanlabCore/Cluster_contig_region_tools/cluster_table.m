function clusters = cluster_table(clusters, varargin)
% Print output of clusters in table
%
% Option to print text labels from Carmack atlas.
%
% Database loading is done from talairach_info.mat which should be in the 
% path.
%
% To speed up performance, load talairach_info.mat in the base workspace or
% calling function and include xyz, L3 and L5 as inputs to cluster_table.
%
% :Example:
% ::
%
%    % create subclusters on the fly, prompt for labels
%    cluster_table(cl);
%
%    % no subclusters, no labels
%    cluster_table(cl, 0, 0);
%
%    % do subclusters, no labels
%    cluster_table(cl, 1, 0);
%
%    create subclusters on the fly, do labels
%    cluster_table(cl, 1, 1);
%
%    % 3 input variables following 'tal_info' are interpreted as xyz, L3,
%    % and L5 from talairach_info.mat.
%    cluster_table(..., 'tal_info', xyz, L3, L5);
%
%    % loads labels from taldata.mat (Talairach database) instead of
%    % talairach_info.mat. Note that you should use the 'tal_info' call
%    % above if xyz, L3, and L5 have already been loaded to theworkspace
%    % from taldata.mat. Also, if the talairach database is being used,
%    % your cl.XYZmm values MUST correspond to the TALAIRACH, NOT MNI,
%    % database, or the labels will be innaccurate.
%    cluster_table(..., 'talairach');
%
%    % print table to ASCII file, 'filename', instead of to the matlab
%    % command window.
%    cluster_table(..., 'writefile','filename');
%
%    % any set of inputs from above, also print clusters.myfield in output
%    cluster_table(..., 'myfield');
%
% ..
%    Tor Wager
%
%    List of edits can be found in programmers' notes (edit this function)
%
%    Programmers' notes
%
%    Edited 10/17/06 by jared to obviate multiple calls to Carmack_get_label
%    should reduce the time it takes to get labels for clusters by 50%.
%    Also edited to allow loading of the Talairach, instead of Carmack,
%    database. 
%
%    Edited 9/07/06 by jared to fix incorrect behavior when using 'tal_info'
%    keyword. Added 'writefile' keyword to allow writing to file instead of
%    the command window.
%
%    Edited 7/13/06 by jared to obviate the use of global variables
%
%    Edited 7/8/06 by tor to handle Z values that are not row vectors
%    previously returned incorrect max stats if Z is column vector
%
%    Edited 7/8/06 to add subclusters from spm_max functionality
%    this is now the default

verbose = 1;
subc = 1;
if isempty(clusters), disp('No results to print.'); return, end
if isfield(clusters, 'M'), M = clusters(1).M; end
if length(varargin) > 0, subc = varargin{1}; end
if length(varargin) > 1, dolabs = varargin{2}; end

if isa(clusters, 'region')
    warning off
    for i = 1:length(clusters)
        cl(i) = struct(clusters(i));
    end
    warning on
    clusters = cl;
end

% -----------------------------------------------------------------------------------
% Set up user-input fields
% -----------------------------------------------------------------------------------

printfields = {};
i=3;
while i <= length(varargin)
    if strcmp(varargin{i}, 'tal_info')
        xyz=varargin{i+1};L3=varargin{i+2};L5=varargin{i+3};
        i=i+3;
    elseif strcmp(varargin{i},'talairach')
        fprintf(1, 'Loading database.');
        load taldata
    elseif strcmp(varargin{i}, 'writefile')
        writefile=1;outfilename=varargin{i+1};
        i=i+1;
    elseif strcmp(varargin{i}, 'noverbose')
        verbose = 0;
    else
        printfields{end+1} = varargin{i};
    end
    i=i+1;
end

if ~exist('writefile','var'),writefile=0;end

% -----------------------------------------------------------------------------------
% Set up subclustering option
% -----------------------------------------------------------------------------------
dosubcluster = 1;   % default

if isstruct(subc)
    % do nothing
elseif subc == 0
    % turn option off
    dosubcluster = 0;
else
    % create
    if verbose, fprintf(1, 'Getting local maxima within 10 mm for subcluster reporting\n'); end
    try
        subc = subclusters_from_local_max(clusters, 10);
    catch
       disp('**************************************************************')
       disp('Error in  subclusters_from_local_max')
       disp('This could have multiple causes.  One cause may be that a continguous region in a cluster')
       disp('Is too big for matlab''s memory to calculate distances among all voxels with pdist.');
       disp('If you get this error with small clusters, you may be missing subfunctions or have other related issues.')
        
       subc = clusters;
       dosubcluster = 0;
    
    end
end


% -----------------------------------------------------------------------------------
% Set up labeling output (see bottom of function for old BA code!!!)
% -----------------------------------------------------------------------------------



% new for Carmack labels, do L3 and L5 (most informative levels).
if ~exist('dolabs', 'var')
    dolabs = input('Do text labels for clusters (current = Carmack labels)?');
end

if dolabs && ~(exist('talairach_info.mat') == 2)
    disp('Cannot find talairach_info.mat, so cannot get Talairach labels.');
    dolabs = 0;
end

if dolabs
    global xyz L3 L5
    fprintf(1, 'Loading database.');
    if ~exist('xyz', 'var') || isempty(xyz), load talairach_info xyz, end
    if ~exist('L3', 'var') || isempty(L3), load talairach_info L3, end
    if ~exist('L5', 'var') || isempty(L5), load talairach_info L5, end

    fprintf(1, 'Done. Getting text labels.');
    fprintf(1, '%03d', 0);
    for i = 1:length(clusters)
        fprintf(1, '\b\b\b%03d', i);
        [name, perc, number, totalnum, stricbm{i}] = Carmack_get_label(clusters(i).XYZmm, {L3,L5}, xyz);
%         [name, perc, number, totalnum, stricbm2{i}] = Carmack_get_label(clusters(i).XYZmm, L5, xyz);
    end
    if dosubcluster
        fprintf(1, '\nGetting text labels for subclusters.');
        fprintf(1, '%03d', 0);
        for i = 1:length(subc)
            fprintf(1, '\b\b\b%03d', i);
            [name, perc, number, totalnum, stricbm2{i}] = Carmack_get_label(subc(i).XYZmm, {L3,L5}, xyz);
        end
    end
    fprintf(1, 'Done.\n');
end



% -----------------------------------------------------------------------------------
% Start print-out; loop through clusters
% -----------------------------------------------------------------------------------

for i = 1:length(clusters)

    if isfield(clusters, 'correl'), 
        if ~isempty(clusters(i).correl)
            try
                cmatx(i, 1) = clusters(i).correl;
            catch
                cmatx(i, 1) = NaN;
            end
        else
            cmatx(i, 1) = NaN;
        end
    else cmatx(i, 1) = NaN;
    end

    cmatx(i, 2) = clusters(i).numVox;

    cmatx(i, 3) = get_maxstat(clusters(i));


    % -----------------------------------------------------------------------------------
    % Define variables to report in an easy-to-access matrix (cmatx)
    % -----------------------------------------------------------------------------------

    if isfield(clusters, 'corr_range'), 
        cmatx(i, 8) = clusters(i).corr_range(1);
        cmatx(i, 9) = clusters(i).corr_range(end);
    end

    if isfield(clusters, 'snr_avgts'), 
        cmatx(i, 10) = clusters(i).snr_avgts;
    end

    if isfield(clusters, 'snr'), 
        cmatx(i, 11) = min(clusters(i).snr);
        cmatx(i, 12) = max(clusters(i).snr);
    end

    if isfield(clusters, 'numpos') && isfield(clusters, 'power80'), 
        cmatx(i, 13) = clusters(i).numpos;
        cmatx(i, 14) = ceil(clusters(i).power80);
    end

    if isfield(clusters, 'numpeaks'), cmatx(i, 15) = clusters(i).numpeaks; end

    x(i) = clusters(i).mm_center(1);
    y(i) = clusters(i).mm_center(2);
    z(i) = clusters(i).mm_center(3);

    if exist('v') == 1
        vox = mm2voxel(clusters(i).mm_center, V);
        cmatx(i, 16) = round(v(vox(1), vox(2), vox(3)));
    end
end

% -----------------------------------------------------------------------------------
% Print table header
% -----------------------------------------------------------------------------------

disp(' ')
if isempty(clusters), disp('No clusters.'), return, end
if isfield(clusters, 'name'), disp(clusters(1).name), end
if isfield(clusters, 'descrip1'), disp(clusters(1).descrip1), end
if isfield(clusters, 'descrip2'), disp(clusters(1).descrip2), end

if writefile
    fid=fopen(outfilename,'w');
else
    fid=1;
end

if isfield(clusters, 'Z_descrip') && verbose
    fprintf('Z field contains: %s (shown in maxstat)\n', clusters(1).Z_descrip);
elseif verbose
    fprintf('Z field contains: Unknown statistic (shown in maxstat).\n');
end

if isfield(clusters, 'center') && exist('M') == 1 && isfield(clusters, 'from_cluster')
    % sort by which cluster its from
    try cmatx = sortrows(cmatx, 4);catch end
    fprintf(fid, 'corr\tvoxels\tmaxstat\tfrom_clust\tmax_coords\n');
    for i = 1:size(cmatx, 1)
        fprintf(fid, '%3.2f\t%3.0f\t%3.2f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t\n', cmatx(i, 1), cmatx(i, 2), cmatx(i, 3), cmatx(i, 4), cmatx(i, 5), cmatx(i, 6), cmatx(i, 7));
    end
else
    disp(' ')
    if isfield(clusters, 'shorttitle'), fprintf(fid, 'Name\t'); end
    fprintf(fid, 'index\tx\ty\tz\tcorr\tvoxels\tvolume_mm3\tmaxstat\t');
    if isfield(clusters, 'numpeaks'), fprintf(fid, 'numpeaks\t'); end
    if isfield(clusters, 'corr_range'), fprintf(fid, 'mincorr\tmaxcorr\t'); end
    if isfield(clusters, 'snr_avgts'), fprintf(fid, 'snr_avgts(d)\t'); end
    if isfield(clusters, 'snr'), fprintf(fid, 'minsnr\tmaxsnr\t'); end
    if isfield(clusters, 'numpos') && isfield(clusters, 'power80'), 
        fprintf(fid, 'numpos\tpower80\t');
    end
    if exist('v') == 1
        fprintf(fid, ('BA\tBA_composition\t'));
    end
    if exist('stricbm') == 1
        %fprintf(fid, ('ICBM_single_subj\t'))
        fprintf(fid, ('Tal_Labels\t'));
    end
% 
%     if exist('stricbm2') == 1
%         %fprintf(fid, ('ICBM_single_subj\t'))
%         fprintf(fid, ('Carmack_Level5\t'));
%     end

    % print additional fields
    for i = 1:length(printfields)
        fprintf(fid, '%s\t', printfields{i});
    end

    fprintf(fid, '\n');

    % -----------------------------------------------------------------------------------
    % Print a row for each cluster
    % -----------------------------------------------------------------------------------

    for i = 1:size(cmatx, 1)
        if isfield(clusters, 'shorttitle'), fprintf(fid, '%s\t', clusters(i).shorttitle);end
        fprintf(fid, '%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.2f\t%3.0f\t%3.0f\t%3.2f\t', i, x(i), y(i), z(i), cmatx(i, 1), cmatx(i, 2), cmatx(i, 2).*prod(clusters(i).voxSize), cmatx(i, 3));
        if isfield(clusters, 'numpeaks'), fprintf(fid, '%3.0f\t', cmatx(i, 15)); end
        if isfield(clusters, 'corr_range'), fprintf(fid, '%3.2f\t%3.2f\t', cmatx(i, 8), cmatx(i, 9)); end
        if isfield(clusters, 'snr_avgts'), fprintf(fid, '%3.2f\t', cmatx(i, 10)); end
        if isfield(clusters, 'snr'), fprintf(fid, '%3.2f\t%3.2f\t', cmatx(i, 11), cmatx(i, 12)); end
        if isfield(clusters, 'numpos') && isfield(clusters, 'power80'), 
            fprintf(fid, '%3.0f\t%3.0f\t', cmatx(i, 13), cmatx(i, 14));
        end
        if exist('v', 'var') == 1
            fprintf(fid, ('%3.0f\t'), cmatx(i, 16));
        end

        if exist('strs', 'var') == 1
            fprintf(fid, '%s\t', strs);
        end

        if exist('stricbm', 'var') == 1
            fprintf(fid, ('%s\t'), stricbm{i});
        end

%         if exist('stricbm2', 'var') == 1
%             fprintf(fid, ('%s\t'), stricbm2{i});
%         end
        
        % print additional fields
        for j = 1:length(printfields)
            if ~isfield(clusters, printfields{j})
                warning([printfields{j} ' is not a field in clusters.']);
            else
                myval = clusters(i).(printfields{j});
                if isempty(myval), myval = NaN; end
                
                % get formatting string based on content
                fmtstring = get_formatstring(myval);
                
                fprintf(fid, fmtstring, myval);
               
                
            end
        end

        fprintf(fid, '\n');

        if dosubcluster
            % print sub-cluster table
            whsc = cat(1, subc.from_cluster) == i;
            whsc = find(whsc)';
            if length(whsc)>1
                for j = whsc
                    if ~dolabs
                        print_row(subc(j), j, clusters, fid)
                    else
                        print_row(subc(j), j, clusters, fid,stricbm2)
                    end
                end
            end
        end

    end
end

if writefile,fclose(fid);end

return





function print_row(clusters, i, bigcl, fid, varargin)
% prints a row for a subcluster (varargin{1}) below its corresponding cluster

if ~isempty(varargin)
    stricbm2=varargin{1};
end

if isfield(clusters, 'correl'), cmatx(i, 1) = clusters(1).correl;
else cmatx(i, 1) = NaN;
end

cmatx(i, 2) = clusters(1).numVox;
cmatx(i, 3) = get_maxstat(clusters(1));

if isfield(clusters, 'corr_range'), 
    cmatx(i, 8) = clusters(1).corr_range(1);
    cmatx(i, 9) = clusters(1).corr_range(end);
elseif isfield(bigcl, 'corr_range'), 
    cmatx(i, 8) = NaN;
    cmatx(i, 9) = NaN;
end

if isfield(clusters, 'snr_avgts'), 
    cmatx(i, 10) = clusters(1).snr_avgts;
elseif isfield(bigcl, 'snr_avgts'), 
    cmatx(i, 10) = NaN;
    cmatx(i, 10) = NaN;
end

if isfield(clusters, 'snr'), 
    cmatx(i, 11) = min(clusters(1).snr);
    cmatx(i, 12) = max(clusters(1).snr);
elseif isfield(bigcl, 'snr'), 
    cmatx(i, 11) = NaN;
    cmatx(i, 12) = NaN;
end

if isfield(clusters, 'numpos') && isfield(clusters, 'power80') 
    cmatx(i, 13) = clusters(1).numpos;
    cmatx(i, 14) = ceil(clusters(1).power80);
elseif isfield(bigcl, 'numpos') && isfield(bigcl, 'power80') 
    cmatx(i, 13) = NaN;
    cmatx(i, 14) = NaN;
end

%if isfield(clusters, 'numpeaks'), 
%    cmatx(i, 15) = clusters(1).numpeaks;
%elseif isfield(bigcl, 'numpeaks'), 
%    cmatx(i, 15) = NaN;
%end

fprintf(fid, '%s\t%3.0f\t%3.0f\t%3.0f\t%3.2f\t%3.0f\t%3.5f\t', ...
    '->', clusters(1).mm_center(1), clusters(1).mm_center(2), clusters(1).mm_center(3), ...
    cmatx(i, 1), cmatx(i, 2), cmatx(i, 3));

%if isfield(clusters, 'numpeaks'), fprintf(fid, '%3.0f\t', cmatx(i, 15)), end
fprintf(fid,'\t');   % skip numpeaks - makes no sense for subcluster
if isfield(bigcl, 'corr_range'), fprintf(fid, '%3.2f\t%3.2f\t', cmatx(i, 8), cmatx(i, 9)); end
if isfield(bigcl, 'snr_avgts'), fprintf(fid, '%3.2f\t', cmatx(i, 10)); end
if isfield(bigcl, 'snr'), fprintf(fid, '%3.2f\t%3.2f\t', cmatx(i, 11), cmatx(i, 12)); end
if isfield(bigcl, 'numpos') & isfield(bigcl, 'power80'), 
    fprintf(fid, '%3.0f\t%3.0f\t', cmatx(i, 13), cmatx(i, 14));
end

if exist('stricbm2', 'var') == 1
    fprintf(fid, ('%s\t'), stricbm2{i});
end

fprintf(fid, '\n');

return




function maxstat = get_maxstat(cl)
% max absolute value
if size(cl.Z, 2) ~= size(cl.XYZmm, 2), cl.Z = cl.Z'; end
[maxabs, whmax] = max(abs(cl.Z(1, :)));
maxstat = cl.Z(1, whmax(1));
return



function fmtstring = get_formatstring(myval)
if isempty(myval)
    fmtstring = '%3.0f\t';
elseif ischar(myval)
    fmtstring = '%s\t';
elseif isinteger(myval) || islogical(myval) || all( (myval - round(myval)) == 0 )
    fmtstring = repmat('%3.0f\t', 1, length(myval));
elseif myval < .01
    fmtstring = repmat('%3.4f\t', 1, length(myval));
else
    fmtstring = repmat('%3.2f\t', 1, length(myval));
end

return
                
% try to load table with BAs.  Brodmann Areas.
% Old: for ICBM and old Talairach atlas.
%
% if isfield(clusters, 'BAstring')
%     strs = str2mat(clusters.BAstring);
% else
% try
%     V = spm_vol(which('Tal_gray.img')); v = spm_read_vols(V);
%     V.M = V.mat;
%
%     % try to get BA composition for all clusters
%     %diary off
%     %disp('Finding composition of all clusters - ctrl c to cancel')
%     %[clusters, strs, clusvec, all_bas, ba_counts] = cluster_ba(clusters, 1:length(clusters));
%     %disp('Success - printing table.')
%     %diary on
%
% catch
% end
% end

% try to get ICBM composition for all clusters
% stricbm = icbm_orthview_label(clusters);
% for i = 1:length(clusters)
%     clusters(i).ICBMstr = stricbm{i};
% end
