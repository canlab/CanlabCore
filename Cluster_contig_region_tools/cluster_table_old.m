function clusters = cluster_table(clusters,varargin)
% function cluster_table(clusters,[opt] subclusters)
% Print output of clusters in table
% Tor Wager
%
% Option to print text labels from Carmack atlas
% Database loading is done from talaraich_info.mat
% which should be in the path.  
% To speed up performance, declare global xyz, L3 and L5
% in the base workspace and re-use them in repeated calls
% to cluster_table.

global xyz L3 L5

if isfield(clusters,'M'),M = clusters(1).M;, end
if length(varargin) > 0, subc = varargin{1};, end

% try to load table with BAs.  Brodmann Areas.
% Old: for ICBM and old Talairach atlas.
%
% if isfield(clusters,'BAstring')
%     strs = str2mat(clusters.BAstring);
% else
% try 
%     V = spm_vol(which('Tal_gray.img')); v = spm_read_vols(V);
%     V.M = V.mat;
%     
%     % try to get BA composition for all clusters
%     %diary off
%     %disp('Finding composition of all clusters - ctrl c to cancel')
%     %[clusters,strs,clusvec,all_bas,ba_counts] = cluster_ba(clusters,1:length(clusters));
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

% new for Carmack labels, do L3 and L5 (most informative levels).
dolabs = 1; %dolabs = input('Do text labels for clusters (current = Carmack labels)?');

if dolabs && ~(exist('talairach_info.mat') == 2)
    disp('Cannot find talairach_info.mat, so cannot get Talairach labels.');
    dolabs = 0;
end

if dolabs
    fprintf(1,'Loading database.');
    if isempty(xyz), load talairach_info xyz, end
    if isempty(L3), load talairach_info L3, end
    if isempty(L5), load talairach_info L5, end

    fprintf(1,'Done. Getting text labels.');
    fprintf(1,'%03d',0);
    for i = 1:length(clusters)
        fprintf(1,'\b\b\b%03d',i);
        [name,perc,number,totalnum,stricbm{i}] = Carmack_get_label(clusters(i).XYZmm,L3,xyz);
        [name,perc,number,totalnum,stricbm2{i}] = Carmack_get_label(clusters(i).XYZmm,L5,xyz);
    end
    
    fprintf(1,'Done.\n');
end 
    
for i = 1:length(clusters)
    
    if isfield(clusters,'correl'), 
        if ~isempty(clusters(i).correl)
             try
                 cmatx(i,1) = clusters(i).correl;,
             catch
                 cmatx(i,1) = NaN;
             end
         else
             cmatx(i,1) = NaN;
         end
    else cmatx(i,1) = NaN;
    end

    cmatx(i,2) = clusters(i).numVox;
    maxZ = clusters(i).Z(abs(clusters(i).Z(1,:)) == max(abs(clusters(i).Z(1,:))));
    cmatx(i,3) = maxZ(1);
    
    % not working yet
    %clusters(i).cor_stat = tor_r2z(a(1,2),length(clusters(i).timeseries-1);
    
    % for use with tor_get_spheres2.m, breaking up cluster into spheres
    % OLD
    %if isfield(clusters,'center') & exist('M') == 1 & isfield(clusters,'from_cluster')
    %    cmatx(i,4) = clusters(i).from_cluster;
    %    centers{i} = (round(voxel2mm(clusters(i).XYZ(:,clusters(i).Z == max(clusters(i).Z)),M)'));
    %    cmatx(i,5) = centers{i}(1);
    %    cmatx(i,6) = centers{i}(2);
    %    cmatx(i,7) = centers{i}(3);
    %end
    
    % -----------------------------------------------------------------------------------
    % Define variables to report in an easy-to-access matrix (cmatx)
    % -----------------------------------------------------------------------------------
    
    if isfield(clusters,'corr_range'),
        cmatx(i,8) = clusters(i).corr_range(1);
        cmatx(i,9) = clusters(i).corr_range(end);
    end
    
    if isfield(clusters,'snr_avgts'),
        cmatx(i,10) = clusters(i).snr_avgts;
    end
    
    if isfield(clusters,'snr'),
        cmatx(i,11) = min(clusters(i).snr);
        cmatx(i,12) = max(clusters(i).snr);
    end  
        
    if isfield(clusters,'numpos') & isfield(clusters,'power80'),
        cmatx(i,13) = clusters(i).numpos;
        cmatx(i,14) = ceil(clusters(i).power80);
    end 
    
    if isfield(clusters,'numpeaks'), cmatx(i,15) = clusters(i).numpeaks;, end
        
    x(i) = clusters(i).mm_center(1);
    y(i) = clusters(i).mm_center(2);
    z(i) = clusters(i).mm_center(3);
    
    if exist('v') == 1
        vox = mm2voxel(clusters(i).mm_center,V);
        cmatx(i,16) = round(v(vox(1),vox(2),vox(3)));
    end
end

    % -----------------------------------------------------------------------------------
    % Print table header
    % -----------------------------------------------------------------------------------
    
    disp(' ')
    if isempty(clusters), disp('No clusters.'), return, end
    if isfield(clusters,'name'),disp(clusters(1).name),end
    
if isfield(clusters,'center') & exist('M') == 1 & isfield(clusters,'from_cluster')
    % sort by which cluster its from
    try,cmatx = sortrows(cmatx,4);,catch,end
    fprintf(1,'corr\tvoxels\tmaxZ\tfrom_clust\tmax_coords\n')
    for i = 1:size(cmatx,1)
        fprintf(1,'%3.2f\t%3.0f\t%3.2f\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t\n',cmatx(i,1),cmatx(i,2),cmatx(i,3),cmatx(i,4),cmatx(i,5),cmatx(i,6),cmatx(i,7))
    end
else
    disp(' ')
    if isfield(clusters,'shorttitle'),fprintf(1,'Name\t'),end
    fprintf(1,'index\tx\ty\tz\tcorr\tvoxels\tvolume_mm3\tmaxZ\t')
    if isfield(clusters,'numpeaks'),fprintf(1,'numpeaks\t'), end
    if isfield(clusters,'corr_range'), fprintf(1,'mincorr\tmaxcorr\t'), end
    if isfield(clusters,'snr_avgts'),fprintf(1,'snr_avgts(d)\t'), end
    if isfield(clusters,'snr'),fprintf(1,'minsnr\tmaxsnr\t'), end
    if isfield(clusters,'numpos') & isfield(clusters,'power80'),
        fprintf(1,'numpos\tpower80\t')
    end 
    if exist('v') == 1
        fprintf(1,('BA\tBA_composition\t'))
    end
    if exist('stricbm') == 1
        %fprintf(1,('ICBM_single_subj\t'))
        fprintf(1,('Carmack_Tal_Labels\t'))
    end
    
    if exist('stricbm2') == 1
        %fprintf(1,('ICBM_single_subj\t'))
        fprintf(1,('Carmack_Level5\t'))
    end
    fprintf(1,'\n')
    
    % -----------------------------------------------------------------------------------
    % Print a row for each cluster
    % -----------------------------------------------------------------------------------
    
    for i = 1:size(cmatx,1)
        if isfield(clusters,'shorttitle'),fprintf(1,'%s\t',clusters(i).shorttitle);,end
        fprintf(1,'%3.0f\t%3.0f\t%3.0f\t%3.0f\t%3.2f\t%3.0f\t%3.0f\t%3.5f\t',i,x(i),y(i),z(i),cmatx(i,1),cmatx(i,2),cmatx(i,2).*prod(clusters(i).voxSize),cmatx(i,3))
        if isfield(clusters,'numpeaks'),fprintf(1,'%3.0f\t',cmatx(i,15)), end
        if isfield(clusters,'corr_range'), fprintf(1,'%3.2f\t%3.2f\t',cmatx(i,8),cmatx(i,9)), end
        if isfield(clusters,'snr_avgts'), fprintf(1,'%3.2f\t',cmatx(i,10)), end
        if isfield(clusters,'snr'), fprintf(1,'%3.2f\t%3.2f\t',cmatx(i,11),cmatx(i,12)), end
        if isfield(clusters,'numpos') & isfield(clusters,'power80'),
            fprintf(1,'%3.0f\t%3.0f\t',cmatx(i,13),cmatx(i,14))
        end
        if exist('v') == 1
            fprintf(1,('%3.0f\t'),cmatx(i,16))
        end
        
        if exist('strs') == 1
            fprintf(1,'%s\t',strs(i,:))
        end
        
        if exist('stricbm') == 1
            fprintf(1,('%s\t'),stricbm{i})
        end
    
        if exist('stricbm2') == 1
            fprintf(1,('%s\t'),stricbm2{i})
        end
    
        fprintf(1,'\n')
        
        if length(varargin) > 0
            % print sub-cluster table
            whsc = cat(1,subc.from_cluster) == i;
            for j = find(whsc)'
                print_row(subc(j),j,clusters)
            end
        end
        
    end
end

return





function print_row(clusters,i,bigcl)
% prints a row for a subcluster (varargin{1}) below its corresponding cluster

    if isfield(clusters,'correl'), cmatx(i,1) = clusters(1).correl;,
    else cmatx(i,1) = NaN;
    end

    cmatx(i,2) = clusters(1).numVox;
    cmatx(i,3) = max(clusters(1).Z);
    
    if isfield(clusters,'corr_range'),
        cmatx(i,8) = clusters(1).corr_range(1);
        cmatx(i,9) = clusters(1).corr_range(end);
    elseif isfield(bigcl,'corr_range'),
        cmatx(i,8) = NaN;
        cmatx(i,9) = NaN;
    end
        
    if isfield(clusters,'snr_avgts'),
        cmatx(i,10) = clusters(1).snr_avgts;
    elseif isfield(bigcl,'snr_avgts'),
        cmatx(i,10) = NaN;
        cmatx(i,10) = NaN;
    end
    
    if isfield(clusters,'snr'),
        cmatx(i,11) = min(clusters(1).snr);
        cmatx(i,12) = max(clusters(1).snr);
    elseif isfield(bigcl,'snr'),
        cmatx(i,11) = NaN;
        cmatx(i,12) = NaN;
    end
        
    if isfield(clusters,'numpos') & isfield(clusters,'power80'),
        cmatx(i,13) = clusters(1).numpos;
        cmatx(i,14) = ceil(clusters(1).power80);
    elseif isfield(bigcl,'numpos') & isfield(bigcl,'power80'),
        cmatx(i,13) = NaN;
        cmatx(i,14) = NaN;
    end
    
    %if isfield(clusters,'numpeaks'), 
    %    cmatx(i,15) = clusters(1).numpeaks;, 
    %elseif isfield(bigcl,'numpeaks'),
    %    cmatx(i,15) = NaN;
    %end
    
fprintf(1,'%s\t%3.0f\t%3.0f\t%3.0f\t%3.2f\t%3.0f\t%3.5f\t', ...
    '->',clusters(1).mm_center(1),clusters(1).mm_center(2),clusters(1).mm_center(3), ...
    cmatx(i,1),cmatx(i,2),cmatx(i,3))

        %if isfield(clusters,'numpeaks'),fprintf(1,'%3.0f\t',cmatx(i,15)), end
        fprintf('\t')   % skip numpeaks - makes no sense for subcluster
        if isfield(bigcl,'corr_range'), fprintf(1,'%3.2f\t%3.2f\t',cmatx(i,8),cmatx(i,9)), end
        if isfield(bigcl,'snr_avgts'), fprintf(1,'%3.2f\t',cmatx(i,10)), end
        if isfield(bigcl,'snr'), fprintf(1,'%3.2f\t%3.2f\t',cmatx(i,11),cmatx(i,12)), end
        if isfield(bigcl,'numpos') & isfield(bigcl,'power80'),
            fprintf(1,'%3.0f\t%3.0f\t',cmatx(i,13),cmatx(i,14))
        end
        
        fprintf(1,'\n')
        
return