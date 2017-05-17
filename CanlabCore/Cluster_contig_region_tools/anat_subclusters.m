function clout = anat_subclusters(cl,varargin)
% Clusters voxels within 'clusters' structure based on anatomical locations
% in space.  Outputs subgroups of smaller clusters.
%
% :Usage:
% ::
%
%    clout = anat_subclusters(cl,[resume at],[output cl to resume])
%
% ..
%    tor wager, 8/6/04
% ..

clout = []; startat = 1;

if length(varargin) > 0
    startat = varargin{1};
end
if length(varargin) > 1
    clout = varargin{2};
end

for i = startat:length(cl)
    
    if ~isfield(cl,'shorttitle'), cl(i).shorttitle = [];,end
    
    cluster_orthviews(cl(i),{[1 0 0]})
    
    fprintf(1,'Cluster %3.0f, with %3.0f voxels.',i,cl(i).numVox);
    if cl(i).numVox > 500, disp('Warning: may be slow! Sil plot not recommended'),end
    %D = squareform(pdist(cl(i).XYZ'));
    %[CLUST]=permnmds(D,100,min(5,round(size(D,2)./10)),3,1);           %FIND N by K best cluster/element solution
    % clas = CLUST.clustermat(CLUST.bestdims,:);
    dos = input('Do silhouette plot? (1/0) ');
    if dos, try, [s] = silhouetteplot(cl(i).XYZ');, catch, disp('Error making sil plot.'), end, end
    
    k = input('Choose number of subclusters: ');
    if max(size(cl(i).XYZ)) > 3, 
        %T = clusterdata(cl(i).XYZ','maxclust',k,'linkage','average');
        
        if k == 1
            T = ones(size(cl(i).XYZ,2),1)';
        else
            %point estimate cluster solution
            Y = pdist(cl(i).XYZ');
            Z = linkage_t(Y,'average');
            T = cluster(Z,k);
        end
    else
        % not enough data!
        T = ones(size(cl(i).XYZ,2),1)';
    end
    
    if size(T,1) > size(T,2), T = T';,end
    if size(cl(i).Z,1) > size(cl(i).Z,2), cl(i).Z = cl(i).Z';,end
             
    for j = 1:k
        
        if any(T == j)
        
            wh = find(T == j);
            
            if isempty(clout), 
                clout = cl(i);, 
            else, 
                clout(end+1) = cl(i);,

            end

            clout(end).XYZ = clout(end).XYZ(:,wh);
            clout(end).XYZmm = clout(end).XYZmm(:,wh);
            clout(end).Z = clout(end).Z(:,wh);
            clout(end).mm_center = mean(clout(end).XYZmm,2)';
            if isfield(clout,'all_data')
                clout(end).all_data = clout(end).all_data(:,wh);
            end
            if isfield(clout,'timeseries'),
                clout(end).timeseries = mean(clout(end).all_data,2);
            end
            
            cluster_orthviews(clout(end),{rand(1,3)},'add')
            clout(end).mm_center=mean(clout(end).XYZmm,2)';
            spm_orthviews('Reposition',clout(end).mm_center)
            
            clout(end).shorttitle = input('Enter short title for this cluster: ','s');
        end
        
    end
    

    disp('Saving clout_tmp for recovery.'), save clout_tmp clout
    
    % done - display output clusters
    %colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};
    %hh = findobj('Tag','Graphics');
    %h1 = get(hh,'Children');
    %f1 = figure('Color','w'); c= copyobj(h1,f1);

    %set(gcf,'Position', [ 210   322   929   407])    %[464 606 672./2 504./2])
    %for i = 1:length(c)
    %    tmp = get(c(i),'Position'); tmp(3) = tmp(3) .* .5;
    %    set(c(i),'Position',tmp);
    %end
    %keyboard
    
end

for ii = 1:length(clout), clout(ii).numVox = size(clout(ii).XYZmm,2);,end
    

return

