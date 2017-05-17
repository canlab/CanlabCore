function [clout] = clusters2roimask(cl)
% The purpose of this function is to facilitate making masks with ROIs for
% future studies, give a clusters structure.  ROIs are constrained to be
% within activation blobs specified by input clusters, and are masked by
% selected ICBM regions. Clusters may be smoothed before or after masking.
%
% :Usage:
% ::
%
%     [clout] = clusters2roimask(cl)

% Option to do 3 things, in order:
%   1) Enlarge selected clusters
%   2) Mask clusters with anatomical regions from ICBM atlas
%   3) Subdivide clusters using hierarchical clustering of voxel coordinates
%
% Output is a clusters file and a mask file in a 2 x 2 x 2 standard brain
% space.
%
% ..
%    Well, no mask file yet. And no shrinking.
%
%    Tor Wager, Aug 2004.
% ..

gom = input('Mask clusters with anatomical regions from ICBM (1/0)? ');
clout = [];

for i = 1:length(cl)
    
    icbm_localize(cl(i));
    [str,num] = icbm_orthview_label(cl(i));

    % -------------------------------------------------------
    % Cluster enlargement
    % -------------------------------------------------------
    
    go = 1;
    while go
        go = input('Enlarge (1) or shrink (-1) cluster, 0 if done? ');
        if go > 0
            cl(i) = enlarge_cluster(cl(i));
        elseif go < 0
            cl(i) = shrink(cl(i));
        else
        end
        [num,name] = icbm_localize(cl(i));    
        [str,num,name,perc] = icbm_orthview_label(cl(i));
    end
    
    if ~gom, clout = cl;, end
    
    % -------------------------------------------------------
    % Cluster masking
    % -------------------------------------------------------
    disp('Here is a list of anatomical regions for this cluster, indexed by number:')
    disp('You can combine anatomical regions by entering a vector of numbers: [3 60]')
    
    gom1 = gom; % save orig all clu
    clxyz = round(cl(i).XYZmm');

        
    while gom
        % ITERATE CHOICE of ANATOMICA REG
        % 
        
        % re-display list
        
        for jj = 1:length(num)
            fprintf('%3.0f\tPerc: %3.0f\t%s\t\n',num(jj),perc(jj),name{jj})
        end
        
        
        % get choice and make sure its valid
        
        stopme = 0;
        while stopme == 0
            stopme = 1;
            gom = input('Mask with: Enter anatomy number, vector to combine, or 0 if done: ');
            if gom ~= 0, 
                for jj = 1:length(gom), if ~any(gom(jj) == num), stopme = 0;, end, end
            end
        end
        
                
        if gom
            % Do the mask for THIS CLUSTER
            
            if isempty(clout), clout = cl(i);, else, clout(end+1) = cl(i);,end
            
            XYZmm = []; XYZ = [];
            for jj = 1:length(gom)
                [XYZmmtmp,XYZtmp] = icbm_reader('coordinates',gom(jj));
                XYZmm = [XYZmm XYZmmtmp];
                XYZ = [XYZ XYZtmp];
            end
            
            % lump together if multiple areas
            
            if iscell(XYZmm), tmp = cat(2,XYZmm{:});, XYZmm = unique(tmp','rows')';, end
            if iscell(XYZ), tmp = cat(2,XYZ{:});, XYZ = unique(tmp','rows')'; end
            

            %[xyztmp,b,c] = intersect(XYZmm',round(cl.XYZmm)','rows');
            [clout(end).XYZmm] = dominance_point_match(round(cl(i).XYZmm'),XYZmm',0)';
            clout(end).XYZ = mm2voxel(clout(end).XYZmm,cl(i),1)';
            clout(end).Z = ones(1,size(clout(end).XYZ,2));
            
            % keep track of all XYZmm coordinates used; remaining final cluster is
            % unused in any mask.
            [dummy,wh] = dominance_point_match(clxyz,XYZmm',0);
            clxyz(wh,:) = [];       % eliminate used points
            
            % Display
            
            cluster_orthviews(clout(end),{rand(1,3)},'add')
            spm_orthviews('Reposition',cl(i).mm_center)
            disp('Saved masked cluster in cl output structure.')
        end
    end
    
    gom = gom1;
    
    % Enter unused voxels in last cluster
    if ~isempty(clxyz)
         if isempty(clout), clout = cl(i);, else, clout(end+1) = cl(i);,end
        
         [clout(end).XYZmm] = dominance_point_match(round(cl(i).XYZmm'),clxyz,0)';
         clout(end).XYZ = mm2voxel(clout(end).XYZmm,cl(i),1)';
         clout(end).Z = ones(1,size(clout(end).XYZ,2));
    end
    
end

gom = input('Further subdivide clusters with anatomy? (1/0)? ');
if gom
    clout = anat_subclusters(clout);
end
    

return
        
        
        
        
