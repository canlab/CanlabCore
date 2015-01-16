function CLU = clusters2CLU(clusters,varargin)
    % function CLU = clusters2CLU(clusters,[opt] M)
    %
    % Inputting an M matrix will transform the coordinates
    % by that M, to convert between voxel sizes, etc.
    %
    % by Tor Wager

    CLU = [];
    if isempty(clusters), return, end

    if ~isfield(clusters(1),'threshold'), clusters(1).threshold = 1; end

    if isfield(clusters(1),'Z_descrip'), CLU.Z_descrip = clusters(1).Z_descrip; end

    CLU.XYZmm = cat(2,clusters.XYZmm);
    CLU.mm_center = mean(CLU.XYZmm, 2)';
    CLU.XYZ = cat(2,clusters.XYZ);
    try
        CLU.Z = cat(2,clusters.Z);
    catch
        %CLU.Z = cat(1,clusters.Z)';
        for i = 1:length(clusters), if size(clusters(i).Z,1) > size(clusters(i).Z,2), clusters(i).Z = clusters(i).Z'; end, end
        CLU.Z = cat(2,clusters.Z);
    end
    CLU.title = clusters(1).title;

    CLU.u = clusters(1).threshold;
    CLU.threshold = clusters(1).threshold;

    CLU.voxSize = clusters(1).voxSize;
    CLU.VOX = clusters(1).voxSize;
    
    try
        CLU.P = strvcat(clusters.P);
        CLU.imP = strvcat(clusters.imP);
    catch
    end

    if isfield(clusters, 'all_data')
        CLU.all_data = cat(2, clusters.all_data);
    end
    
    if nargin > 1
        CLU.M = varargin{1};
        CLU = transform_coordinates(CLU,CLU.M);
        
    elseif ~isfield(clusters,'M')
        disp('Choose SPM.mat file with volume info (affine mat file)')
        [SPM,VOL,xX,xCon,xSDM] = spm_getSPM;
        CLU.M = VOL.M;
    else
        CLU.M = clusters(1).M;
    end

    CLU.numVox = size(CLU.XYZmm,2);

    if isempty(CLU.Z)
        disp('clusters Z field is empty.  Filling with ones as a placeholder.')
        CLU.Z = ones(1,size(CLU.XYZmm,2));
    end

    if size(CLU.Z,1) > 1, CLU.allZ = CLU.Z; CLU.Z = CLU.Z(1,:); end
end