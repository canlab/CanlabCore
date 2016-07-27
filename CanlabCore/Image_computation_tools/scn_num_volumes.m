% :Usage:
% ::
%
%     [n, V] = scn_num_volumes(V)
%
% :Input:
%
%   **V:**
%        is spm_vol structure, or image filename
%
% :Output:
%
%   **V:**
%        is spm_vol structure
%
% counts number of images n stored in a 4-D file
% .img files can contain multiple volumes (indexed by .n in spm_vol structure).
% how many are in this volume?
%
% ::
%
%    % Adapted (copied) from code by Tom Nichols.
%     function n = scn_num_volumes(V)
% 
% 
%    if ~isstruct(V)
%         V = spm_vol(V);
%        spm_close_vol(V);
%    end
% 
%    fp   = fopen(V.fname);
%    fseek(fp,0,'eof');
%    Len  = ftell(fp);
%    fclose(fp);
%    switch(spm('Ver'))
% %     case 'SPM2'
% %         mydt = V.dim(4);
% %     case 'SPM5'
% %         mydt = V.dt(1); 
% %     otherwise
% %         error('Unknown SPM version "%s": neuroscientists of the future, fix me!', spm('Ver'));
% %     end
% %   n = Len/(prod(V.dim(1:3))*spm_type(mydt,'bits')/8);
% % end
%
% search index words: number of volumes, nvol
%
% .img files can contain multiple volumes (indexed by .n in spm_vol structure).
% how many are in this volume?

function [n, V] = scn_num_volumes(V)
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

        case {'SPM5' 'SPM8' 'SPM12'}
            n = length(V);

        otherwise
            error('Unknown SPM version "%s": neuroscientists of the future, fix me!', spm('Ver'));
    end

end
