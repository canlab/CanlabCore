function [pcl, ncl] = posneg_separate(cl, varargin)
% Separate a region object (cl) into clusters with positive and negative
% peak values, based on values in .val or .Z field (default = val)
% If a region has both positive and negative values, it will be included in
% both sets pcl (positive) and ncl (negative), with only positive-valued and
% negative-valued subsets included in each, respectively.
%
% :Usage:
% ::
%
%    [pcl, ncl] = posneg_separate(cl, ['Z'])
%
% Returns pcl and ncl, region structures with positive- and negative-valued
% peaks, respectively, copied from the original cl input.
%
% :Optional Input:
%
%   **Z:**
%        To use .Z field
%
% :Note: You may have to use reparse_continguous to get this to work right.
% ::
%
%    r = reparse_continguous(r);
%    [pcl, ncl] = posneg_separate(r);
%

reparseflag = 0;  % need to reparse if mixed clusters
myfield = 'val';

if any(strcmp(varargin, 'Z'))
    myfield = 'Z';
end

pindx = 1;
nindx = 1;
pcl = [];
ncl = [];

for i = 1:length(cl)
    
    vals = cl(i).(myfield);
    
    [dummy, wh] = max(abs(vals));
    
    % wh voxels in region have positive/negative effects
    sv = sign(vals);
    whp = sv > 0;
    whn = sv < 0;
    
    if all(whp)
        % positive values
        if isempty(pcl)
            pcl = cl(i);
        else
            pcl(pindx) = cl(i);
        end
        pindx = pindx + 1;
        
    elseif all(whn)
        % negative values
        if isempty(ncl)
            ncl = cl(i);
        else
            ncl(nindx) = cl(i);
        end
        nindx = nindx + 1;
        
    else  % mixed.  add to both. remove the appropriate voxels from each.
        
        reparseflag = 1;
        
        if isempty(pcl)
            pcl = cl(i);
        else
            pcl(pindx) = cl(i);
        end
        
        if isempty(ncl)
            ncl = cl(i);
        else
            ncl(nindx) = cl(i);
        end
        
        whrem = whn;
        
        if ~isempty(pcl(pindx).XYZ), pcl(pindx).XYZ(:, whrem) = []; end
        if ~isempty(pcl(pindx).XYZmm), pcl(pindx).XYZmm(:, whrem) = []; end
        if ~isempty(pcl(pindx).val), pcl(pindx).val(whrem, :) = []; end
        if ~isempty(pcl(pindx).Z), pcl(pindx).Z(whrem) = []; end
        if ~isempty(pcl(pindx).all_data), pcl(pindx).all_data(:, whrem) = []; end
        
        whrem = whp;
        
        if ~isempty(ncl(nindx).XYZ), ncl(nindx).XYZ(:, whrem) = []; end
        if ~isempty(ncl(nindx).XYZmm), ncl(nindx).XYZmm(:, whrem) = []; end
        if ~isempty(ncl(nindx).val), ncl(nindx).val(whrem, :) = []; end
        if ~isempty(ncl(nindx).Z), ncl(nindx).Z(whrem) = []; end
        if ~isempty(ncl(nindx).all_data), ncl(nindx).all_data(:, whrem) = []; end
        
        
        pindx = pindx + 1;
        nindx = nindx + 1;
    end
    
end

if reparseflag
    
    pcl = reparse_continguous(pcl);
    ncl = reparse_continguous(ncl);
    
end

end



% 
% function Z = get_signed_max(cl, myfield)
% % Returns a table with var "MaxZ" nregions x 1, or empty if cl.Z is empty
% 
% %maxZ = @(i) max(double(cl(i).Z));
% 
% smax = @(i) signedmax(cl(i).(myfield));
% 
% for i = 1:length(cl)
%     
%     if ~isempty(cl(i).(myfield))
%         myZ(i, 1) = smax(i);
%     else
%         myZ(i, 1) = NaN;
%     end
%     
% end
% 
% if all(isnan(myZ))
%     Z = [];
% end
% 
% end % function
% 
% function val = signedmax(vals)
% 
% vals = double(vals);
% [maxabs, wh] = max(vals);
% 
% val = sign(vals(wh)) .* maxabs;
% 
% end
% 
% 
