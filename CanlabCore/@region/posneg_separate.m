function [pcl, ncl] = posneg_separate(cl, varargin)
% Separate a region object (cl) into clusters with positive and negative
% peak values, based on max (peak) value in .val or .Z field (default =
% val)
%
% [pcl, ncl] = posneg_separate(cl, ['Z'])
%
% Returns pcl and ncl, region structures with positive- and negative-valued
% peaks, respectively, copied from the original cl input.
%
% Optional input: 'Z', to use .Z field
%
% Note: You may have to use reparse_continguous to get this to work right.
% r = reparse_continguous(r);
% [pcl, ncl] = posneg_separate(r);

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
        
    else  % mixed.  remove the appropriate values from each.
        
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
