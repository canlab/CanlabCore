function clout = reparse_continguous(cl)
% Re-define regions in region object based on contiguous blobs
%
% :Usage:
% ::
%
%    clout = reparse_continguous(cl)
%
% ..
%    NEEDS SOME ADDITIONAL WORK/CHECKING
% ..

ivec = region2imagevec(cl);

clout = region(ivec, 'contiguous_regions', 'noverbose');

% ----------------------------------------------
% Add other fields
% ----------------------------------------------

clindx = ivec.volInfo.cluster;  % wh are in new clusters
nvox = length(clindx);

dat = cat(2, cl.all_data);

if size(dat, 2) == nvox
    
    % matching field; get dat for new clout
    for j = 1:length(clout)
        clout(j).all_data = dat(:, clindx == j);
        
        clout(j).dat = nanmean(clout(j).all_data', 1)';
    end
    
elseif isempty(dat)
    % do nothing
else
    disp('Warning!  all_dat field is wrong size.');
    
end

end



