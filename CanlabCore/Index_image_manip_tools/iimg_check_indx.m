% Checks the index of an image
function [dattype,dat] = iimg_check_indx(dat,volInfo,varargin)
% :Usage:
% ::
%
%     [dattype,dat] = iimg_check_indx(dat,volInfo,[outtype])
%
% :Inputs:
%
%   **dattype:**
%        is 'full' or 'masked'
%
% :Optional Input:
%
%   **outtype:**
%        'full' or 'masked' output
%
% ..
%    tor wager
% ..


% make sure data is correct format

if isstr(dat), error('iimg_check_indx does not work with image files.  use iimg_read_img first'); end

% make sure volInfo fields exist

if ~isfield(volInfo,'nvox') || ~isfield(volInfo,'n_inmask')
    error('iimg_check_indx requires nvox and n_inmask fields for volInfo.  use iimg_read_img with extended output = 1 or 2')
end

% check length of dat and return dattype

if size(dat,1) == volInfo.nvox
    dattype = 'full';
elseif size(dat,1) == volInfo.n_inmask
    dattype = 'masked';
else error('data vector does not match size of image in volInfo!')
end

% format output
if length(varargin) > 0 && nargout > 1
    outtype = varargin{1};  % full or masked

    if ~strcmp(dattype, outtype)
        switch outtype
            case 'masked'
                dat = dat(volInfo.wh_inmask);
            case 'full'
                dat2 = zeros(volInfo.nvox,1);
                dat2(volInfo.wh_inmask) = dat;
                dat = dat2;
            otherwise
                error('Unknown output type: choose ''masked'' or ''full''')
        end
    end
end

return
