function P = check_valid_imagename(P, varargin)
% :Usage:
% ::
%
%     P = check_valid_imagename(P, [return error])
%
% :Optional Input:
%
%   0 to return empty output, 1 to break with an error message,
%
%   2 to return list of bad/missing images, or
%
%   nothing to select missing filenames graphically

if isempty(P)
    if varargin{1}
        error(sprintf('No images to check.\n'));
    elseif ~isempty(varargin)
        return
    end
end

wasbad = false(size(P, 1), 1);

for i = 1:size(P,1)
    
    img = deblank(P(i, :));
    
    % take off trailing commas (SPM notation for multiple vols)
    wh = find(img == ',');
    if ~isempty(wh)
        wh = wh(end);
        img(wh:end) = [];
    end
    
    if ~(exist(img, 'file'))
        
        if ~isempty(varargin) && varargin{1} == 2
            wasbad(i, 1) = 1;
            
        elseif ~isempty(varargin) && varargin{1} == 1
            error(sprintf('Cannot find image:\n%s\n', P(i, :)));
            
        elseif ~isempty(varargin)
            % do nothing - return empty P
            P = [];
            return
            
        else
            fprintf(1,'Cannot find image:\n%s\nPlease select correct directory.\n',P(i,:));
            
            dd = spm_select(1,'dir','Select directory.');
            dd = [dd filesep];
            
            [myd] = fileparts(P(i,:));
            len = length(myd);
            
            P = [repmat(dd,size(P,1),1) P(:,len+1:end)];
            
            disp(P(i,:));
            
        end
        
    end
end % list of images

if ~isempty(varargin) && varargin{1} == 2
    fprintf('%3.0f OK images, %3.0f bad/missing images\n', sum(~wasbad), sum(wasbad));
    if any(wasbad)
        fprintf('Missing images:\n');
        fprintf('%3.0f ', find(wasbad))
        fprintf('\n');
    end
end

end
