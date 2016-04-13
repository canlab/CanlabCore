function [nanvec, varargout] = nanremove(varargin)
% :Usage:
% ::
%
%     varargout = nanremove(varargin)
%     [nanvec, x1, x2, etc.] = nanremove(x1,x2,etc...)
%
% removes cases that are NaN in any column of any variable x1...xn
%
% ..
%    tor wager
%    update: 7/10/2007 to handle multi-column matrix inputs
% ..

    if nargin == 0, error('nanremove: no input args.'), end

    m = size(varargin{1}, 1);
    nanvec = false(m, 1);

    % get vector of any nans in any columns of inputs
    for i = 1:length(varargin)
        if isempty(varargin{i})
            % skip

        else
            if m ~= size(varargin{i}, 1), error('Data vectors are not equal length.'); end
            nanvec = nanvec | any(isnan(varargin{i}), 2);
        end
    end

    % remove from all outputs
    for i = 1:length(varargin)
        if isempty(varargin{i})
            varargout{i} = [];
        else
            varargout{i} = varargin{i}(~nanvec, :);
        end
    end
end
