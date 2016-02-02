function actcolor = map_data_to_colormap(datavaluesets, poscm, negcm, varargin)
% Usage
% ::
%
%    actcolor = map_data_to_colormap(datavaluesets, poscm, negcm, varargin)
%
% Given sets of data values (each cell is a row vector of data values,
% e.g., z-scores) and color maps for positive and negative values,
% returns mapped colors for each data value in order.
% These colors can be used for direct plotting.
%
% :Inputs:
%
%   input 1:
%        data values (e.g., z-scores). k data value sets, in cells.  Each cell contains row vector of data values
%
%   input 2/3:
%        color maps [n x 3] for positive and negative values
%
%   input 4:
%        optional: fixed range of data defining max and min colors
%
% :Examples:
% ::
%
%    poscm = colormap_tor([0 0 0], [1 1 0]);
%    negcm = colormap_tor([0 0 1], [0 0 0]);
%    Z = randn(40, 1)';
%    actcolors = map_data_to_colormap({Z}, poscm, negcm)
%    [Z' actcolors{1}]
%
% ..
% Tor Wager, Sept. 2007
% ..

    nposcolors = size(poscm, 1);
    nnegcolors = size(negcm, 1);

    % -------------------------------------------------------------
    % determine overall data range
    % -------------------------------------------------------------

    % k data value sets, in cells.  Each cell contains row vector of data
    % values
    datavalues = cat(2,datavaluesets{:});

    % exactly zero values will give wrong length, so make pos. small number
    datavalues(datavalues == 0) = 100 * eps;
    
    posvalues = datavalues(datavalues > 0);
    negvalues = datavalues(datavalues < 0);
   
    if ~isempty(posvalues)

        zrange = [min(posvalues) max(posvalues)];

        % input fixed data range for max and min colors
        if ~isempty(varargin) && ~isempty(varargin{1})
            zrange = varargin{1}(1:2);

        end

        zh =  linspace(zrange(1),zrange(2),nposcolors);

        if isempty(zh), zh = [1 1 0]; end   % only one element?
    end

    if ~isempty(negvalues)
        zrangec = [min(negvalues) max(negvalues)];

        % input fixed data range for max and min colors
        if ~isempty(varargin) && ~isempty(varargin{1})
            zrangec = varargin{1}(3:4);
        end

        zhc =  linspace(zrangec(1),zrangec(2),nnegcolors);

        if isempty(zhc), zhc = [0 0 1]; end   % only one element?
    end

    % -------------------------------------------------------------
    % find color for each xyz
    % -------------------------------------------------------------
    
    if sum(isnan(datavalues))
        disp('Warning! NaNs in data values mapped to colors.  These will be mapped to black [0 0 0].');
    end
    
    for k = 1:length(datavaluesets)
        datavalues = datavaluesets{k};
        
        actcolor{k} = zeros(length(datavalues), 3);

        for i = 1:length(datavalues)

            dv = datavalues(i);

            if dv < 0
                [mydistance, wh] = min(abs(zhc - dv), [], 2);
                actcolor{k}(i,:) = negcm(wh, :);
                
            elseif dv >= 0

                [mydistance, wh] = min(abs(zh - dv), [], 2);
                actcolor{k}(i,:) = poscm(wh, :);

            else
                % could be NaN; leave as 0

            end
        end
    end

end
