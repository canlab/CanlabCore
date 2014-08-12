function subplots(data, varargin)
    figureTitle = '';
    numRows = [];
    numCols = [];
    xAxis = [];
    linkingAxes = 1;


    for i=1:length(varargin)
        if(ischar(varargin{i}) || isscalar(varargin{i}))
            switch(varargin{i})
                case 'Title'
                    figureTitle = varargin{i+1};
                case 'NumRows'
                    numRows = varargin{i+1};
                case 'NumCols'
                    numCols = varargin{i+1};
                case 'LinkingAxes'
                    linkingAxes = varargin{i+1};
                case {'XAxis', 'xAxis'}
                    xAxis = varargin{i+1};
            end
        end
    end

    if(isempty(numCols) && isempty(numRows))
        numRows = size(data,2);
        numCols = 1;
    elseif(~isempty(numCols) && isempty(numRows))
        numRows = ceil(size(data,2)/numCols);
    end

    figure('Name',figureTitle);
    %subplot(numRows, numCols, 1);
    if(isempty(xAxis))
        xAxis = (1:length(data));
    end
    ah = [];
    for i=1:size(data,2)
        ah(i) = subplot(numRows, numCols, i);
        plot(xAxis, data(:,i));
    end
    if(linkingAxes)
        linkaxes(ah);
    end
end