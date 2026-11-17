function plotRegionalHrfSummaries(HRF, at, varargin)
    % data: Cell array with the data to be processed
    % conditions: Cell array with condition names
    % signalNames: Cell array with signal names
    % specificRegions: Cell array with specific regions for each signal and condition

    % Get specific regions for this signal and condition
    % regs = getSpecificRegions(specificRegions, i, c);

    % regs = getSpecificRegions(regs, i, c);

    if numel(varargin)>0
        if iscell(varargin{1})
            regs=varargin{1};
        else
            regs=HRF.region;
        end
    else
        regs=HRF.region;
    end

    % regs={'ACC'};
    % Generate a set of maximally different colors
    % colormap('jet');  % Set colormap to parula
    % colors = colormap;   % Get the colormap matrix
    
    % Get the colors
    colors = cbrewer2('qual', 'Accent', numel(regs)+3);
    
    % Define a threshold to identify yellow colors (you may need to adjust this)
    yellowThreshold = 0.8;
    
    % Identify the rows in 'colors' that correspond to yellow shades
    yellowRows = all(colors(:,1:2) > yellowThreshold, 2);
    
    % Remove the yellow colors
    colors(yellowRows, :) = [];

    % 
    % % Select maximally different colors
    nColors = numel(regs);
    indices = round(linspace(1, size(colors, 1), nColors));
    maxDifferentColors = colors(indices, :);
   
    % plot time-to-peak, height, width, start_time, and end_time of every
    % peaks_voxnormed and trough_voxnormed for a region

    % create_figure(['Regional HRF summaries for ', HRF.params.CondNames{c}]);
    figure;
    for c = 1:numel(HRF.params.CondNames)
        

        handles = [];
        subplot(2, 2, c)
        
        for r=1:numel(regs)

            meant=[];
            meanh=[];
            minw=[];
            maxw=[];

            for p=1:numel(HRF.fit{1}{r}.(HRF.params.CondNames{c}).peaks_voxnormed)

                % meant=mean(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).t);
                meant=[meant, HRF.fit{1}{r}.(HRF.params.CondNames{c}).peaks_voxnormed(p).time_to_peak];

                % standard error time-to-peak
                % stet=std(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).t)/sqrt(numel(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).t));

                % mean height
                % meanh=mean(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).h);
                meanh=[meanh, HRF.fit{1}{r}.(HRF.params.CondNames{c}).peaks_voxnormed(p).height];

                % height standard error
                % steh=std(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).h)/sqrt(numel(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).h));

                % minimum width
                % minw=meant-mean(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).w_times(:,1));
                % maximum width
                % maxw=mean(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).w_times(:,2))-meant;
                
                minw=[minw, HRF.fit{1}{r}.(HRF.params.CondNames{c}).peaks_voxnormed(p).start_time];
                    maxw=[maxw, HRF.fit{1}{r}.(HRF.params.CondNames{c}).peaks_voxnormed(p).end_time];

            end

            for t=1:numel(HRF.fit{1}{r}.(HRF.params.CondNames{c}).troughs_voxnormed)

                % meant=mean(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).t);
                meant=[meant, HRF.fit{1}{r}.(HRF.params.CondNames{c}).troughs_voxnormed(t).time_to_peak];

                % standard error time-to-peak
                % stet=std(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).t)/sqrt(numel(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).t));

                % mean height
                % meanh=mean(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).h);
                meanh=[meanh, HRF.fit{1}{r}.(HRF.params.CondNames{c}).troughs_voxnormed(t).height];

                % height standard error
                % steh=std(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).h)/sqrt(numel(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).h));

                % minimum width
                % minw=meant-mean(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).w_times(:,1));
                % maximum width
                % maxw=mean(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).w_times(:,2))-meant;
                
                minw=[minw, HRF.fit{1}{r}.(HRF.params.CondNames{c}).troughs_voxnormed(t).start_time];
                    maxw=[maxw, HRF.fit{1}{r}.(HRF.params.CondNames{c}).troughs_voxnormed(t).end_time];

            end

            % handles(end+1)=errorbar(meant, meanh, steh, steh, stet, stet, 'o', 'Color', maxDifferentColors(r,:));
            [~,vox,~,~]=at.select_atlas_subset(regs(r), 'exact').get_region_volumes;
            % vox=vox(r);
            se=barplot_get_within_ste(HRF.fit{1}{r}.(HRF.params.CondNames{c}).models);
            se=se/vox;

            se=repmat(se, 1, numel(meant));

            % Sort and plot
            
            % Combine the data into a single matrix for sorting
            data = [meant', meanh', minw', maxw', se'];
            
            % Sort the data based on meant values (1st column)
            sortedData = sortrows(data, 1);
            
            % Extract the sorted data
            meantSorted = sortedData(:, 1)';
            meanhSorted = sortedData(:, 2)';
            minwSorted = sortedData(:, 3)';
            maxwSorted = sortedData(:, 4)';
            seSorted = sortedData(:, 5)';
            
            % Plot the sorted data
            % handles=[handles, errorbar(meant, meanh, se, se, repmat(0, 1, numel(meant)), repmat(0, 1, numel(meant)), 'o-', 'Color', maxDifferentColors(r,:))];
            handles = [handles, errorbar(meantSorted, meanhSorted, seSorted, seSorted, zeros(1, numel(meantSorted)), zeros(1, numel(meantSorted)), 'o-', 'Color', maxDifferentColors(r,:))];

            if numel(meant) > 1
                label(handles(end), format_strings_for_legend(regs{r}), 'location', 'left', 'FontWeight', 'bold', 'Margin', 3, 'HorizontalAlignment', 'right', 'FontSize', 14);
                label(handles(end), format_strings_for_legend(regs{r}), 'location', 'right', 'FontWeight', 'bold', 'Margin', 3, 'HorizontalAlignment', 'right', 'FontSize', 14);
            else                
                label(handles(end), format_strings_for_legend(regs{r}), 'location', 'right', 'FontWeight', 'bold', 'Margin', 3, 'HorizontalAlignment', 'right', 'FontSize', 14);
            end

            % % Calculate the position and size of the rectangle
            % rectX = meant - minw;         % left boundary of the rectangle
            % % rectY = meanh - max(steh);       % bottom boundary of the rectangle
            % rectY = meanh - meanh/10;           % placeholder for now       
            % rectWidth = (meant+maxw)-(meant-minw); % width of the rectangle
            % % rectHeight = 2*steh;   % height of the rectangle
            % rectHeight = 2*(meanh/10);   % placeholder for now
            % 
            % patch([rectX, rectX+rectWidth, rectX+rectWidth, rectX], [rectY, rectY, rectY+rectHeight, rectY+rectHeight], maxDifferentColors(r,:), 'FaceAlpha', 0.2);

            % Loop through each rectangle and plot
            for i = 1:numel(minwSorted)
                rectangleHeight=seSorted(i)*2;
                x = minwSorted(i);
                y = meanhSorted(i) - rectangleHeight / 2;  % Adjust y to center the rectangle on meanh
                width = maxwSorted(i) - minwSorted(i);
                height = rectangleHeight;
                
                % Plot the rectangle
                rectangle('Position', [x, y, width, height], 'EdgeColor', maxDifferentColors(r,:));
            end

            hold on

        end

        % hold off

        xlim([0,45]);
    
        hlin=refline(0,0);
        hlin.Color='k';
        
        vline(13/0.46, 'r--');        

        legend(handles, format_strings_for_legend(regs), 'Location', 'best');
        ylabel('Similarity Amplitude')
        xticks([0:5:45])
        xticklabels([0:5:45]*0.46)
        xlabel('Time to Peak (seconds)');
        % title({[signames{i}, ' Condition: ', conds{c}], 'HRF Summary Statistics', 'With Standard Errors, Stimulus offset is demarcated'})
        title({['Condition: ', HRF.params.CondNames{c}], 'HRF Summary Statistics', 'With Standard Errors, Stimulus offset is demarcated'})

    end
end

function regs = getSpecificRegions(specificRegions, signalIdx, conditionIdx)
    % Extract and/or calculate the specific regions for the given signal and condition.
    % specificRegions: Cell array defining specific regions for each signal and condition
    % signalIdx: Index of the current signal
    % conditionIdx: Index of the current condition
    
    % Example logic - adapt as per your actual requirements
    regs = specificRegions{signalIdx, conditionIdx};
end
