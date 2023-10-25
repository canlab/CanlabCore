function [model, models]=plotHRF(HRF, t, varargin)
    % Generates a plot from an HRF Structure given a specified fit type and a region name.
    % Michael Sun, Ph.D.
    % - Takes HRF Structure object generated from EstimateHRF_inAtlas()
    % - t is a cellstring for fit type e.g., 'FIR', 'IL', or 'CHRF'
    % - r is a cellstring for region label from atlas.labels. e.g., 'ACC', or a cell-array of regions to compare relative to each other e.g., {'ACC', 'DLPFC'}  
    %
    % *Usage:
    % ::
    %    plotHRF(HRF_structure, 'FIR', 'CA2_Hippocampus_')

    % Check if the specified region 'r' and type 't' exist in HRF
    % Requires the MATlab Signal Processing Toolbox to run
    % detectPeaksandTroughs()

    % Flags to keep track of whether a cell array or atlas object is found
    isCellArrayFound = false;
    isAtlasFound = false;
    r=[];

    for k = 1:length(varargin)
        if ischar(varargin{k})
            r=varargin{k};

            if ~ismember(r, HRF.region)
                error('Invalid region specified.');
            end
    
        elseif iscell(varargin{k})
            r=varargin{k};
            isCellArrayFound = true;
        elseif isa(varargin{k}, 'atlas')  % assuming 'atlas' is a class you're checking for
            at=varargin{k};
            isAtlasFound = true;
        else
            disp(['Input argument ' num2str(k) ' unknown.'])
        end
    end

    if isempty(r) && ~isCellArrayFound
        r=HRF.region;
    end

    if ~isAtlasFound
        if isempty(HRF.atlas)
            at=load_atlas('canlab2018_2mm');
        else
            at=HRF.atlas;
        end
    end


    if ~ismember(t, HRF.types)
        error('Invalid type specified.');
    end

    if iscell(r)
        for i = 1:numel(r)
            if ~ismember(r{i}, HRF.region)
                error('Invalid region specified.');
            end
        end

        % disp(r)
        % at.get_region_volumes
    else
        r={r};
    end



    % Find the indices of the specified region and type
    reg = find(ismember(HRF.region, r));
    typ = find(ismember(HRF.types, t));

    % disp(reg)
    % disp(typ)

    % Get the list of conditions from HRF_PARAMS
    conds = HRF.params.CondNames;

    % % Initialize a cell array to store the model matrices
    % array3D = cell(1, length(conds));
    % 
    % % Populate the cell array with model matrices for each condition
    % for c = 1:numel(conds)
    %     for i = 1:numel(r)
    %         if isfield(HRF.fit{typ}{i}, conds{f}) && isfield(HRF.fit{typ}{i}.(conds{f}), 'model')
    %             array3D{c}{i} = HRF.fit{typ}{i}.(conds{f}).model;
    %         else
    %             error('Model data not found for condition: %s', conds{f});
    %         end
    %     end
    % end
    % 
    % % Check if any model matrices were found
    % if all(cellfun(@isempty, array3D))
    %     error('No model data found for the specified type and region.');
    % end
    % 
    % % Concatenate the model matrices along the third dimension
    % array3D = cat(3, array3D{:});

    % Initialize color map and legend entries
    % colors = jet(numel(conds));

    colors = cbrewer2('qual', 'Accent', numel(reg)+1);
    if numel(reg)>=4
        colors(4,:)=[]; % Remove yellow, too hard to see
    end
    model={};
    models={};

    % Create a figure
    figure;

    % Loop through each condition and plot
    for cond = 1:numel(conds)
        subplot(numel(conds), 1, cond);
        
        if numel(reg)==1
            model=HRF.fit{typ}{reg}.(conds{cond}).model;
            se=HRF.fit{typ}{reg}.(conds{cond}).wse;

            % plot(HRF.fit{typ}{reg}.(conds{cond}).model, '-')

            % detectPeaksTroughs(squeeze(array3D(:, :, cond))', true);
            detectPeaksTroughs(HRF.fit{typ}{reg}.(conds{cond}).model', true);
            hold on;

            % Plot error shade:
            x = 1:length(model);
            upper_bound = model + se;
            lower_bound = model - se;
            % Create x values for fill (concatenate forward and reverse x values)
            x_fill = [x, fliplr(x)];
            % Create y values for fill (concatenate upper_bound and reverse of lower_bound)
            y_fill = [upper_bound, fliplr(lower_bound)];
            fill_color = [0.7, 0.7, 0.7]; % Change to desired color
            fill_alpha = 0.3; % Transparency, change to desired value
            fill(x_fill, y_fill, fill_color, 'FaceAlpha', fill_alpha, 'EdgeColor', 'none');
            hold on;

            hline(0);

            if isfield(HRF.fit{typ}{reg}, conds{cond}) && isfield(HRF.fit{typ}{reg}.(conds{cond}), 'models')
                models=HRF.fit{typ}{reg}.(conds{cond}).models;

                % Plot all underlying sublines
                % for m = 1:height(HRF.fit{typ}{reg}.(conds{cond}).models)
                % 
                %     plot(HRF.fit{typ}{reg}.(conds{cond}).models(m,:), '-', 'Color', [0.9,0.9,0.9,0.2], 'LineWidth', 0.2);
                %     hold on;
                % end
   
            end

            region=format_strings_for_legend(r(1));
            region=region{1};
            title({['Condition ', conds{cond}, ' Fit-type: ', t, ' Region: ', region]}, 'Interpreter', 'none');
                

        else
            h=[]; % Linehandles for legend and labels
            for i = 1:numel(reg)
                [~, regionVoxNum, ~, ~]=at.select_atlas_subset(r(i), 'exact').get_region_volumes;
                if isfield(HRF.fit{typ}{i}, conds{cond}) && isfield(HRF.fit{typ}{i}.(conds{cond}), 'models')
                    models{i}=HRF.fit{typ}{i}.(conds{cond}).models/regionVoxNum;
                end

                model{i}=HRF.fit{typ}{reg(i)}.(conds{cond}).model/regionVoxNum;
                se{i}=HRF.fit{typ}{reg(i)}.(conds{cond}).wse/regionVoxNum;

                % Plot error shade:
                x = 1:length(model{i});
                upper_bound = model{i} + se{i};
                lower_bound = model{i} - se{i};
                % Create x values for fill (concatenate forward and reverse x values)
                x_fill = [x, fliplr(x)];
                % Create y values for fill (concatenate upper_bound and reverse of lower_bound)
                y_fill = [upper_bound, fliplr(lower_bound)];
                fill_color = colors(i,:); % Change to desired color
                fill_alpha = 0.3; % Transparency, change to desired value
                fill(x_fill, y_fill, fill_color, 'FaceAlpha', fill_alpha, 'EdgeColor', 'none');
                hold on;

                h(i)=plot(HRF.fit{typ}{reg(i)}.(conds{cond}).model/regionVoxNum, '-', 'Color', colors(i,:), 'DisplayName', char(format_strings_for_legend(r(i))));
                hline(0);
                title({['Condition ', conds{cond}, ' Fit-type: ', t, ' Regions: ', strjoin(r)], ['Error: ', 'within-subject SE']}, 'Interpreter', 'none');
                
                label(h(end), format_strings_for_legend(r(i)), 'location', 'left', 'FontWeight', 'bold', 'Margin', 3, 'HorizontalAlignment', 'right', 'FontSize', 14);
                label(h(end), format_strings_for_legend(r(i)), 'location', 'right', 'FontWeight', 'bold', 'Margin', 3, 'HorizontalAlignment', 'right', 'FontSize', 14);
                label(h(end), format_strings_for_legend(r(i)), 'location', 'center', 'FontWeight', 'bold', 'Margin', 3, 'HorizontalAlignment', 'right', 'FontSize', 14);
                label(h(end), format_strings_for_legend(r(i)), 'location', 'top', 'FontWeight', 'bold', 'Margin', 3, 'HorizontalAlignment', 'right', 'FontSize', 14);

                hold on;
            end
            % region_labels = format_strings_for_legend(r);
            % legend(region_labels);
            legend(h, format_strings_for_legend(r));
        end
    end


end






function plotRegionalHrfSummaries(HRF, regs)
    % data: Cell array with the data to be processed
    % conditions: Cell array with condition names
    % signalNames: Cell array with signal names
    % specificRegions: Cell array with specific regions for each signal and condition

    % Get specific regions for this signal and condition
    regs = getSpecificRegions(specificRegions, i, c);

    regs=HRF.region;

    % regs={'ACC'};
    % Generate a set of maximally different colors
    colormap('jet');  % Set colormap to parula
    colors = colormap;   % Get the colormap matrix
    % 
    % % Select maximally different colors
    nColors = numel(regs);
    indices = round(linspace(1, size(colors, 1), nColors));
    maxDifferentColors = colors(indices, :);
   
    % plot time-to-peak, height, width, start_time, and end_time of every
    % peaks_voxnormed and trough_voxnormed for a region

    create_figure(['Regional HRF summaries for ', HRF.params.CondNames{c}]);
    for c = 1:numel(HRF.params.CondNames)
        

        handles = [];
        subplot(2, 2, c)
        
        for r=1:numel(regs)
            
            
            % if c==3
            % 
            %     % mean time-to-peak
            %     meant=[mean(dat{i}(dat{i}.condition==[conds{c},'-cue'] & dat{i}.region==regs{r}, :).t), ...
            %                             mean(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).t)+20];
            % 
            %     % standard error time-to-peak
            %     stet=[std(dat{i}(dat{i}.condition==[conds{c}, '-cue'] & dat{i}.region==regs{r}, :).t)/sqrt(numel(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).t)), ...
            %                 std(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).t)/sqrt(numel(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).t))];
            % 
            %     % mean height
            %     meanh=[mean(dat{i}(dat{i}.condition==[conds{c},'-cue'] & dat{i}.region==regs{r}, :).h), ...
            %                             mean(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).h)];
            %     % height standard error
            %     steh=[std(dat{i}(dat{i}.condition==[conds{c},'-cue'] & dat{i}.region==regs{r}, :).h)/sqrt(numel(dat{i}(dat{i}.condition==[conds{c},'-cue'] & dat{i}.region==regs{r}, :).h)), ...
            %                             std(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).h)/sqrt(numel(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).h))];
            % 
            %     % minimum width
            %     minw=meant-[mean(dat{i}(dat{i}.condition==[conds{c},'-cue'] & dat{i}.region==regs{r}, :).w_times(:,1)), mean(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).w_times(:,1))+20];
            %     % maximum width
            %     maxw=[mean(dat{i}(dat{i}.condition==[conds{c},'-cue'] & dat{i}.region==regs{r}, :).w_times(:,2)), mean(dat{i}(dat{i}.condition==conds{c} & dat{i}.region==regs{r}, :).w_times(:,2))+20]-meant;
            % 
            %     handles(end+1)=errorbar(meant, meanh, steh, steh, stet, stet, 'o-', 'Color', maxDifferentColors(r,:));
            % 
            %     % Calculate the position and size of the rectangle
            %     rectX = meant(1) - minw(1);         % left boundary of the rectangle
            %     rectY = meanh(1) - max(steh(1));       % bottom boundary of the rectangle
            %     rectWidth = (meant(1)+maxw(1))-(meant(1)-minw(1)); % width of the rectangle
            %     rectHeight = 2*steh(1);   % height of the rectangle
            % 
            %     % Draw the rectangle
            %     % rectangle('Position', [rectX, rectY, rectWidth, rectHeight], 'FaceColor', colors(r,:), 'LineWidth', 1);
            %     patch([rectX, rectX+rectWidth, rectX+rectWidth, rectX], [rectY, rectY, rectY+rectHeight, rectY+rectHeight], maxDifferentColors(r,:), 'FaceAlpha', 0.2);
            % 
            %     % Calculate the position and size of the rectangle
            %     rectX = meant(2) - minw(2);         % left boundary of the rectangle
            %     rectY = meanh(2) - max(steh(2));       % bottom boundary of the rectangle
            %     rectWidth = (meant(2)+maxw(2))-(meant(2)-minw(2)); % width of the rectangle
            %     rectHeight = 2*steh(2);   % height of the rectangle
            % 
            %     % Draw the rectangle
            %     % rectangle('Position', [rectX, rectY, rectWidth, rectHeight], 'EdgeColor', 'r', 'LineWidth', 1);
            %     patch([rectX, rectX+rectWidth, rectX+rectWidth, rectX], [rectY, rectY, rectY+rectHeight, rectY+rectHeight], maxDifferentColors(r,:), 'FaceAlpha', 0.2);
            % 
            %     label(handles(end), format_strings_for_legend(regs{r}), 'location', 'right', 'FontWeight', 'bold', 'Margin', 3, 'HorizontalAlignment', 'right', 'FontSize', 14);
            %     label(handles(end), format_strings_for_legend(regs{r}), 'location', 'left', 'FontWeight', 'bold', 'Margin', 3, 'HorizontalAlignment', 'right', 'FontSize', 14);
            %     % label(handles(end), format_strings_for_legend(regs{r}), 'location', 'right', 'slope', 'FontWeight', 'bold', 'Margin', 3, 'HorizontalAlignment', 'right', 'FontSize', 14);
            %     % label(handles(end), format_strings_for_legend(regs{r}), 'location', 'top', 'slope', 'FontWeight', 'bold', 'Margin', 3, 'HorizontalAlignment', 'right', 'FontSize', 14);
            %     % label(handles(end), format_strings_for_legend(regs{r}), 'location', 'bottom', 'slope', 'FontWeight', 'bold', 'Margin', 3, 'HorizontalAlignment', 'right', 'FontSize', 14);
            % 
            % 
            % else
                % mean time-to-peak

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
                [~,vox,~,~]=get_region_volumes(at);
                vox=vox(r);
                se=barplot_get_within_ste(HRF.fit{1}{r}.(HRF.params.CondNames{c}).models);
                se=se/vox;

                se=repmat(se, 1, numel(meant));

                handles=[handles, errorbar(meant, meanh, se, se, repmat(0, 1, numel(meant)), repmat(0, 1, numel(meant)), 'o-', 'Color', maxDifferentColors(r,:))];

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

                
            % end

            hold on
        end

        hold off

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


% 
% barplot_columns()
% 
% model3d=cat(3, WASABIROIs_HRF_3.fit{1}{1}.hot.models, WASABIROIs_HRF_3.fit{1}{1}.warm.models, WASABIROIs_HRF_3.fit{1}{1}.imagine.models)
% 
% [se_within, stats]=barplot_get_within_ste(WASABIROIs_HRF_3.fit{1}{1}.hot.models)





