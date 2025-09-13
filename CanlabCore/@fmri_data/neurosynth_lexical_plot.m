function [image_by_feature_correlations, top_feature_tables]=neurosynth_lexical_plot(fmri_data)
    % NEUROSYNTH_LEXICAL_PLOT Plots the correlations between fMRI data and lexical features.
    %
    % This function takes fMRI data and calculates the correlation between
    % the data and lexical features using the Neurosynth toolbox. It then
    % creates plots to visualize these correlations.
    %
    % USAGE:
    %   [image_by_feature_correlations, top_feature_tables] = neurosynth_lexical_plot(fmri_data, varargin)
    %
    % INPUTS:
    %   fmri_data : fMRI data to be analyzed.
    %
    % OUTPUTS:
    %   image_by_feature_correlations : Correlations between images and features.
    %   top_feature_tables            : Tables containing the top features and their scores.
    %
    % EXAMPLES:
    %   [correlations, tables] = neurosynth_lexical_plot(fmri_data);
    %
    % Authors: Michael Sun, Ph.D. 
    % Date: 05/30/2024


    % Future: Add a between option for separate image analysis instead of
    % replicates
    % if ~isempty(varargin) & strcmp(varargin{1}, 'between')
    %     [image_by_feature_correlations, top_feature_tables] = neurosynth_feature_labels(fmri_data, true, 'noverbose', 'display_output');
    % else
    %     [image_by_feature_correlations, top_feature_tables] = neurosynth_feature_labels(fmri_data, 'images_are_replicates', true, 'noverbose');
    % end

    [image_by_feature_correlations, top_feature_tables] = neurosynth_feature_labels(fmri_data, 'images_are_replicates', true, 'noverbose');
    % Aggregate results for plot
    lowwords = [top_feature_tables{1}.words_low(:)];
    % disp(lowwords);
    
    highwords = [top_feature_tables{1}.words_high(:)];
    % disp(highwords);
    
    r_low = top_feature_tables{1}.t_low;
    r_high = top_feature_tables{1}.t_high;
    
    r_to_plot = [r_high; r_low];
    textlabels = [highwords; lowwords];
    textlabels = strcat(textlabels, '  (', num2str(round(r_to_plot,1)), ')'); 
    
    figure;
    subplot(1,3,1)

    % Old version used Tor's Wedge Plot, but this doesn't display well most
    % of the time.
    % tor_wedge_plot(r_to_plot, textlabels, 'outer_circle_radius', .3, 'nofigure', 'colors', {[1 .7 0] [.4 0 .8]}, 'bicolor', 'labelstyle', 'close');

    hitbl=table(highwords, r_high, 'VariableNames', ["Word", "t"]);
    lotbl=table(lowwords, -1*r_low, 'VariableNames', ["Word", "t"]);
    
    % Horizontal Bar plot looks better:
    % Sort the values in descending order and reorder the text labels accordingly
    [sorted_values, sort_idx] = sort(r_to_plot, 'ascend');
    sorted_labels = textlabels(sort_idx);
    
    % Create a horizontal bar plot
    h = barh(sorted_values, 'FaceColor', 'flat');
    
    % Color the bars: green for positive values, red for negative values
    colors = repmat([0, 1, 0], length(sorted_values), 1);  % Initialize all as green
    colors(sorted_values < 0, :) = repmat([1, 0, 0], sum(sorted_values < 0), 1);  % Change to red for negative values
    h.CData = colors;
    
    % Remove the y-axis tick labels
    set(gca, 'YTick', []);
    
    % Add text labels at the end of the bars
    for i = 1:length(sorted_values)
        if sorted_values(i) >= 0
            text(sorted_values(i) + 0.1, i, sorted_labels{i}, 'VerticalAlignment', 'middle');
        else
            text(sorted_values(i) - 0.1, i, sorted_labels{i}, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
        end
    end
    
    % Optionally, you can adjust other properties like axis labels, title, etc.
    xlabel('t-score');
    % title('Horizontal Bar Plot with Sorted Bars and Custom Labels');

    subplot(1,3,2);
    wordcloud(hitbl,'Word','t', 'Color', 'red', 'HighlightColor', 'red');
    title("High Word Cloud");
    
    subplot(1,3,3);
    wordcloud(lotbl,'Word','t', 'Color', 'blue', 'HighlightColor', 'blue');
    title("Low Word Cloud");
    
    % Control Sizing
    set(gcf, 'position', [10, 10, 2000, 1600]);
    drawnow; snapnow;
end