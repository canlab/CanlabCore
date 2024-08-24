function wordcloud_stat(words, t_stats)
% wordcloud_stat Generates a word cloud with words colored based on their t-statistics.
%
% Usage:
%   wordcloud_stat(words, t_stats)
%
% Inputs:
%   words    - A cell array of words to be displayed in the word cloud.
%   t_stats  - A numeric array of t-statistics (or any statistic) corresponding to each word.
%
% Outputs:
%   A figure displaying the word cloud where the color of each word represents
%   the magnitude and sign of its corresponding t-statistic. Negative t-statistics
%   are colored from blue to white, positive t-statistics are colored from white to red.
%   Requires wordcloud() from the specgraph toolbox
%
% Example:
%   words = {'word1', 'word2', 'word3'};
%   t_stats = [2.3, -1.5, 0.7];
%   wordcloud_stat(words, t_stats);
%
% Notes:
%   - The function requires that the length of words and t_stats be the same.
%   - Words are sized proportionally to the absolute values of their t-statistics.
%   - If all t-statistics are negative, words will be colored blue.
%   - If all t-statistics are positive, words will be colored red.
%
% Created by: Michael Sun, Ph.D.
% Date: 08/08/2024
% Version: 1.0
%
% See also: wordcloud


    % Check if input lengths match
    if length(words) ~= length(t_stats)
        error('The length of words and t_stats must be the same.');
    end
   
    % Normalize t-statistics to the range [0, 1] for color mapping
    min_t = min(t_stats);
    max_t = max(t_stats);
    abs_t_stats = abs(t_stats);
    
    figure;
    if min_t < 0 & max_t > 0
        % Convert t-statistics to absolute values for the word sizes

        normalized_t = (t_stats - min_t) / (max_t - min_t);
    

        % Create a custom colormap from blue (for negative t-stats) to red (for positive t-stats)
        n = length(t_stats);
        colors = zeros(n, 3); % Initialize RGB color matrix
        for i = 1:n
            if t_stats(i) < 0
                % Interpolate from blue to white for negative t-stats
                colors(i, :) = [0, 0, 1] * (1 - normalized_t(i)) + [1, 1, 1] * normalized_t(i);
            else
                % Interpolate from white to red for positive t-stats
                colors(i, :) = [1, 1, 1] * (1 - normalized_t(i)) + [1, 0, 0] * normalized_t(i);
            end
        end

        % Create the word cloud with custom colors
        wordcloud(words, abs_t_stats, 'Color', colors, 'SizePower', 1);
    elseif min_t < 0 & max_t <0
        
        wordcloud(words, abs_t_stats, 'Color', 'blue', 'HighlightColor', 'blue', 'SizePower', 1);  
    elseif min_t > 0 & max_t > 0 
        
        wordcloud(words, abs_t_stats, 'Color', 'red', 'HighlightColor', 'red', 'SizePower',1);
    end

end