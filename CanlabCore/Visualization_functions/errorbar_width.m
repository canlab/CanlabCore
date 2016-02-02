function errorbar_width(h, x, interval)
% Work with errorbar.m: Adjust the width of errorbar
%
% :Usage:
% ::
%
%    errorbar_width(h, x, interval)
%
% :Inputs:
%
%   **h:**
%        errorbar graphic handle
%
%   **x:**
%        vector x, which is used in errorbar
%
%   **interval:**
%        e.g., [-.1 .1] or [0 0] 
%
% :Examples: you can see this output in 
% http://wagerlab.colorado.edu/wiki/doku.php/help/core/figure_gallery
% ::
%
%    x = 1:5; % x values
%    y = [32 40 55 84 130]; % mean
%    e = [6 6 6 6 6]; % standard error of the mean
%
%    create_figure(y_axis);
%    set(gcf, 'Position', [1   512   268   194]);
%    col = [0.3333    0.6588    1.0000];
%    markercol = col-.2;
%
%    h = errorbar(x, y, e, 'o', 'color', 'k', 'linewidth', 1.5, 'markersize', 7, 'markerfacecolor', col);
%    hold on;
%    sepplot(x, y, .75, 'color', col, 'linewidth', 2);
%    errorbar_width(h, x, [0 0]); % here
%
%    set(gca, 'xlim', [.5 5.5], 'linewidth', 1.5);
%
%    try
%       pagesetup(gcf);
%       saveas(gcf, 'example.pdf');
%    catch
%       pagesetup(gcf);
%       saveas(gcf, 'example.pdf');
%    end
%
% ..
%    Copyright (C) 2014  Wani Woo
% ..

xdata = [];
for i = x
    xdata = [xdata repmat(i,1,2) NaN (repmat(i,1,2) + interval) NaN (repmat(i,1,2) + interval) NaN];
end

hh = get(h, 'children'); set(hh(2), 'XData', xdata);

end
