function histogram(obj)


% ---------------------------------------------------------------
% Histogram
% ---------------------------------------------------------------
Xtmp = obj.dat(:);

[h, x] = hist(Xtmp, 100);
han = bar(x, h);
set(han, 'FaceColor', [.3 .3 .3], 'EdgeColor', 'none');
axis tight;
xlabel('Values'); ylabel('Frequency');
title('Histogram of values');
drawnow


end

