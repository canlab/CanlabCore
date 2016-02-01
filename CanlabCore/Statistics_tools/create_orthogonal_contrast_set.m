function contrasts = create_orthogonal_contrast_set(nconditions)
% Create an orthogonal contrast set across nconditions conditions
%
% :Usage:
% ::
%
%     create_orthogonal_contrast_set(nconditions)
%
% Contrasts sum to zero, are orthogonal, and positive values and negative
% values each sum to 1.
%
% This is useful in evaluating designs, and for F-tests across all pairwise
% differences (e.g., one-way ANOVA contrast).
%
% ..
%    Tor Wager, Aug 2015
% ..

[row1, col1] = deal(zeros(1, nconditions));
col1(:) = 1;
row1(1:2) = [1 -1];

contrasts = toeplitz(col1, row1);

contrasts(end, :) = [];

% adjust values to sum to zero
for i = 2:nconditions-1
 mycon = contrasts(i, :);
 mycon(mycon ~= 0) = scale(mycon(mycon ~= 0)', 1)';  % mean-center non-zero vals
 mycon = mycon ./ sum(mycon(mycon > 0)); % normalize to sum to 1 on each side - "average" diff
 
 contrasts(i, :) = mycon;
end

% sum(contrasts, 2)  % contrasts sum to zero
% corr(contrasts') % contrasts are orthogonal
%contrasts = orth(contrasts')'


end


