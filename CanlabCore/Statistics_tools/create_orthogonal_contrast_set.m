function contrasts = create_orthogonal_contrast_set(nconditions)
% Create an orthogonal contrast set across nconditions conditions
%
% :Usage:
% ::
%
%     contrasts = create_orthogonal_contrast_set(nconditions)
%
% Contrasts sum to zero, are orthogonal, and positive values and negative
% values each sum to 1.
%
% This is useful in evaluating designs, and for F-tests across all pairwise
% differences (e.g., one-way ANOVA contrast).
%
% contrasts:  (nconditions-1) x (nconditions) matrix; each ROW is a contrast
%
% sum(contrasts, 2)  % contrasts sum to zero
% corr(contrasts') % contrasts are orthogonal
% contrasts = orth(contrasts')';                % Make orthonormal
%
% % Create orthogonal set of codes to span space of 10 subjects
% % Start with a condition function, subjectno, with [observations x 1],
% % integers indicate subject identity
% subjectno = repmat(1:10, 1, 20)'; % Create 20 observations each for 10 simulated subjects
% indic = condf2indic(subjectno);
% contrasts = create_orthogonal_contrast_set(size(indic, 2))
% subject_basis_set = indic * contrasts'
% % or (subject_basis_set = orth(indic * contrasts'))
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




end


