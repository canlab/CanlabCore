function R = local_rank_transform(M)
% local_rank_transform  Fallback for rsa.rdm.rankTransform when rsatoolbox absent.
%
% Rank-transforms an RDM/RSM into [0, 1]:
%   - Pulls the upper triangle (excluding diagonal).
%   - tiedrank() those values.
%   - Scale to [0, 1] via (rank - 1) / (n_pairs - 1).
%   - Re-symmetrize and set diagonal to 1.
%
% Matches the visual output of the Kriegeskorte rsatoolbox's
% rsa.rdm.rankTransform when used for heatmap display.
%
% Input
%   M  [k x k] numeric, symmetric (NaN diagonal allowed)
%
% Output
%   R  [k x k] numeric in [0, 1]

if isempty(M), R = M; return; end

k = size(M, 1);
mask = triu(true(k), 1);
vals = M(mask);

% Exclude NaN values from ranking but keep their positions
nan_pos = isnan(vals);
ranks   = nan(size(vals));
finite_vals = vals(~nan_pos);
if numel(finite_vals) > 1
    r = tiedrank(finite_vals);
    r = (r - 1) ./ (numel(finite_vals) - 1);
    ranks(~nan_pos) = r;
elseif ~isempty(finite_vals)
    ranks(~nan_pos) = 0.5;
end

R = zeros(k);
R(mask) = ranks;
R = R + R';
R(1:k+1:end) = 1;   % diagonal = 1 (matches workflow showRSM convention)

end
