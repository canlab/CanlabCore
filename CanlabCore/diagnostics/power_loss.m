function power_loss(y, ons, X)

%  'true' model fit

% assume 'true' is FIR estimate
DX = tor_make_deconv_mtx3(ons, 20, 1);


% get residual variance and rss from best model
[bdx, fdx, rdx, rss, r_var, N, kdx] = fit_model(DX, y);

% if you had specified the right canonical shape at the start
Xbest = conv(ons, bdx); Xbest = Xbest(1:N);
Xbest(:, end+1) = 1;

% get residual variance and rss from best model
[bdx, fdx, rdx, rss, r_var, N, kdx] = fit_model(DX, y);


% E(t), expected t-value if you had gotten the HRF exactly right
% --------------------------------------------------------
c = eye(2);
c(:, end) = []; % remove intercept; assumes intcpt is at end!

[bbest, fbest, rbest, rss, r_var, N, kbest] = fit_model(Xbest, y);

[tbest, pbest] = get_t_values(c, bbest, Xbest, r_var);

delta = tbest;
alpha = .001;
df = N - kbest;
t_thresh = tinv(1 - alpha ./ 2, df);

power_best = 1 - nctcdf(t_thresh, df, delta);


% E(t), expected t-value with canonical HRF
% --------------------------------------------------------
hrf = spm_hrf(2);
CX = conv(ons, hrf); CX = CX(1:N);
CX(:, end+1) = 1;

% get residual variance and rss from canonical model
[bc, fc, rc, rss_c, r_var_c, N, kc] = fit_model(CX, y);

[tc, pc] = get_t_values(c, bc, CX, r_var_c);

delta = tc;
df = N - kc;
power_c = 1 - nctcdf(t_thresh, df, delta);


% Difference: Power loss
power_loss_val = 100 * (power_best - power_c) ./ power_best;

% significance for FIR model: full vs. reduced
% --------------------------------------------------------
Xred = ones(N, 1);

[F, p, dfb, dfe] = full_vs_reduced(y, DX, Xred, rss, N, kdx);

f_thresh = finv(1 - alpha, dfb, dfe);
fdelta = F;
power_fir = ncfcdf(f_thresh, dfb, dfe, fdelta);

power_loss_fir = 100 * (power_best - power_fir) ./ power_best;

% Plot
% --------------------------------------------------------
create_figure('plot'); plot(y, 'k'); hold on;
plot(fbest, 'r');
plot(fc, 'b');

fprintf(1,'With best (ideal HRF) model: t = %3.2f, p = %3.4f, power = %3.4f\n', tbest, pbest, power_best);
fprintf(1,'With canonical HRF model: t = %3.2f, p = %3.4f, power = %3.4f\n', tc, pc, power_c);
fprintf(1,'With FIR model F-test: F = %3.2f, p = %3.4f, power = %3.4f\n', F, p, power_fir);

fprintf(1, 'Power loss: Canonical vs. ideal: %3.2f%%\n', power_loss_val);
fprintf(1, 'Power loss: FIR vs. ideal: %3.2f%%\n', power_loss_fir);

end

%
% full vs. reduced
function [F, p, dfb, dfe] = full_vs_reduced(y, Xfull, Xred, rss, N, kfull)

[bred, fred, rred, rss_red] = fit_model(Xred, y);

dfb = (kfull - 1);
dfe = (N - kfull);

F = ( (rss_red - rss) ./ dfb ) ./ ( rss ./ dfe );

p = 1 - fcdf(F, dfb, dfe);


end


function [b, f, r, rss, r_var, N, k] = fit_model(X, y)
N = size(y, 1);
k = size(X, 2);


b = pinv(X) * y;
f = X * b;
r = y - f;

rss = r' * r;
r_var = rss ./ (N - k);

end



function [t, p, cb, var_cbeta] = get_t_values(c, b, X, r_var)


var_cbeta = diag( sqrt(r_var) * c' * inv(X' * X) * c );

% effect
cb = c' * b;

% t-values
t = cb ./ var_cbeta;

p = 2 * ( 1 - tcdf(abs(t), size(X,1) - size(X,2)) );

end
