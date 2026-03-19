function [sim_values, Nvox, d, low_agreement] = jackknife_similarity(obj, varargin)
% jackknife_similarity: Jackknife spatial similarity / agreement with leave-one-image-out reference
%
% Computes similarity (or agreement) between each image j and a leave-one-out
% reference computed from the remaining images (N − j).
%
% The reference image for image j is:
%
%     ref_j = median(dat(:, setdiff(1:k, j)), 2)
%
% where dat is the voxel × image matrix (or mean if specified).
%
% Similarity is then computed between image j and ref_j using the selected
% metric.
%
% :Usage:
% ::
%
%     [sim_values, Nvox, d, low_agreement] = jackknife_similarity(obj, 'similarity_metric', 'correlation')
%
% :Inputs:
%
%   **obj:**
%        An image_vector object or numeric matrix [voxels × images]
%
% :Optional Inputs:
%
%   **'similarity_metric'** (string):
%        Similarity or agreement metric to compute. Larger values indicate
%        better agreement.
%
%        Options:
%            'correlation' (default)
%            'cosine_similarity'
%            'dot_product'
%            'dice'
%            'absolute_agreement'
%            'concordance_correlation'
%            'standardized_abs_deviation'
%            'mean_shift_z'
%            'scale_shift_z'
%
%   **'treat_zero_as_data'** (logical):
%        Treat zeros as valid data instead of missing values.
%        Default = false.
%
%   **'complete_cases'** (logical):
%        Restrict calculations to voxels valid for all images.
%        Default = false.
%
%   **'verbose'** (logical):
%        Display progress and summary information.
%        Default = true.
%
%   **'doplot'** (logical):
%        Plot similarity values.
%        Default = false.
%
%   **'plot'**
%        Convenience keyword equivalent to 'doplot', true.
%
% :Outputs:
%
%   **sim_values:**
%        [k × 1] similarity/agreement values.
%
%   **Nvox:**
%        [k × 1] number of voxels used in each calculation.
%
%   **d:**
%        Standardized group-separation index (mean / std of sim_values).
%
%   **low_agreement:**
%        Logical vector flagging images >3 MAD below the median similarity.
%
% :Author:
%   2026, Tor Wager. GPL v3.
%

% -------------------------------------------------------------------------
% Handle object
% -------------------------------------------------------------------------

if isa(obj, 'image_vector')
    dat = obj.dat;
else
    dat = obj;
end

% -------------------------------------------------------------------------
% Allow 'plot' shorthand
% -------------------------------------------------------------------------

doplot = false;
plot_idx = strcmpi(varargin, 'plot');
if any(plot_idx)
    doplot = true;
    varargin(plot_idx) = [];
end

% -------------------------------------------------------------------------
% Input parsing
% -------------------------------------------------------------------------

p = inputParser;
p.addRequired('obj');
p.addParameter('similarity_metric','correlation', ...
    @(x) ismember(x, {'correlation','cosine_similarity','dot_product', ...
                      'dice','absolute_agreement','concordance_correlation', ...
                      'standardized_abs_deviation','mean_shift_z','scale_shift_z'}));
p.addParameter('treat_zero_as_data',false,@(x) islogical(x) || isnumeric(x));
p.addParameter('complete_cases',false,@(x) islogical(x) || isnumeric(x));
p.addParameter('verbose',true,@(x) islogical(x) || isnumeric(x));
p.addParameter('doplot',doplot,@(x) islogical(x) || isnumeric(x));

p.parse(obj,varargin{:});

sim_metric = p.Results.similarity_metric;
treat_zero_as_data = logical(p.Results.treat_zero_as_data);
complete_cases = logical(p.Results.complete_cases);
verbose = logical(p.Results.verbose);
doplot = logical(p.Results.doplot);

% -------------------------------------------------------------------------
% Preprocess
% -------------------------------------------------------------------------

if strcmp(sim_metric,'dice')
    treat_zero_as_data = true;
end

if ~treat_zero_as_data
    dat(dat==0) = NaN;
end

[nvox, k] = size(dat);

% -------------------------------------------------------------------------
% Stability constant and summaries
% -------------------------------------------------------------------------

dat_valid = dat(~isnan(dat));
global_abs_scale = median(abs(dat_valid));
if isempty(global_abs_scale) || global_abs_scale == 0
    global_abs_scale = 1;
end
epsilon = 1e-6 * global_abs_scale;

img_mean = mean(dat,1,'omitnan')';
img_centered = dat - img_mean';
img_scale = median(abs(img_centered),1,'omitnan')';

% -------------------------------------------------------------------------
% Initialize
% -------------------------------------------------------------------------

sim_values = nan(k,1);
Nvox = zeros(k,1);

idx = true(1,k);

% -------------------------------------------------------------------------
% Main loop
% -------------------------------------------------------------------------

for j = 1:k

    idx(j) = false;
    others = find(idx);

    ref = median(dat(:,others),2,'omitnan');
    x = dat(:,j);

    if strcmp(sim_metric,'standardized_abs_deviation')
        ref_scale = row_mad_omitnan(dat(:,others));
    else
        ref_scale = [];
    end

    valid = ~isnan(x) & ~isnan(ref);

    if complete_cases
        valid = valid & all(~isnan(dat(:,others)),2);
    end

    if strcmp(sim_metric,'standardized_abs_deviation')
        valid = valid & ~isnan(ref_scale);
    end

    a = x(valid);
    b = ref(valid);

    if strcmp(sim_metric,'standardized_abs_deviation')
        s = ref_scale(valid);
    end

    Nvox(j) = numel(a);

    if Nvox(j) < 2
        sim_values(j) = NaN;
        idx(j) = true;
        continue
    end

    switch sim_metric

        case 'correlation'
            sim_values(j) = corr(a,b);

        case 'cosine_similarity'
            den = norm(a)*norm(b);
            if den==0
                sim_values(j)=NaN;
            else
                sim_values(j)=dot(a,b)/den;
            end

        case 'dot_product'
            sim_values(j)=dot(a,b);

        case 'dice'
            a = a~=0; b=b~=0;
            denom = sum(a)+sum(b);
            if denom==0
                sim_values(j)=NaN;
            else
                sim_values(j)=2*sum(a & b)/denom;
            end

        case 'absolute_agreement'
            den = sum(abs(a)+abs(b));
            if den==0
                sim_values(j)=NaN;
            else
                sim_values(j)=1 - sum(abs(a-b))/den;
            end

        case 'concordance_correlation'
            ma=mean(a); mb=mean(b);
            va=var(a); vb=var(b);
            cab=cov(a,b); cab=cab(1,2);
            den = va+vb+(ma-mb)^2;
            if den==0
                sim_values(j)=NaN;
            else
                sim_values(j)=2*cab/den;
            end

        case 'standardized_abs_deviation'
            shift = mean(abs(a-b)./(s+epsilon));
            sim_values(j)=1/(1+shift);

        case 'mean_shift_z'
            others_mu = img_mean(others);
            mu_center = mean(others_mu,'omitnan');
            mu_spread = scalar_std_omitnan(others_mu);
            shift = abs(img_mean(j)-mu_center)/(mu_spread+epsilon);
            sim_values(j)=1/(1+shift);

        case 'scale_shift_z'
            others_scale = img_scale(others);
            scale_center = mean(others_scale,'omitnan');
            scale_spread = scalar_std_omitnan(others_scale);
            shift = abs(img_scale(j)-scale_center)/(scale_spread+epsilon);
            sim_values(j)=1/(1+shift);

    end

    idx(j) = true;

end

% -------------------------------------------------------------------------
% Effect size and outliers
% -------------------------------------------------------------------------

if std(sim_values,'omitnan')==0
    d = Inf;
else
    d = mean(sim_values,'omitnan')/std(sim_values,'omitnan');
end

low_agreement = sim_values < median(sim_values,'omitnan') - 3*mad(sim_values,1);

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------

if doplot
    figure;
    plot(sim_values,'ko-','LineWidth',1.5); hold on
    low_idx = find(low_agreement);
    plot(low_idx,sim_values(low_idx),'ro','LineWidth',2);
    yline(median(sim_values,'omitnan') - 3*mad(sim_values,1),'--');
    title(sprintf('Jackknife %s agreement (d=%.2g)',strrep(sim_metric,'_',' '),d));
    xlabel('Image'); ylabel('Agreement');
    grid on
end

% -------------------------------------------------------------------------
% Verbose
% -------------------------------------------------------------------------

if verbose
    fprintf('Jackknife computed for %d images\n',k);
    fprintf('d = %.2g\n',d);
    if any(low_agreement)
        fprintf('Low agreement images: %s\n',num2str(find(low_agreement)'));
    end
end

end

% -------------------------------------------------------------------------
% Helpers
% -------------------------------------------------------------------------

function m = row_mad_omitnan(X)
medx = median(X,2,'omitnan');
m = median(abs(X - medx),2,'omitnan');
end

function s = scalar_std_omitnan(x)
x = x(~isnan(x));
if numel(x)<2
    s = NaN;
else
    s = std(x);
end
end
