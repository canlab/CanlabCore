function [xBF_hires, xBF] = fmri_spline_basis(TR, varargin)
% :Usage:
% ::
%
%     xBF = spline_hrf_basis(TR, optional args)
%
% :Inputs:
%
%   **imTR:**
%        repetition time; sampling resolution of data
%
%   **'plot':**
%        optional: plot basis set
%
%   **'nbasis':**
%        optional: number of knot points
%
%   **'order':**
%        optional: order of spline model (# matched derivatives)
%
%   **'length':**
%        optional: length of window to model, in seconds
%
% Outputs:
% :Inputs:
%
%   **xBF.dt:**
%        time bin length {seconds}
%
%   **xBF.name:**
%        description of basis functions specified
%
%   **xBF.length:**
%        window length (seconds)
%
%   **xBF.order:**
%        order
%
%   **xBF.bf:**
%        Matrix of basis functions
%
%        32 second long spline basis set for fmri model
%
%   **xBF_hires:**
%        Sampled at high resolution, TR * 16
%
%   **xBF:**
%        Sampled at TR
%
% :Examples:
% ::
%
%    [xBF_hires, xBF] = fmri_spline_basis(TR, varargin)
%    [xBF_hires, xBF] = fmri_spline_basis(2, 'length', 12, 'nbasis', 3, 'order', 3, 'plot');


doplot = 0;
nBas = 8; % number of basis functions; higher = more spatial resolution. Must be at least k.
k = 4;  % spline order; higher = spline fit curvature matched on more derivatives; 2 = linear; 3 = quadratic; 4 = cubic (default)
windowlen = 32; % in sec

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            % functional commands
            case 'plot', doplot = 1;
            case 'nbasis', nBas = varargin{i+1};
            case 'order', k = varargin{i+1};
            case 'length', windowlen = varargin{i+1};
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


n_elements = ceil(windowlen ./ TR); % how many elements to span with basis


bf = Bspline(1:n_elements,k,[ones(1,k-1) linspace(1,n_elements+0.1,nBas-1) (n_elements+0.1)*ones(1,k-1)]);

xBF = struct('name', 'bspline_basis', 'dt', TR, 'bf', bf, 'length', windowlen, 'order', '8 bf, 4th order');

% high-res one

n_elements = ceil(16 .* windowlen ./ TR); % how many elements to span with basis

bf = Bspline(1:n_elements,k,[ones(1,k-1) linspace(1,n_elements+0.1,nBas-1) (n_elements+0.1)*ones(1,k-1)]);

xBF_hires = struct('name', 'bspline_basis', 'dt', TR ./ 16, 'bf', bf, 'length', windowlen, 'nbasis', nBas, 'order', k);

if doplot
    
    plot_inline1()

end


    function plot_inline1
            
        create_figure('basis functions', 2, 2);
        plot(bf, 'LineWidth', 2)
        title('Basis set');
        
        subplot(2, 2, 2);
        title('Simulated true HRF')
        hold on; 
        
        h = spm_hrf(TR./16) ./  max(spm_hrf(TR./16));
        plot((1:length(h)) ./ 16, h, 'k', 'LineWidth', 2)
        
        subplot(2, 2, 3);
        title('Example data (blue) and fits (red), hi-res')
        hold on; 
        
        h = h + .3 * noise_arp(length(h), [.5 .1]);
        plot((1:length(h)) ./ 16, h)
        
        f = xBF_hires.bf * pinv(xBF_hires.bf) * h(1:size(xBF_hires.bf, 1));
        plot((1:length(f)) ./ 16, f, 'Color', [1 .5 0])
        
        subplot(2, 2, 4);
        title('Example data (blue) and fits (red), sampled at TR')
        hold on; 
        
        y_at_tr = h(1:16:end);
        plot((1:length(y_at_tr)) - 1, y_at_tr, 'bo-', 'LineWidth', 2)
        
        f_at_tr = xBF.bf * pinv(xBF.bf) * h(1:16:size(xBF_hires.bf, 1));
        plot((1:length(f_at_tr)) - 1, f_at_tr, 'ro-', 'LineWidth', 2)
        
    end

end  % main fcn
   
