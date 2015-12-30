function [myfft, freq, handle] = fft_plot_scnlab(dat, TR, varargin)
% :Usage:
% ::
%
%     [myfft, freq, handle] = fft_plot_scnlab(dat, TR, varargin)
%
% :Inputs:
%
%   **dat:**
%        is a data vector (column)
% 
%   **TR:**
%        is the sampling rate of the data you put in
%        in seconds / sample, or 1/Hz
%
% :Optional inputs:
%   - 'samefig'
%   - 'color, ['b'] or other color
%   - 'bar'
%   - 'linebar': both line and bar
%
% :Examples:
% ::
%
%    % plot effects of filtering on a difference
%    % between two regressors
%    spm_hplength = SPM.xX.K.HParam;
%    d = SPM.xX.X * SPM.xCon(mycon).c(:, 1);
%    create_figure('Contrast'); plot(d) % contrast we care about
%    px = pinv(SPM.xX.K.X0);           % pinv of the filtering matrix
%    y = d;
%    y = y - SPM.xX.K.X0 * px * y;     % residuals after filtering
%
%    [myfft, freq, handle] = fft_plot_scnlab(d, 2);
%    hold on;
%    [myfft2, freq2, handle] = fft_plot_scnlab(y, 2); set(handle,'Color','r')
%    plot_vertical_line(1/spm_hplength)
%    set(ans, 'Color', 'b', 'LineWidth', 3)
%

fftfig = 1;
ptype = 'line';
color = 'k';

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'samefig', fftfig = 0;

            case 'bar', ptype = 'bar';
            case 'color', color = varargin{i+1}; varargin{i+1} = [];

            case 'linebar', ptype = 'both';
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end


n = length(dat);

% from matlab example
% power = abs(Y(1:floor(n/2))).^2;
% nyquist = 1/2;
% freq = (1:n/2)/(n/2)*nyquist

nyq = 1 ./ (2 * TR);

timepts = floor(n ./ 2);

freq = (0:timepts-1)/timepts * nyq;
%freq = linspace(0, nyq, timepts)';

myfft = fft(dat); %real(abs(fft(dat)));
myfft = abs(myfft(1:timepts)) .^ 2;  % power

myfft = myfft ./ sum(myfft);

if fftfig
    create_figure('fftplot',1, 1, 1);
else
    hold on
end

%myfft = myfft(1:timepts);
switch ptype
    case 'line'
        handle = plot(freq, myfft, '-', 'Color', color, 'LineWidth', 2);

    case 'bar'
        handle = bar(freq, myfft);
        set(handle, 'FaceColor', color);

    case 'both'
        handle = bar(freq, myfft);
        set(handle, 'FaceColor', color, 'EdgeColor', 'none');
        handle = plot(freq, myfft, '-', 'Color', color, 'LineWidth', 2);
        
    otherwise
        error('unknown plot type');
end

end

