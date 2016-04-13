function [y, x] = resample_scnlab(data, p, q, varargin)
% Resample : Uses matlab's resample.m, but pads ends to avoid edge
% artifacts
%
% :Usage:
% ::
%
%     [y, x] = resample_scnlab(data, p, q)
%     % OR
%     [y, x] = resample_scnlab(data, p, q, origHz, targetHz)
%
% Y = RESAMPLE(X,P,Q) resamples the sequence in vector X at P/Q times
%     the original sample rate using a polyphase implementation.  Y is P/Q 
%     times the length of X (or the ceiling of this if P/Q is not an integer).  
%     P and Q must be positive integers.
%
% Other features:
%
% Returns x values for resampled data in original index scale
%
% IF two additional args are entered (origHz and targetHz), 
% p and q are determined automatically, based on your desired sampling rate
% (targetHz)
%
% :Example:
% ::
%
%    create_figure('test'); plot(y);
%    [y2, x] = resample_scnlab(y, 1, 5);
%    plot(x, y2, 'r');
%
%    %Use target Hz...take 100 Hz vector and resample at 20 Hz
%    create_figure('test');
%    plot(y); [y2, x] = resample_scnlab(y, [], [], 100, 20);
%    plot(x, y2, 'r');

if ~isempty(varargin)
    if length(varargin) < 2, error('Enter both origHz and targetHz or neither.'); end
    if ~isempty(p) || ~isempty(q), disp('Warning: p and q values entered will not be used.  Using origHz and targetHz.  Enter empty p and q to avoid this warning.'); end
    
    origHz = varargin{1};
    targetHz = varargin{2};
    
    p = 1; q = origHz / targetHz;
    
    while (abs(p - round(p)) > eps  || abs(q - round(q)) > eps)
        p = p * 10, q = q * 10;
    end
end
    
% origHz = 1000; targetHz = 50;
% origHz / targetHz
% p = 1, q = origHz / targetHz
% q = origHz / targetHz
% origHz = 925; targetHz = 50;
% p = 1, q = origHz / targetHz
% p = p * 10, q = q* 10


y = [];
x = [];

if isempty(data), disp('No data to resample.'); return, end
    
[nobs,ncols] = size(data);

if nobs == 1 && ncols > 1
    disp('Transposing data...seems like you''ve entered a row vector. Pass in column vectors.')
    data = data';
    
    [nobs,ncols] = size(data);
end
    
nshift = ceil(nobs ./ 10);

nshiftds = ceil(nshift * p / q); % downsample


% filter
% -------------------------------------------


if nobs < nshift, error('Not enough observations to support kernel.'); end
y = zeros(ceil(nobs * p / q) , ncols);


for i = 1:ncols
    
    % pad data at end to avoid edge artifacts
    padstart = data(nshift:-1:1, i);
    padend = data(end:-1:end-nshift, i);
    
    tmpy = resample([padstart; data(:,i); padend],p, q); 
    
    
    
    y(:,i) = tmpy(nshiftds+1:nshiftds + ceil(nobs * p / q));
    
end

x = linspace(1, nobs, (nobs) * p / q + 1)';

% quick fix
if length(x) > nobs
    x = x(1:nobs);
end


return
