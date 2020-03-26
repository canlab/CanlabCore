function y = time_varying_estimate(meth,data,varargin)
% Performs correlation operation (default) or any function you pass in (as a function handle) on symmetric moving average of data
%
% :Usage:
% ::
%
%     y = time_varying_estimate(meth,data,[window width],[function handle])
%
% - Works on all columns of the data input matrix together, so function handles can
%   operate on multivariate data (e.g. cond())
% - For moving averages, use the Matlab built-in smoothdata.m
%
% :Inputs:
%
%   **meth:**
%   Window-type options, either:
%   - 'gaussian'
%   - 'tukey'
%
% :Optional Inputs: 
%
%   **stepby:**
%
%   Sometimes input data is very high resolution, and it would take too much
%   time to work element by element across the inputs.  You can enter an
%   option here to compute estimates at every n-th lag.
%
% 	**Any function handle**
%   Any function handle to a function to evaluate on local rows of [data]
%
%   e.g., @corr
%
% :Examples:
% ::
%
%    y =  time_varying_estimate('gaussian',data,20);
%
%    % Generate sample data:
%    x = mvnrnd([0 0], [1 .6; .6 1], 100);
%
%    % Correlation between columns of x:
%    r = time_varying_estimate('tukey', x, 20);
% 
%    % St. deviation of first column
%    mystd = time_varying_estimate('tukey', x(:, 1), 20, @(y) std(y));
%
% ..
%    By Tor Wager
%    Last updated: Dec 2008
% ..

nshift = 0;
stepby = 1;                 % step size for shift
center_local_data = true;   % mean-center local data: good for corr, bad for moving average

% set up the kernel
% -------------------------------------------
switch meth
    
    case 'gaussian'
        
         ntrials = varargin{1};
         kern = normpdf(-3:6/ntrials:3); 
         kern = kern./max(kern);            % norm to max of 1
         
         mymax = find(kern == max(kern));
         nshift = mymax - 1;  % kernel shifts by n points; adjust
         
         kern = kern';  % column

    case 'tukey'
        % Window length is the zero-influence to zero-influence time
        kern = tukeywin(varargin{1});
        nshift = round(varargin{1} ./ 2);
        
        % kludgey adjust in case we need extra element
        %if length(i - shift : i + nshift) > length(kern)
            kern = [kern; 0];
        %end
    otherwise error('Unknown method.')
        
end

% set up function handle
% -------------------------------------------
if length(varargin) > 1
    fhandle = varargin{2};
else
    fhandle = @(y) my_corrcoef(y);
end

if length(varargin) > 2
    stepby = varargin{3};
end

% set up data
% -------------------------------------------
[nobs,ncols] = size(data);

% replicate kernel for each column
kern = repmat(kern,1,ncols);

if nobs < nshift, error('Not enough observations to support kernel.'); end
y = zeros(nobs,1);

% pad data at ends to avoid edge artifacts
% -------------------------------------------
paddat = data(end:-1:end-nshift,:);

data = [data; paddat];

paddat = data(nshift:-1:1,:);

data = [paddat; data];


% execute
% -------------------------------------------
for i = [(nshift+1):stepby:(nobs + nshift) (nobs + nshift)]
    % start at nshift to avoid ends; data is padded
    
    % set up windowed data
        
    dati = data(i - nshift : i + nshift,:);
    
    if center_local_data
        dati = scale(dati,1) .* kern;   % center and multiply by kern so data taper towards mean
    else
        dati = dati .* kern;
    end
    
    % execute: IN DEVELOPMENT: THIS is hard-coded for weighted correlation
    %r = weighted_corrcoef(dati,kern(:,1));
    %r = corrcoef(dati);
    
    y(i, :) = fhandle(dati);  %r(1,2);
    
    %y(:,i) = tmpy(nshift+1:nshift+nobs);
    
end

if stepby == 1
    y = y( (nshift+1):(nobs + nshift) );

else
    indx = [(nshift+1):stepby:(nobs + nshift) (nobs + nshift)];
    y = interp1(indx, y(indx), (nshift+1):(nobs + nshift));

end

end



function y = my_corrcoef(dati)
r = corrcoef(dati);
y = r(1,2);
end