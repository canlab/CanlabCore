function y = moving_average(meth,data,varargin)
% Symmetrical moving average filters
% (matlab's internal moving avg filter functions are asymmetric)
%
% Works on each column of the data input matrix
%
% :Examples:
% ::
%
%    y = moving_average('gaussian',data,fwhm)
%
%    % fwhm is not actually fhwm now...it's related to width though (temporary)
%    y = moving_average('gaussian',data,20);

nshift = 0;

% set up the kernel
% -------------------------------------------
switch meth
    
    case 'gaussian'
        
         ntrials = varargin{1};
         kern = normpdf(-3:6/ntrials:3); 
         kern = kern./sum(kern);
         
         mymax = find(kern == max(kern));
         nshift = mymax - 1;  % kernel shifts by n points; adjust
         

    otherwise error('Unknown method.')
        
end

% filter
% -------------------------------------------
[nobs,ncols] = size(data);

if nobs < nshift, error('Not enough observations to support kernel.'); end
y = zeros(nobs,ncols);


for i = 1:ncols
    
    % pad data at beginning and end to avoid edge artifacts
    padbegin = data(nshift:-1:1, i);
    paddat = data(end:-1:end-nshift,i);
    tmpy = conv([padbegin; data(:,i); paddat],kern); 
    
    y(:,i) = tmpy(nshift + nshift+1:nshift + nshift + nobs);
    
end



return
