function model = onsets2parametric_mod_X(ons, pm_vals, nscan, basisset, varargin)
% :Examples:
% ::
%
%    ons = ons{1}; %obj.Sess(1).U(1).ons ./ TR;
%    pm_vals = obj.Sess(1).U(1).P.P;
%    nscan = obj.nscan(s);

if ~iscell(ons)
    error('ons must be cell array.')
end

if ~iscell(pm_vals)
    error('pm_vals must be cell array.')
end

nconds = length(ons);

delta = double(onsets2delta(ons, nscan)) ;  % double from logical


% make sure basisset is column vector, if single vector
if size(basisset,1) == 1, basisset = basisset'; end
nbf = size(basisset, 2);


% get modulated delta function
for i = 1:nconds
    
    wh = find(delta(:, i));
    pm_u = scale(pm_vals{i}, 1);
    
    if length(wh) ~= length(pm_u)
        disp('Number of events in parametric modulator does not match number of events in delta function.')
        disp('Bad input?  Stopping in debugger so you can check.')
        keyboard
    end
    
    delta(wh, i) = pm_u;
end

% Convolve
for i = 1:nconds
    
    lendelta = size(delta, 1);
    X{i} = zeros(lendelta, nbf);
    
    for j = 1:size(basisset,2)
        tmp = conv(double(delta(:,i)), basisset(:,j)); % Changed for Matlab 7.9+ compatibility - Thanks, Liane
    
        X{i}(:, j) = tmp(1:lendelta);
    end
    
end

model = cat(2, X{:});


% downsample, if necessary
if ~isempty(varargin)
    
    dsrate = varargin{1};
    [n, k] = size(model);
    nt = n ./ dsrate;
    t = (1:n)';             % minor change to save time, tor: 10/12/10
    
    if nt ~= round(nt), error('Length of stimList is not evenly divisible by downsampling factor.'); end
    
    modeli = zeros(nt , k);
    
    for i = 1:size(model, 2)
        
        xi = 1:varargin{1}:n;       % downsample rate
        
        modeli(:, i) = interp1(t, model(:, i), xi, 'linear', 'extrap'); % tor: 10/12/10 : minor change: add extrap
    end
    
    model = modeli;
    %model = model(1:varargin{1}:end,:); % equivalent to resample(X,1,varargin{1},0)
end


model = model(1:nscan,:);      	% eliminate extra values



end

