function x = intercept_model(nvols_per_run, varargin)
% Build design matrix X for intercepts
% given vector of session lengths [s1 s2 s3] in images
%
% :Usage:
% ::
%
%     x = intercept_model(nvols_per_run, [indx of dummy scans in each session])
%
% :Examples:
% ::
%
%    nvols_per_run = [166 166 144 137];
%    x = intercept_model(nvols_per_run);
%
%    x = intercept_model(repmat(166, 1, 5));
%
%    Xi = intercept_model(EXPT.FIR.nruns, 1:2);
%
% ..
%    tor modified april 07: separate column for each run
% ..

    if ~isvector(nvols_per_run), error('nvols_per_run must be a vector of the number of volumes for each run'); end
    nvols_per_run = nvols_per_run(:)'; %ensure row vector
    
    npoints = sum(nvols_per_run); % number of time points total
    nruns = length(nvols_per_run); % number of runs

    x = zeros(npoints, nruns);

    st = cumsum([1 nvols_per_run]);
    en = st(2:end) - 1; % ending values
    st = st(1:end-1); % starting values

    for i = 1:nruns
        x(st(i):en(i), i) = 1;
    end

    if length(varargin) > 0 && ~isempty(varargin{1})
        % Model dummy regressors
        x = [x model_dummy(npoints, st', varargin{1})];
    end
end



function x = model_dummy(npoints, st, wh)
    % dummy regressors for first scan or two (at least, that's the intended
    % use)
    % separate column for each run!

    k = length(wh) .* length(st);
    x = zeros(npoints, k);

    dosamecol = 0;

    if dosamecol
        for i = 1:k
            ind = st(r) + wh(i) - 1; % which to filter
            ind(ind < 1) = []; % if negative numbers (end of runs), this makes it ok
            x(ind, i) = 1;
        end

    else

        colindx = 1;
        
        for r = 1:length(st) % for each run
            for i = 1:length(wh)
        
                ind = st(r) + wh(i) - 1; % which to filter
                ind(ind < 1) = []; % if negative numbers (end of runs), this makes it ok
                
                x(ind, colindx) = 1;
                colindx = colindx + 1;
            end
        end

    end

end
