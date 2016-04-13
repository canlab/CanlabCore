function [p,errval,fit,linkfun,fhan] = nonlin_fit(y,x,varargin)
% Multi-purpose nonlinear fitting to data
%
% :Usage:
% ::
%
%     [p,errval,fit] = nonlin_fit(y,x,['link',linktype],['err',objtype],['start',startp], ...
%                                ['plot'],['verbose'],[noverbose'],['quickstart'],['smartstart'])
%
% (Optional arguments are in [ ])
%
% :Inputs:
%
%   **y:**
%        data vector
%
%   **x:**
%        vector of x-values
%
%   **linktype:**
%        functional form to fit to data.  options are:
%
%        'sigmoid' 'exp0' 'exp2' 'exp3' 'exps'
%
%   **objtype:**
%        error measures to be minimized.  options are:
%
%        'sse' : sums of squared errors
%
%        'mad' : median absolute deviation (more robust to outliers)
%
%   **'plot':**
%        Create plot of data and fit
%
%   **startp:**
%        Vector of starting parameters for minimization
%
% :Functions:
%
%   **'sigmoid':**
%
%
%   **'exp1', 'exp2', 'exp3':**
%        Exponential function
%          - 'exp1'        y = e^(ax)          with free parameter a
%          - 'exp2'        y = b*e^(ax)        with free a, b
%          - 'exp3'        y = c + b*e^(ax)    with free a, b, c
%
%   **'pow0' 'pow2' 'pow3' 'pows1' 'pows2' 'pows':**
%       Power functions: 
%          - 'pows'        a + -(bx)^(-c)     with params a, b, c
%             simplifies to 1 - (1/x) with all p = 1
%             range: (for + b and c) -Inf at x = 0, to asymptote = a at x = Inf.
%
% :Usage:
%
% Default behavior: sigmoid fit, SSE objective:
% ::
%
%     [p,sse,fit] = nonlin_fit(y,x,'plot');
%
% Specify sigmoid function and SSE objective, starting at params [2 1 1]:
% ::
%
%     [p,sse,fit] = nonlin_fit(y,x,'start',[2 1 1],'link','sigmoid','err','sse');
%
% Specify median absolute deviation error minimization (more robust):
% ::
%
%     [p,sse,fit] = nonlin_fit(y,x,'err','mad','plot');
%
% Fit and add a fit line on a graph
% ::
%
%     [p,errval,fit,linkfun,fhan] = nonlin_fit(y,x,'linktype','exps','start',[1.2 .5 .8]);
%     fitx = 1:.1:5; fitline = fhan(p,fitx);
%     hold on; plot(fitx,fitline,'r');
%
% :Examples:
% ::
%
%    % sigmoid fit to y
%    % Generate fake data:
%    x = -5:.01:5;
%    sigmoid = inline('p(1) .* ( 1 ./ (1 + p(2)*exp(-p(3)*x)) )','p','x');
%    y = sigmoid([2 1 1.5],x);
%    y = y + .3 .* randn(1,size(y,2)) ;
%
%    % Fit, starting at param values [1 1 1]:
%    [p,sse,fit] = nonlin_fit(y,x,'plot','start',[1 1 1]);
%
%    % Bit more complicated data
%    x2 = -5:.01:5;
%    y = sigmoid([2 1 .5],x2);
%    y = [y 2*ones(1,100)];
%    y = y + .3 .* randn(1,size(y,2)) ;
%    x = 1:length(y);
%
%    % Exponential function (see also exp1, exp2, exp3 keywords)
%    expfun = inline('p(1)*exp(p(2)*x)','p','x');
% 
%    % simulated sample data
%    x = -1:.01:5;
%    y = expfun([.5 .7], x);                               
%    y = y + 1 .* randn(1,size(y,2));
%
%    % function handle to pass in
%    funhandle = @(p, x) expfun(p, x);
%    [p,sse,fit] = nonlin_fit(y,x,'linktype', funhandle, 'plot','start',[1 1]);
%
%    % Then, try to fit one-parameter exponential model to the same data:
%    [p,sse,fit] = nonlin_fit(y,x,'linktype', 'exp1', 'plot','start', 1);
%
% ..
%    Author and date / revision info:
%    Tor Wager, Dec. 6, 2006   % created
%    Tor        Dec. 10        % last modified
%
%    to-do ideas:
%    rsquare and other output stats in output structure
%    convergence flag : run 10 times w/random starting points and assess
%    convergence
%    sobol : sample space of starting params with sobol sequence
%    verbose output table with verbose flag
% ..


[linktype,objtype,startp,doplot,doverbose,dosmartstart] = setup_inputs; % setup defaults and inputs
startp;


% specify link function, transform of x to fit to y
% _________________________________________________________________
[fhan, linkfun] = get_link_function();


% specify error (objective) function
% function to be minimized
% _________________________________________________________________

objfun = get_objective_function();



% get starting parameter estimates:
% input, default, or 'smartstart' (specific for function)

startp = get_start_estimates(startp);

% run minimization of errfun
% find parameters p
% _________________________________________________________________
options = optimset('Display','off','MaxIter',100000); % ,%'LevenbergMarquardt','on');

PROBLEM = struct('objective',objfun,'x0',startp,'options',options,'solver','fminsearch');

[p,errval,exitflag,output] = fminsearch(PROBLEM);
% Note: output could be used for verbose reporting and problem checking


% get additional outputs: fits, errors, etc.
% _________________________________________________________________
fit = fhan(p,x);

if ~isreal(objfun(p)), disp('Error is not a real number!'); end
if isinf(objfun(p)), disp('Error is infinite!'); end
if exitflag == 0, disp('Max function evaluations reached!'); end

if doplot

    figure('Color','w');
    plot(x,y,'o','MarkerFaceColor',[.5 .5 .5],'MarkerSize',6);

    % this only works for some functions
    fitx = min(x):.1:max(x); fitline = fhan(p,fitx);
    hold on; plot(fitx,fitline,'k','LineWidth',2);

end


%%% END MAIN FUNCTION %%%













% Nested functions


% _________________________________________________________________

% setup defaults and inputs
% _________________________________________________________________
    function [linktype,objtype,startp,doplot,doverbose,dosmartstart] = setup_inputs

        doplot = 0;
        doverbose = 1;
        dosmartstart = 1;
        linktype = 'sigmoid';
        objtype = 'sse';
        startp = [1 1 1];

        for i = 1:length(varargin)
            if isstr(varargin{i})
                switch lower(varargin{i})
                    % reserved keywords
                    case 'plot', doplot = 1;
                    case {'smart','smartstart'}, dosmartstart = 1;
                    case {'quick','quickstart'}, dosmartstart = 0;
                    case {'verb','verbose'}, doverbose = 1;
                    case {'noverb','noverbose'}, doverbose = 0;

                        % functional commands
                    case {'link','linktype'}, linktype = varargin{i+1}; varargin{i+1} = [];
                    case 'err', objtype = varargin{i+1}; varargin{i+1} = [];
                    case 'start', startp = varargin{i+1};


                    otherwise, warning(['Unknown input string option:' varargin{i}]);
                end
            end
        end

    end

if doverbose
    % check Matlab version
    ok = checkMatlabVersion('8/3/2006');
end


% _________________________________________________________________

% specify link function, transform of x to fit to y
% _________________________________________________________________
%
    function [fhan, linkfun] = get_link_function

        if strcmp(class(linktype),'function_handle')
            % it's already a custom function handle.
            fhan = linktype;
            linkfun = linktype;

        else
            % p is a vector of parameters; x is the x-data

            switch linktype

                case 'sigmoid'
                    % sigmoid, 3-parameter: amplitude, scale (shift), exponent (steepness)
                    %linkfun = inline('p(1) .* ( 1 ./ (1 + p(2)*exp(-p(3)*x)) )','p','x');
                    % altered: ( 1 ./ (1 + exp(-p(3)*(x - p(2)))) )

                    % time-shifted sigmoid: params are amplitude, shift, and exponent
                    % (steepness)
                    linkfun = inline('p(1) .* ( 1 ./ (1 + exp(-p(3)*(x - p(2)))) )','p','x');

                    
                    
                case {'pow0','power0'}
                    linkfun = inline('x.^p(1)','p','x');

                case {'pow2','power2'}
                    linkfun = inline('p(1) .* x.^p(2)','p','x');

                case {'pow3','power3'}
                    linkfun = inline('p(1) + p(2) .* x.^p(3)','p','x');

                case {'pows1'}
                    linkfun = inline('3.5 + -(.005.*x).^-p(1)','p','x');

                case {'pows2'}
                    linkfun = inline('1 + -(p(1).*x).^-p(2)','p','x');

                case {'pows'}
                    linkfun = inline('p(1) + -(p(2).*x).^-p(3)','p','x');
                    % simplifies to 1 - (1/x) with all p = 1

                    
                case {'exp1','exponential1'}
                    linkfun = inline('exp(p(1)*x)','p','x');

                case {'exp2','exponential2'}
                    linkfun = inline('p(2) .* exp(p(1)*x)','p','x');

                case {'exp3','exponential3'}
                    linkfun = inline('p(3) + p(2) .* exp(p(1)*x)','p','x');
          

                otherwise
                    error('Unknown string command for objective function.')
            end

            % create a handle to the link function, to pass into error function
            %fhan = @(p) linkfun(p,x);
            fhan = @(p,x) linkfun(p,x);

        end
    end



% _________________________________________________________________

% specify error (objective) function
% function to be minimized
% _________________________________________________________________
    function objfun = get_objective_function

        switch objtype

            case 'sse'
                % sums of squared errors (standard)
                % this would be a way to call this function directly
                % sse = objfun(y,x,fhan,p);

                % Command-line call form:
                % objfun = @(y,x,fhan,p) sum((y - fhan(p,x)).^2);

                % fminsearch call form:
                objfun = @(p) sum((y - fhan(p,x)).^2);
                %objfun = @(p) sum((y - fhan(p,x)).^2);

                % Notes:
                % fminsearch needs a single-argument function
                % but if your objfun has multiple arguments, fminsearch
                % will use values of variables in existing workspace
                % as long as their names in the definition of objfun
                % are the same as those in the workspace.
                % that is why there are different command-line and fminsearch call forms
                % for objfun.

                % old
                %myerrfun = inline( 'sum((y - myfun(p,x)).^2)','y','x','myfun','p' );
                %objfun = @(p) myerrfun(y,x,fhan,p);

            case 'mad'
                objfun = @(p) mad(y - fhan(p,x));

            otherwise
                error('Unknown string command for objective function.')
        end

    end


% _________________________________________________________________

% get starting parameter estimates:
% input, default, or 'smartstart' (specific for function)
% _________________________________________________________________


    function startp = get_start_estimates(startp)

        % return if user-input
        if exist('startp','var')
            return
        end

%         if ~ischar(linktype)
%             startp = [];
%             return
%         end
        
        switch linktype

            case 'sigmoid'

                % quick-start, default
                startp = [1 1 1];

                if dosmartstart

                    len = length(y);
                    % get central 50% of x-values
                    iq = iqr(x)./2; med = median(x);
                    st = floor(med - iq); en = ceil(med + iq);

                    for i = st:en
                        dy(i) = sum(y(1:i))./i - sum(y(i+1:end))./(len-i);
                    end
                    [mmymax,wh] = max(abs(dy));


                    % shift param
                    startp(2) = x(wh);

                    % amplitude param
                    startp(1) = abs( sum(y(1:wh))./wh - sum(y(wh+1:end))./(len-wh) );

                end

            case {'exp3','exponential3','exps'}
                startp = [1 1 1];

            case {'exp2','exponential2','exps2'}
                startp = [1 1];

            case {'exp0','exponential0','exps1'}
                startp = [1];

            otherwise
                error('Unknown Link Type in get_start_estimates subfcn.')
        end
    end



end  % end main function





