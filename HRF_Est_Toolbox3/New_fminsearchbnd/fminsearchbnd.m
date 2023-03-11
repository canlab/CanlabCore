function [x,fval,exitflag,output]=fminsearchbnd4(fun,x0,LB,UB,options,varargin)
% FMINSEARCHBNDNEW: FMINSEARCH, but with bound constraints by transformation
%
% Changes from fminsearchbnd:
% 1) in options structure, user may pass an 'output function' and 'plot function' to fminsearch.
% Original fminsearchbnd handled the output function via a nested wrapper function.  I have extended
% this to the plot function too. 
% 2) I have moved the 'intrafun' and 'xtransform' functions and wrappers to be nested functions 
% (INSIDE the fminsearchbnd function), so they do not need to pass the params structure around 
% (into fminsearch) - but have access to it directly.  This maintains the integrity of the varargin, 
% which the user may be passing thru fminsearch to their optmization funciton (fminsearchbnd had 
% passed the params structure to fminsearch, thus ruining any varargin that the user passed in).
% This also obviates the params.(whatever) structure the author had, so I've eliminated it so things
% are simpler.
% 3) I have created a test example so the user can see not only how fminseachbnd works, but also how
% the OutputFn and PrintFns functions work, which were heretofore poorly documented by MathWorks.
% Many thanks to the original author, John D'Errico, for excellent work - very useful!
%
%   Modifications by: Ken Purchase
%              Email: kpurchase at yahoo
%               Date: 2007-Nov-29
%
%
% usage: x=FMINSEARCHBND(fun,x0)
% usage: x=FMINSEARCHBND(fun,x0,LB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options)
% usage: x=FMINSEARCHBND(fun,x0,LB,UB,options,p1,p2,...)
% usage: [x,fval,exitflag,output]=FMINSEARCHBND(fun,x0,...)
% 
% arguments:
%  fun, x0, options - see the help for FMINSEARCH
%
%  LB - lower bound vector or array, must be the same size as x0
%
%       If no lower bounds exist for one of the variables, then
%       supply -inf for that variable.
%
%       If no lower bounds at all, then LB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
%  UB - upper bound vector or array, must be the same size as x0
%
%       If no upper bounds exist for one of the variables, then
%       supply +inf for that variable.
%
%       If no upper bounds at all, then UB may be left empty.
%
%       Variables may be fixed in value by setting the corresponding
%       lower and upper bounds to exactly the same value.
%
% Notes:
%
%  If options is supplied, then TolX will apply to the transformed
%  variables. All other FMINSEARCH parameters should be unaffected.
%
%  Variables which are constrained by both a lower and an upper
%  bound will use a sin transformation. Those constrained by
%  only a lower or an upper bound will use a quadratic
%  transformation, and unconstrained variables will be left alone.
%
%  Variables may be fixed by setting their respective bounds equal.
%  In this case, the problem will be reduced in size for FMINSEARCH.
%
%  The bounds are inclusive inequalities, which admit the
%  boundary values themselves, but will not permit ANY function
%  evaluations outside the bounds. These constraints are strictly
%  followed.
%
%  If your problem has an EXCLUSIVE (strict) constraint which will
%  not admit evaluation at the bound itself, then you must provide
%  a slightly offset bound. An example of this is a function which
%  contains the log of one of its parameters. If you constrain the
%  variable to have a lower bound of zero, then FMINSEARCHBND may
%  try to evaluate the function exactly at zero.
%
%
% Example:
% rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
%
% fminsearch(rosen,[3 3])     % unconstrained
% ans =
%    1.0000    1.0000
%
% fminsearchbnd(rosen,[3 3],[2 2],[])     % constrained
% ans =
%    2.0000    4.0000
%
% See test_main.m for other examples of use.
%
%
% See also: fminsearch, fminspleas
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 4
% Release date: 7/23/06

% size checks
xsize = size(x0);
x0 = x0(:);
xLength=length(x0);

if (nargin<3) || isempty(LB)
  LB = repmat(-inf,xLength,1);
else
  LB = LB(:);
end
if (nargin<4) || isempty(UB)
  UB = repmat(inf,xLength,1);
else
  UB = UB(:);
end

if (xLength~=length(LB)) || (xLength~=length(UB))
  error 'x0 is incompatible in size with either LB or UB.'
end

% set default options if necessary
if (nargin<5) || isempty(options)
  options = optimset('fminsearch');
end


% 0 --> unconstrained variable
% 1 --> lower bound only
% 2 --> upper bound only
% 3 --> dual finite bounds
% 4 --> fixed variable
BoundClass = zeros(xLength,1);
for i=1:xLength
  k = isfinite(LB(i)) + 2*isfinite(UB(i));
  BoundClass(i) = k;
  if (k==3) && (LB(i)==UB(i))
    BoundClass(i) = 4;
  end
end

% transform starting values into their unconstrained
% surrogates. Check for infeasible starting guesses.
x0u = x0;
k=1;
for i = 1:xLength
  switch BoundClass(i)
    case 1
      % lower bound only
      if x0(i)<=LB(i)
        % infeasible starting value. Use bound.
        x0u(k) = 0;
      else
        x0u(k) = sqrt(x0(i) - LB(i));
      end
      
      % increment k
      k=k+1;
    case 2
      % upper bound only
      if x0(i)>=UB(i)
        % infeasible starting value. use bound.
        x0u(k) = 0;
      else
        x0u(k) = sqrt(UB(i) - x0(i));
      end
      
      % increment k
      k=k+1;
    case 3
      % lower and upper bounds
      if x0(i)<=LB(i)
        % infeasible starting value
        x0u(k) = -pi/2;
      elseif x0(i)>=UB(i)
        % infeasible starting value
        x0u(k) = pi/2;
      else
        x0u(k) = 2*(x0(i) - LB(i))/(UB(i)-LB(i)) - 1;
        % shift by 2*pi to avoid problems at zero in fminsearch
        % otherwise, the initial simplex is vanishingly small
        x0u(k) = 2*pi+asin(max(-1,min(1,x0u(k))));
      end
      
      % increment k
      k=k+1;
    case 0
      % unconstrained variable. x0u(i) is set.
      x0u(k) = x0(i);
      
      % increment k
      k=k+1;
    case 4
      % fixed variable. drop it before fminsearch sees it.
      % k is not incremented for this variable.
  end
  
end
% if any of the unknowns were fixed, then we need to shorten
% x0u now.
if k<=xLength
  x0u(k:xLength) = [];
end

% were all the variables fixed?
if isempty(x0u)
  % All variables were fixed. quit immediately, setting the
  % appropriate parameters, then return.
  
  % undo the variable transformations into the original space
  x = xtransform(x0u);
  
  % final reshape
  x = reshape(x,xsize);
  
  % stuff fval with the final value
  fval = feval(fun,x,varargin);
  
  % fminsearchbnd was not called
  exitflag = 0;
  
  output.iterations = 0;
  output.funcount = 1;
  output.algorithm = 'no call (all variables fixed)';
  output.message = 'All variables were held fixed by the applied bounds';
  
  % return with no call at all to fminsearch
  return
end


% Add the wrapper function to the user function right here inline:
    intrafun = @(x, varargin) fun(xtransform(x), varargin{:});

% Added code:  Add wrappers to output function(s) and plot function(s) - you can specify multiple
% output and/or print functions if you use a cell array of function handles.
    if ~isempty(options)
        % Add a wrapper to the output function(s) 
        % fetch the output function and put it(them) into a cell array:
        OutputFcn = createCellArrayOfFunctions(optimget(options,'OutputFcn',struct('OutputFcn',[]),'fast'),'OutputFcn');
        for ii = 1:length(OutputFcn)
            %stop = firstOutputFunction(OutStructure, optimValues, state, varargin)
            OutputFcn{ii} = @(x, varargin) OutputFcn{ii}(xtransform(x), varargin{:});
        end
        % store the "wrapped" output function back into the options.
        options = optimset(options, 'OutputFcn', OutputFcn);
        
        % Add a wrapper to the plot function(s) 
        % fetch the plot function and put it(them) into a cell array:
        PlotFcn = createCellArrayOfFunctions(optimget(options,'PlotFcns',struct('PlotFcns',[]),'fast'),'PlotFcns');
        for ii = 1:length(PlotFcn)
            %stop = firstOutputFunction(OutStructure, optimValues, state, varargin)
            PlotFcn{ii} = @(x, varargin) PlotFcn{ii}(xtransform(x), varargin{:});
        end
        % store the "wrapped" output function back into the options.
        options = optimset(options, 'PlotFcns', PlotFcn);
        % Add a wrapper to the print function(s) 
    end

    
   
% now we can call fminsearch, but with our own
% intra-objective function.
[xu,fval,exitflag,output] = fminsearch(intrafun,x0u,options,varargin);
output.algorithm = [output.algorithm ' bounded using fminsearchbnd'];

% undo the variable transformations into the original space
x = xtransform(xu);

% final reshape
x = reshape(x,xsize);


    % ======================================
    % ========= begin NESTED subfunctions =========
    % ======================================
        function xtrans = xtransform(x)
        % converts unconstrained variables into their original domains

        xtrans = zeros(xsize); %zeros(xLength, 1);   % I changed this to make it same dimension as the x in fminsearch
                                       % was zeros(1, params.xLength)
        % k allows some variables to be fixed, thus dropped from the
        % optimization.
        k=1;
        for i = 1:xLength
          switch BoundClass(i)
            case 1
              % lower bound only
              xtrans(i) = LB(i) + x(k).^2;

              k=k+1;
            case 2
              % upper bound only
              xtrans(i) = UB(i) - x(k).^2;

              k=k+1;
            case 3
              % lower and upper bounds
              xtrans(i) = (sin(x(k))+1)/2;
              xtrans(i) = xtrans(i)*(UB(i) - LB(i)) + LB(i);
              % just in case of any floating point problems
              xtrans(i) = max(LB(i),min(UB(i),xtrans(i)));

              k=k+1;
            case 4
              % fixed variable, bounds are equal, set it at either bound
              xtrans(i) = LB(i);
            case 0
              % unconstrained variable.
              xtrans(i) = x(k);

              k=k+1;
          end
        end

        end % sub function xtransform end




end % mainline end





