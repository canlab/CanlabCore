function testFMinSearchNew
% Test function to show two things:
%
% - how the newly-modified fminsearchbnd works
% - how to use output functions and plot functions in fminsearch or fminsearchbnd (or other Matlab
% optimization routines) - heretofore documentation and examples for these have been sparse.
%
% note: use the fminsearchbnd that was modified 2007-Nov-29 by Ken Purchase, which handles plot and
% output functions properly.
%


    % set up the function that returns a sim structure (the simulation function)
    optFn = @(x, varargin) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
    
    x0 = [-1.2 1];
    mins = [-2 0]; %[2 -inf];
    maxs = [Inf Inf]; %[inf 3];

    % I'm going to pass in a useless extra parameter, just to show that fminsearchbnd and fminsearch
    % pass this back thru to your optimization function.  Your function could make use if this if
    % you watned.
    extraParams = 1;
    
    
    % Set up optimization options - you can leave any of these blank and fminsearch will use
    % defaults.
    searchOptions = struct(...
        'Display','none',...
        'MaxIter','200*numberOfVariables',...
        'MaxFunEvals','200*numberOfVariables',...
        'TolX',1e-6,...
        'TolFun',1e-6, ...
        'FunValCheck','off',...
        'OutputFcn', @firstOutputFunction,...  
        'PlotFcns',@firstPlotFunction);
    % NOTE: you could add several output or plot functions by incluing a cell array of function
    % handles, such as {@firstOutputFunction, @secondOutputFunction}
    
      
    % Run the optimization:
    [outX,fval,exitflag,output] = fminsearchbnd(optFn, x0, mins, maxs, searchOptions, extraParams);
    

    % Finally, re-run the best case thru the simulation function and display result:
    outX
    finalValue = optFn(outX, extraParams)
    
    
    
    
    % Define output and print functions.  These functions are nexted WITHIN the overall routine so they 
    % have access to variables in the above code if needed.
    %

    function stop = firstOutputFunction(xOutputfcn, optimValues, state, varargin)
        % create an output function for the fMinSearch
        %
        % inputs:
        % 1) xOutputfcn = the current x values
        % 2) optimValues - structure having:
        %         optimValues.iteration = iter;  % iteration number
        %         optimValues.funccount = numf;  % number of function eval's so far
        %         optimValues.fval = f;          % value of the function at current iter.
        %         optimValues.procedure = how;   % how is fminsearch current method (expand, contract, etc)
        % 3) State = 'iter','init' or 'done'     % where we are in the fminsearch algorithm.
        % 4) varargin is passed thru fminsearch to the user function and can be anything.
        %

        stop = false;

        % NOTE: this makes a bit of a messy display, but shows what you can do with an output
        % function.  You can get much of the same information using the fminsearch input option
        % 'Display', 'iter'
        disp(sprintf('Iteration: %d,  Evals: %d,  Current Min Value: %d', ...
            optimValues.iteration, optimValues.funccount, optimValues.fval));
        disp(['Best x so far: [' sprintf('%g ', xOutputfcn) ']']);
        
        % you could place plotting code here if you didn't want the automatic figure handling of the
        % plot functions.
        
        % you can also modify the value of 'stop' here to true if you want fminseach to terminate
        % based on any criteria you'd put here.
    end
        


    function stop = firstPlotFunction(xOutputfcn, optimValues, state, varargin)
        % create an print function for the fMinSearch
        %
        % NOTE: The plot functions do their own management of the plot and axes - if you want to
        % plot on your own figure or axes, just do the plotting in the output function, and leave
        % the plot function blank.  
        %
        % One thing the plot function DOES have it that it installs STOP and PAUSE buttons on the
        % plot that allow you to interrupt the optimization to go in and see what's going on, and 
        % then resume, or stop the iteration and still have it exit normally (and report output 
        % values, etc).  
        %
        % inputs:
        % 1) xOutputfcn = the current x values
        % 2) optimValues - structure having:
        %         optimValues.iteration = iter;  % iteration number
        %         optimValues.funccount = numf;  % number of function eval's so far
        %         optimValues.fval = f;          % value of the function at current iter.
        %         optimValues.procedure = how;   % how is fminsearch current method (expand, contract, etc)
        % 3) State = 'iter','init' or 'done'     % where we are in the fminsearch algorithm.
        % 4) varargin is passed thru fminsearch to the user function and can be anything.
        %
        
        stop = false;
        
        hold on; 
        % this is fun - it simply plots the optimization variable (inverse figure of merit) as it 
        % goes along, so you can see it improving, or stop the iterations if it stagnates.
        rectangle('Position', ...
            [(optimValues.iteration - 0.45) optimValues.fval, 0.9, 0.5*optimValues.fval]);
        set(gca, 'YScale', 'log');
        
        % when you run this, try pressing the 'stop' or 'pause' buttons on the plot.
        
        % you can add any code here that you desire.
        
    end
        
    
    
    
    
end % end of test code



