function varargout = matrix_eval_function(Y,fhandle,varargin)
% Evaluate any arbitrary function fhandle on each column of an input matrix Y
%
% :Usage:
% ::
%
%     varargout = matrix_eval_function(Y,fhandle,varargin)
%
% evaluate fhandle on paired columns of X and Y
%
% fhandle is a function handle:
% ::
%
%    fhandle = @(variable inputs) fit_gls(variable_input,fixed_inputs);
%    fhandle = @(y) fit_gls(y,X,c,p,PX);
%
% Function is fhandle should return (multiple possible) outputs, each
% in a row vector.  Matrix outputs of this function are Y-cols x output
% values for each output argument.
%
% :Example: Generalized least squares fitting on 100 Y-variables, same X
% ::
%
%    y = rand(100,100); X = y + rand(100,1); X(:,end+1) = 1; c = [1 0]'; p = 2; PX = pinv(X);
%    fhandle = @(y) fit_gls(y,X,c,p,PX);
%    [t, df, beta, Phi, sigma,stebeta, F] = fhandle(y);
%
% ..
%    tor wager, jan 31, 2007
% ..

    varargout = {};

    if nargout == 0, disp('No outputs requested.'); return, end
    
    % % -------------------------------------------------------------------
    % % Setup variable X and check
    % define n and v, # observations and # variables
    % % -------------------------------------------------------------------
    setup_inputs_and_check_sizes;
    
    whichGood = find(whichGood);
    
    if isempty(whichGood), error('No valid variables with non-zero variance to be analyzed.'); end
    
    t1 = clock;
    
    % do this to make sure these are returned so they may be used in later
    % inline
    n; v;

    % % -------------------------------------------------------------------
    % % build string to evaluate that defines outputs and inputs
    % % -------------------------------------------------------------------
    build_eval_string(nargout);

    % % -------------------------------------------------------------------
    % % initialize outputs with "Good" data vector (those that have
    % variability)
    %  we need to do this to get sizes of each output, which vary depending
    %  on the input function.
    %  then we can skip some vars that have no Y variance.
    %  also make sure we can store everything in memory.
    % % -------------------------------------------------------------------
    i = whichGood(1);
    eval(fstr);
    for i = 1:nargout
        varargout{i} = NaN .* zeros(v,size(varargout{i},2));
    end
    

    % % -------------------------------------------------------------------
    % % Evaluate vectors with "Good" data
    % % -------------------------------------------------------------------
       
    if doverbose
        % print banner and define text update points
        print_banner;
    end

    for i = whichGood
        % only run for Y-columns that have some variance; otherwise return
        % NaN

        % evaluate the function for this variable (column)
        eval(fstr);

        if doverbose
            % print string
            update = ( updateiterations == i );
            if any(update), fprintf(1,'\b\b\b\b%3d%%',updateperc(update)); end
        end

    end

    
    if doverbose
        fprintf(1,'\n_______________________________\n')
        endt = etime(clock,t1);
        endt = endt ./ 60;
        fprintf(1,'\nDone in %3.0f hrs, %3.0f min, with %3.3f s per variable \n',floor(endt./60),rem(endt,60),endt.*60./v); 
    end


    %END OF MAIN FUNCTION CODE


    % -------------------------------------------------------------------
    %
    %
    % INLINE FUNCTIONS
    %
    %
    %
    % -------------------------------------------------------------------

    % % -------------------------------------------------------------------
    % % Setup inputs, including variable X and check
    % % -------------------------------------------------------------------
    function setup_inputs_and_check_sizes
        dochk = 1;
        doverbose = 1;
        do_variable_X = 0;
        if length(varargin) > 0
            for i = 1:length(varargin)
                if ischar(varargin{i})
                    switch varargin{i}
                        case {'verbose','verb'}, doverbose = 1;
                        case {'noverb','noverbose'}, doverbose = 0;
                        case {'nochk'}, dochk = 0;
                        otherwise
                            disp('Warning: unrecognized string input argument.');
                    end
                else
                    X = varargin{i};
                    do_variable_X = 1;
                end
            end
        end

        [n,v] = size(Y);
        whichGood = true(1,v);  % OK analysis voxels
        
        if dochk 
            [val,whichGood] = no_variance(Y);
            if val, disp('Warning: Some Y vectors have no variability!'); end
        end
        
        if do_variable_X
            
            if isempty(X), error('X should be omitted or be matrix of size(Y)'); end
            [n2,k] = size(X);

            if dochk
                if n ~= n2, error('data and model sizes do not match.'); end
                if v ~= k, error('data and model sizes do not match.'); end
            end


        end
        
        

    end


    % %     fstr = ['[t, df, beta, Phi, sigma,stebeta, F] = fhandle(y);'];
    % %     eval(fstr)


    % % -------------------------------------------------------------------
    % % build string to evaluate that defines outputs and inputs
    % % -------------------------------------------------------------------
    function build_eval_string(numargs)
        fstr = '[';
        for arg = 1:numargs
            fstr = [fstr 'varargout{' num2str(arg) '}(i,:)'];
            if arg ~= numargs
                fstr = [fstr ','];
            else
                fstr = [fstr '] = fhandle(double(Y(:,i)));'];
            end
        end


    end


    % % -------------------------------------------------------------------
    % % print banner and define update points for text output
    % % -------------------------------------------------------------------
    function print_banner

        fprintf(1,'matrix_eval_function.m')
        fprintf(1,'\n_______________________________\n')
        fprintf(1,'Evaluating this function on %3.0f variables:\n',v);
        disp(fhandle);
        fprintf(1,'...using this command:\n%s\n',fstr);
        fprintf(1,'_______________________________\n')
        str = sprintf('Running ... Done  %03d%%',0); fprintf(1,str);
        updateiterations = 1:round(v ./ 100):v;
        updateperc = round(linspace(0,100,length(updateiterations)));
    end

end % END MAIN FUNCTION




% % -------------------------------------------------------------------
% % Sub-functions
% % -------------------------------------------------------------------

function val = no_intercept(X)

    val = all(any(diff(X,1,1),1));

end

function [val,whichGood] = no_variance(Y)

    whichGood = any(diff(Y,1,1),1);
    val = ~all(whichGood);

end

function erase_string(str1)
    fprintf(1,repmat('\b',1,length(str1))); % erase string
end

