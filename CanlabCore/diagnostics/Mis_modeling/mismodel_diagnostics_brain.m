%function my_function_name(imgs, fixed_inputs, varargin)
%
% imgs         image names, in string matrix
% fixed_inputs one or more inputs that should always be entered
% varargin     optional inputs. Entered as 'name', value pairs of input
%               arguments
% 
% Will write the following files into the current directory:
% 
% 
function mismodel_diagnostics_brain(imgs, X, maskimg, varargin)

    % -------------------------------------------------------------
    % Define Optional Inputs Here
    % -------------------------------------------------------------
    
    % Default values for all optional inputs
    % ---------------------------------------
    nobs = size(imgs, 1);
    weights = ones(nobs, 1);        % default weights
    contrasts = [];
    contrastnames = {};

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % functional commands
                case 'contrasts', contrasts = varargin{i+1};

                case 'weights', weights = varargin{i+1};

                case 'contrastnames', contrastnames = varargin{i+1};

                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end
    
    % Check weights
    wtol = max(sqrt(weights))*eps^(1/3);
    if any(weights == 0) || any(weights < wtol)
        fprintf('Some weights are too low or badly scaled. Setting zero weights to %3.4f\n.', wtol);
        weights(weights < wtol) = wtol;
    end

    W = diag(weights);

    % [t, df, beta/contrast, Phi, sigma,stebeta, F] =
    % fit_gls(y,X,c,p,[PX, equal to pinv(X), for speed])


    % Map mask image into space of imgs, if not already, and write
    % 'mask.img'
    scn_map_image(maskimg, imgs(1, :), 'write', 'mask.img');
    maskimg = 'mask.img';
    
    % Setup inputs and function handle
    % -------------------------------------------------------------
    wh = intercept(X, 'which');
    if isempty(wh), disp('Warning! No intercept found in model.'); end

    %PX = pinv(X);
    PX = inv(X' * W * X) * X' * W;  % General weighted
    
    % Check if X is OK
    condx = cond(X);
    fprintf('Condition number of design matrix: %3.2f\n', condx);
    if condx > 10000, error('Design matrix is not estimable!'); end
    
    fhandle = @(y) fit_gls(y, X, contrasts, arorder, PX, weights);


    % Setup output image names
    % -------------------------------------------------------------
    % [beta, t, pvals, convals, con_t, con_pvals, sigma, Phi, df, stebeta, conste, F]
    
    %if isempty(contrasts), betaconname = 'beta'; else betaconname = 'contrast'; end

    names = cell(1, 12);
    names{7} = 'sigma.img';
    names{8} = 'phi.img';
    names{9} = 'df.img';
    names{12} = 'omnibus_F.img';

    for i = [1:6 10:11]  % Potential Multi-image outputs
        for j = 1:length(conditionnames)
            switch i
                case 1
                    names{i} = char(names{i}, ['beta_' conditionnames{j} '.img']);
                case 2
                    names{i} = char(names{i}, ['t_' conditionnames{j} '.img']);
                case 3
                    names{i} = char(names{i}, ['p_' conditionnames{j} '.img']);



                case 10
                    names{i} = char(names{i}, ['ste_beta_' conditionnames{j} '.img']);

            end
        end

        for j = 1:length(contrastnames)
            switch i
                case 4
                    names{i} = char(names{i}, ['con_beta_' contrastnames{j} '.img']);
                case 5
                    names{i} = char(names{i}, ['con_t_' contrastnames{j} '.img']);
                case 6
                    names{i} = char(names{i}, ['con_p_' contrastnames{j} '.img']);
                case 11
                    names{i} = char(names{i}, ['ste_con_' contrastnames{j} '.img']);
            end
        end

    end
    
    for i = [1:6 10:11]
        names{i} = names{i}(2:end, :);
    end
    
    save fit_gls_brain_SETUP imgs X arorder contrasts names conditionnames contrastnames maskimg

    % Run
    % -------------------------------------------------------------
    [beta, t, pvals, convals, con_t, con_pvals, sigma, Phi, df, stebeta, conste, F] = ...
        image_eval_function(imgs, fhandle, 'mask', maskimg, 'outimagelabels' , names);
    
    % other optional args: 'outimagelabels' , 'connames'

end
