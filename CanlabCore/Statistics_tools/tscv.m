function [trIdx, teIdx] = tscv(vectorlen, varargin)
% Create a crossvalidation test and train index for time series data.
% Larger h will ensure less dependence.  Larger v creates larger test sets
% Number of folds = vectorlen - (2*h + 2*v)
%
% -See http://robjhyndman.com/hyndsight/tscvexample/ for more info about rolling cv
% -See Racine, J. (2000). Consistent cross-validatory model-selection for dependent data: hv-block cross-validation. Journal of Econometrics, 99(1), 39-61.
%
% :Usage:
% ::
%
%    [trIdx, teIdx] = tscv(vectorlen, stepsize)
%
% :Inputs:
%
%   **vectorlen:**
%        length of vector to create holdout cross-validation set
%
% :Optional Inputs with their default values:
%
%   **'hvblock' = [h,v]:**
%        use hvblock cross-validation with a block
%        size of 'h' (0 reduces to v-fold xval)and
%        number of test observations 'v' (0 reduces
%        to h-block xval)
%
%   **'rolling' = [h,v,g]:**
%        use hvblock cross-validation with g training steps
%        surrounding hv block.  Akin to Rolling
%        crossval. Same properties as hvblock.
%
% :Outputs:
%
%   **trIdx:**
%        structure with training label index
%
%   **teIdx:**
%        structure with test label index
%
% :Examples:
% ::
%
%    [trIdx, teIdx] = tscv(100, 'hvblock',[5,2]); % use hvblock with h=5 and v=2
%    [trIdx, teIdx] = tscv(100, 'rolling',[5,2,10]); % use hvblock with h=5, v=2 and g=10
%
% ..
%    Original version: Copyright Luke Chang & Hedwig Eisenbarth 11/2013
%
%    Programmer's Notes:
%    LC 11/28/13:
%          -changed input and documentation
%          -Don't use matlab functions as variable names (e.g., median) -
%          median > mid
%          rewrote the hv block so that it loops though all available data -
%          need to finish coding the rollingcv option
%    LC & HE 12/16/13:
%          -added rollingcv option
%
%    HE & LC 11/19/14: 
%          -increased the test data used by adjusting how the training blocks work at the ends
% ..

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'hvblock'}
                xval_type = varargin{i};
                inval = varargin{i + 1};
                h = inval(1);
                v = inval(2);
                varargin{i} = [];
                varargin{i + 1} = [];
                
            case {'rolling'}
                xval_type = varargin{i};
                inval = varargin{i + 1};
                h = inval(1);
                v = inval(2);
                g = inval(3);
                varargin{i} = [];
                varargin{i + 1} = [];
        end
    end
end

switch xval_type
    case 'hvblock'
        %hv cross validation: leave completely out h steps around the test interval
        %See Racine, J. (2000). Consistent cross-validatory model-selection for dependent data: hv-block cross-validation. Journal of Econometrics, 99(1), 39-61.
        
        stepsize = 2*v + 2*h + 1;
        
        if stepsize > vectorlen;
            error('stepsize is too large, please decrease')
        end
        
        start = 1;
        xval_index = 1;
        while start <= vectorlen - stepsize + 1
            trIdx{xval_index} = true(vectorlen,1);
            teIdx{xval_index} = false(vectorlen,1);
            trIdx{xval_index}(start:(start + 2*h + 2*v)) = false; %train set = everything - 2*v + 2*h + 1
            teIdx{xval_index}((start + h):(start + h + 2*v)) = true; %test set = 2*v + 1
            start = start + 1 + 2*v;
            xval_index = xval_index + 1;
        end
        
    case 'rolling' % this needs to be fixed.
        %rolling cross validation: leave completely out h steps around the
        %test interval with g training steps
        %See http://robjhyndman.com/hyndsight/tscvexample/ for more info
        
        if g == 0
            error('g must be greater than 0')
        end
        
        stepsize = 2*g + 2*v + 2*h + 1;
        
        if stepsize > vectorlen;
            error('stepsize is too large, please decrease')
        end
        
        start = 1;
        xval_index = 1;
        while start <= (vectorlen - (2*v))
            
            %Initialize vector as zeros
            trIdx{xval_index} = false(vectorlen,1);
            teIdx{xval_index} = false(vectorlen,1);
            
            if start <= g+h+1 % adjustment for the g+h first - add points to training on other side of test
                if start <= h+1
                    trIdx{xval_index}((start + 2*v+1 + h) : (start + 2*v + h + 2*g)) = true;
                    teIdx{xval_index}(start:(start + 2*v)) = true; %test set = 2*v + 1
                else
                    trIdx{xval_index}((start - (start-1)) : (start - h-1)) = true;
                    trIdx{xval_index}((start + 2*v+1 + h) : (h + 2*v+1 + h + 2*g)) = true;
                    teIdx{xval_index}(start:(start + 2*v)) = true; %test set = 2*v + 1
                end
                
            elseif start >= (vectorlen - ((2*v) + g+h)) % adjustment for the g+h last timepoints  - add points to training on other side of test
                if start > (vectorlen - (2*v + h))
                    trIdx{xval_index}((start - (2*g + h)) : (start - h - 1)) = true; %(start - (g + h)) : (start - h - 1)
                    teIdx{xval_index}(start:(start + 2*v)) = true;
                else
                    trIdx{xval_index}((vectorlen - stepsize + 1) : (start - h -1)) = true;
                    trIdx{xval_index}((start + (2*v+1) + h) : vectorlen) = true;
                    teIdx{xval_index}(start:(start + 2*v)) = true;
                end
                
            else %Normal hvg block
                trIdx{xval_index}((start- (g+h)):(start - h - 1)) = true; %train set = everything - 2*v + 2*h + 1
                trIdx{xval_index}((start + (2*v+1) + h): (start + (2*v+1) + h + g -1)) = true; %train set = everything - 2*v + 2*h + 1
                teIdx{xval_index}(start:(start + 2*v)) = true;
                %teIdx{start}((start + g + h):(start + g + h + 2*v)) = true; %test set = 2*v + 1
            end
            start = start + 1 + 2*v;
            xval_index = xval_index + 1;
        end

end

