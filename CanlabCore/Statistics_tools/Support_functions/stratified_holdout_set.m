function cvpart = stratified_holdout_set(Y, varargin)
% cvpart = stratified_holdout_set(Y, [nfolds or integer vector for holdout set ID], [optional keywords])
%
% CANlab holdout set identification 
% Simple stratified holdout set, stratifying on level of Y
%
% Use cases:
% -------------------------------------------------------------------
% stratified_holdout_set(Y, [vector of integers length Y])
% Custom holdout set
% Returns folds specified in integers
% fold_indicator is a copy of the input vector of integers
% nfolds is based on unique elements in vector of integers
%
% stratified_holdout_set(Y)
% 5-fold stratified holdout set; see below
%
% stratified_holdout_set(Y, 1)
% Use all obs for train and test; good for bootstrapping weights
%
% stratified_holdout_set(Y, nfolds)
% Returns number of folds specified in nfolds 
% If Y is categorical: Uses cvpartition (matlab) to stratify on class ID
% specified in Y
% If Y is continuous: holdout sets are every nth entry in an ordered
% (sorted) list of Y values
%
% stratified_holdout_set(Y, 'rolling', [h v g])
% Returns rolling HVG block holdout set
% h = window length for test set
% v = window length for buffer
% g = window length for training set
%
% stratified_holdout_set(Y, 'hvblock', [h v])
% Returns rolling HVG block holdout set
% h = window length for test set
% v = window length for buffer
% (all other observations are in training set)
% needs 2 * h + 2 * v observations
%
% Examples:
% 
% 3-fold on random variable:
% cvpart = stratified_holdout_set(rand(100, 1), 3)
%
% Create and visualize rolling hv block
% cvpart = stratified_holdout_set(rand(100, 1), 'hvblock', [5 5]);
% cat(2, cvpart.teIdx{:}) - cat(2, cvpart.trIdx{:})
% testvstrain = cat(2, cvpart.teIdx{:}) - cat(2, cvpart.trIdx{:})
% figure; imagesc(testvstrain);


% ---------------------------------------------------------------------
% Defaults
% ---------------------------------------------------------------------

nfolds = 5;
tsxval_hvblock = 0;
tsxval_rolling = 0;
fold_indicator = ones(size(Y));

for i = 1:length(varargin)
    if ~ischar(varargin{i})
        % number: nfolds or holdout set integer vector (parsed below)
        nfolds = varargin{i};  
        
    else
        % keywords
        switch varargin{i}
                
            case 'rolling'
                tsxval_rolling = 1;
                inval = varargin{i + 1};
                h = inval(1);
                v = inval(2);
                g = inval(3);
                varargin{i} = [];
                varargin{i + 1} = [];
                
            case 'hvblock'
                tsxval_hvblock = 1;
                inval = varargin{i + 1};
                h = inval(1);
                v = inval(2);
                varargin{i} = [];
                varargin{i + 1} = [];
                
            otherwise
                warning(['Unknown input: ' varargin{i}])     
                
        end
    end
end

% opt = statset('crossval');

% ---------------------------------------------------------------------
% stratified partition, or custom holdout set
% ---------------------------------------------------------------------
   
if numel(nfolds) > 1 
    % ---------------------------------------------------------------------
    % a vector; assume integers for holdout sets
    % Custom holdout set
        
    fold_indicator = nfolds;    % nfolds IS the indicator integer vector
    u = unique(fold_indicator);
    nfolds = length(u);
    
    cvpart = struct('NumTestSets', nfolds);
    [trIdx, teIdx] = deal(cell(1, nfolds));
    
    for i = 1:length(u)
        teIdx{i} = fold_indicator == u(i);
        trIdx{i} = ~teIdx{i};
    end
    
elseif nfolds == 1 
    % ---------------------------------------------------------------------
    % special for 1 fold: use all obs for train and test; good for bootstrapping weights
    
    [trIdx, teIdx] = deal(cell(1, nfolds));
    trIdx={ones(size(Y))};
    teIdx={ones(size(Y))};
    cvpart = struct('NumTestSets', length(trIdx)); %LC: 4/3/14: added this as it was missing 
    % fold_indicator = ones(size(Y)); default
    
elseif tsxval_hvblock == 1 
    %special case for timeseries CV using HVBlock : added LC: 11/28/13
    % ---------------------------------------------------------------------
    [trIdx, teIdx] = tscv(length(Y), 'hvblock',[h,v]);
    cvpart = struct('NumTestSets', length(trIdx));
    %fold_indicator = ones(size(Y));
    
elseif tsxval_rolling == 1 
    %special case for timeseries CV using HVBlock : added LC: 12/16/13
    % ---------------------------------------------------------------------
    [trIdx, teIdx] = tscv(length(Y), 'rolling',[h,v,g]);
    cvpart = struct('NumTestSets', length(trIdx));
    
else
    % Continuous or categorical data
    % ---------------------------------------------------------------------
   
    [trIdx, teIdx] = deal(cell(1, nfolds));
    
    % Classification: Stratified holdout set
    % ---------------------------------------------------------------------

    if length(unique(Y)) < length(Y) / 2
        
        %cvpart = cvpartition(Y,'k',nfolds); %changed by LC 2/26/13 to account for multiclass svm
        %cvpart = cvpartition(length(Y),'k',nfolds);
        cvpart = cvpartition(Y,'k',nfolds);   % Y must be entered to use stratification
        
        for i = 1:cvpart.NumTestSets
            
            trIdx{i} = cvpart.training(i);
            teIdx{i} = cvpart.test(i);
            
            fold_indicator(teIdx{i}) = i;
        end
        
        % re-create to avoid conflicts with cvpart special class
        cvpart = struct('NumTestSets', nfolds);
        
    else
        % ---------------------------------------------------------------------

        % Regression:  stratification based on continuous values
        % cvpartition object will not stratify continuous values
        % do our own
        
        cvpart = struct('NumTestSets', nfolds);
        
        [ys, indx] = sort(Y);
        for k = 1:nfolds
            
            wh_holdout = indx(k:nfolds:end);
            if isempty(wh_holdout), error('Holdout set construction error: Dataset too small?'); end
            
            teIdx{k} = false(size(Y));
            teIdx{k}(wh_holdout) = true;
            
            trIdx{k} = ~teIdx{k};
            
            fold_indicator(teIdx{k}) = k;
        end
    end
end

cvpart.teIdx = teIdx;
cvpart.trIdx = trIdx;
cvpart.fold_indicator = fold_indicator;

end % function

