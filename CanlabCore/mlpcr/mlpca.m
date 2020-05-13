% function [CM, SCW, labels] = mlpca(X, options)
% 
% Performs multilevel PCA
%
% Usage ::
%
%   [CM, SCW, labels] = mlpca(X,'dim1',dim1_id,n_pc1,'dim2',dim2_id,n_pc2, ...)
%
%   CM and SCW are cell arrays. CM{i} are the weights of PC i, and SCW{i}
%       are the loadings of PC i.
%
%   'dim1', 'dim2', etc are arbitrary labels.
%
%   dim1_id, dim2_id, etc. are a series of numeric labels which identify
%   random effects blocks.
%
%   Blocks will be nested in the order in which it's provided. The top most
%   level is all fixed effects, so the 'dimi_id' label should be a unit
%   vector, all zeros, all ones, etc.
%
% Dependencies:
%   get_cntrng_mat          (required by mlpca)
%
% References ::
%
%   Timmerman, M. (2006) Multilevel component analysis. British Journal of 
%       Mathematical and Statistical Psychology.

function [CM, SCW, dLabels] = mlpca(X, varargin)

    dLabels = cell(1,1);
    id = dLabels;
    n_d = dLabels;
    
    %sanity check on dimensions still needs to be implemented
    idx = 1;
    for i = 1:length(varargin)
        if ischar(varargin{i})
            dLabels{idx} = varargin{i}; % the name you're calling this level
            id{idx} = varargin{i+1}(:); % group labels
            n_d{idx} = varargin{i+2}; % number of PC dimensions to retain
            idx = idx+1;
        end
    end
    n_lvls = length(dLabels);
    
    if length(unique(id{1})) > 1
        warning('Random effects blocks were specified for model top level. Forcing unit block and ignoring group membership.');
        n = length(id{1});
        id{i} = ones(n,1);
    end

    n = length(id{1});
    for i = 1:n_lvls
        if n ~= length(id{i})
            fprintf('Random effects label vector for level %d is %d, but label level %d has %d.\n',i,length(id{i}),i-1,length(id{i-1}));
            error('Vectors of random effect block labels must all have the same length.');
        end
    end
    
    % sorting labels
    labels = cell2mat(id);
	original_order = 1:size(X,1);
    for i = fliplr(2:n_lvls)
        [labels,new_ord] = sortrows(labels,i);
        X = X(new_ord,:);
		original_order = original_order(new_ord);
        for j = 1:n_lvls
            id{j} = id{j}(new_ord);
        end
    end
    
    block = cell(n_lvls,1);
    block_n = block;
    bid = block;
    for i = 1:n_lvls
        [~, block0] = unique(labels(:,1:i),'rows');
        block_n{i} = diff(block0);
        block_n{i}(end + 1) = length(id{i}) - block0(end) + 1;
        for j = 1:length(block_n{i})
            block{i} = [block{i}(:); j*ones(block_n{i}(j),1)];
        end
        bid{i} = unique(block{i});
    end
    
    X_wi = cell(n_lvls,1);
    X_bt = X_wi;
    X_wi{1} = X - repmat(mean(X,1),size(X,1),1); 

    [n_obs, vx] = size(X);
    clear X;

    %X_bt{1} = X_wi{1};
    for i = 2:n_lvls
        [X_wi{i}, X_bt{i}] = get_X_bt_X_wi(X_wi{i-1},'levelId',block{i});
        X_wi{i-1} = [];
    end
    
    % don't center, we've already done that above. Any residual mean is
    % due to machine precision issues. Should be ~10e-5 if type=single
    CM = cell(n_lvls,1);
    SCW = cell(n_lvls,1);
    
    if n_d{n_lvls} > 0
        [CM{n_lvls}, SCW{n_lvls}] = ...
            pca(X_wi{n_lvls}, 'Centered', false,'NumComponents', n_d{n_lvls}, 'Economy', true);
        X_wi{n_lvls} = [];
    end
    
    for i = 2:n_lvls
        if n_d{i-1} > 0
            [CM{i-1}, pre_SCW] = ...
                pca(X_bt{i}, 'Centered', false, 'NumComponents', n_d{i-1}, 'Economy', true);
            X_bt{i} = [];

            SCW{i-1} = zeros(n_obs,size(pre_SCW,2));

            for j = 1:length(block_n{i})
                idx = find(bid{i}(j) == block{i});
                SCW{i-1}(idx,:) = repmat(1/sqrt(block_n{i}(j))*pre_SCW(j,:),length(idx),1);
            end
        end
    end

    [~, b] = sort(original_order);
    for i = 1:n_lvls
        if n_d{i} > 0
            SCW{i} = SCW{i}(b,:);
        end
    end
end

% This function was purpose built for use with mlpca. It rescales X_bt
% components, weighting each by the sqrt of the number of observations at a
% given level. So for instance, if the level is "subject" and there are 100 
% within subject observations, each row of X_bt will have been weighted by 
% sqrt(100) so that it is accurately handled by mlpca. If this function is 
% invoked in other contexts, be aware that a manual rescaling of the output 
% may be necessary. See Timmerman, et al. for details
function [X_wi, X_bt] = get_X_bt_X_wi(X, varargin)
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                case 'levelId'
                    labels = varargin{i+1};
            end
        end
    end
    
    meanDat = mean(X);
    X = X - repmat(meanDat,size(X,1),1);
    
    [sid, sid0] = unique(labels);
    
    X_wi = zeros(size(X));
    X_bt = zeros(length(sid),size(X,2));
    [c,sqrtn] = get_cntrng_mat(labels);
    X_wi = c*X;
    X_bt = sqrtn(sid0).*(X(sid0,:) - X_wi(sid0,:));
end

% this is an old version of the above function that doesn't rely on
% centering matrices. The centering matrix version should run faster
% function [X_wi, X_bt] = get_X_bt_X_wi(X, varargin)
%     for i = 1:length(varargin)
%         if ischar(varargin{i})
%             switch varargin{i}
%                 case 'levelId'
%                     labels = varargin{i+1};
%             end
%         end
%     end
%     
%     meanDat = mean(X);
%     X = X - repmat(meanDat,size(X,1),1);
%     
%     [sid, sid0] = unique(sort(labels));
%     sid_n = diff(sid0);
%     sid_n(end + 1) = length(labels) - sid0(end) + 1;
%     
%     X_wi = zeros(size(X));
%     X_bt = zeros(length(sid),size(X,2));
%     for i = 1:length(sid)
%         sidx = find(sid(i) == labels);
%         
%         subjX = X(sidx,:);
%         X_wi(sidx,:) = subjX - repmat(mean(subjX),sid_n(i),1);
%         
%         % rescale by level measure count
%         X_bt(i,:) = sqrt(sid_n(i))*mean(subjX);
%     end
% end