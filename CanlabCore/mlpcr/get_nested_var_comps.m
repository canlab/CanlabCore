% function vcomps = get_nested_varcomps(X, options)
% 
% Performs variance decomposition of X
%
% Usage ::
%
%   vcomps = get_nested_varcomps(X,'dim1',dim1_id,n_pc1,'dim2',dim2_id,n_pc2, ...)
%
%   X       - n x p matrix to decompose
%
%   'dim1', 'dim2', etc are arbitrary labels.
%
%   dim1_id, dim2_id, etc. are a series of numeric labels which identify
%   random effects blocks. Must be n x 1.
%
%   n_pc1, n_pc2, etc. are not used, but including them means you can use
%   the same arguments for mlpca calls as get_nested_varcomps calls, which
%   is convenient.
%
%   Blocks will be nested in the order in which it's provided. The top most
%   level is all fixed effects, so the 'dimi_id' label should be a unit
%   vector, all zeros, all ones, etc.
%
% Notes ::
%
%   Uses similar approach as mlpca uses internally, with minor mods, etc. 
%   no sqrt(n) rescaling. 
%   
%   This fxn can be sped up by using centering matrices instead of 
%   iterative demeaning
%
function varComp = get_nested_var_comps(X,varargin)
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
    
    % this line differs from mlpca, mlpca uses what's commented out instead
    %X_wi{1} = X - repmat(mean(X,1),size(X,1),1);  % demean within voxel
    X_wi{1} = X;

    [n_obs, vx] = size(X);
    clear X;

    %This is iterative demeaning
    for i = 2:n_lvls
        [X_wi{i}, X_bt{i}] = get_X_bt_X_wi_unweighted(X_wi{i-1},'levelId',block{i});
        X_wi{i-1} = [];
    end
    
    varComp = [X_bt(2:end); X_wi{end}];
    
    % using centering matrices eliminates the need for this re-expansion
    % step, which is time consuming
    %for i = fliplr(1:n_lvls-1)
    %    ids = unique(block{i+1});
    %    tmpComp = zeros(length(block{i+1}),size(varComp{i},2));
    %    for j = 1:length(ids)
    %        this_id = ids(j);
    %        this_idx = find(block{i+1} == this_id);
    %        tmpComp(this_idx,:) = repmat(varComp{i}(this_id,:),length(this_idx),1);
    %    end
    %    varComp{i} = tmpComp;
    %end
end

function [X_wi, X_bt] = get_X_bt_X_wi_unweighted(X, varargin)
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
        
    c = get_cntrng_mat(labels);
    X_wi = c*X;
    X_bt = X - X_wi;
end

% this is the original version of the above, but it's slower than using
% centering matrices by almost a factor of 2, and that's not even counting
% the reexpansion step it subsequently requires
% function [X_wi, X_bt] = get_X_bt_X_wi_unweighted(X, varargin)
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
%         X_bt(i,:) = mean(subjX);
%     end
% end