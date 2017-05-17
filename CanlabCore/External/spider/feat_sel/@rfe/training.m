function [dat,a] =  training(a,dat)

  dat = set_name(dat,[get_name(dat) ' -> ' get_name(a)]);
  [numEx,vDim,oDim]= get_dim(dat);
  fprintf('training %s...\n', get_name(a))
  
  gi = a.feature_group;
  if isempty(gi), gi = 1:vDim; end
  
  if length(unique(gi)) == vDim
	  fstr = 'feature';
  else
	  fstr = 'feature group';
  end
  vDim = length(unique(gi)); % vDim is now number of groups---in default case, also number of features
  
  force = a.force;
  if ~isempty(force)
  	if ~isequal(unique(gi(:)), [1:vDim]'), error('group ids must be contiguous valid indices if the ''force'' option is set'),end
	if length(force(:)) ~= vDim, error('''force'' property, if supplied, must have as many entries as there are distinct features or feature groups'), end
  end

  feat=a.feat; % number of groups desired 
  rank=[];
  
  untrained = a.child;
  keep = [];
    
  %% if we use the previous trainings, .. (faster to train)
  if a.algorithm.use_prev_train==1 & a.algorithm.trained==1,
      sel = find(ismember(gi, a.rank(1:a.feat)));
	  dat = get(dat,[],sel);
  else
    while 1  
        [res,a.child]=train(untrained,dat);
		[a, testerrstr] = evaluate(a, vDim, keep);
		if vDim <= feat, break, end
		
        w = get_w(a.child);
        gid = unique(gi);
        Ct = [];
        for i = 1:length(gid)
          sel = (gi==gid(i));
% choose elimination criterion:

%          Ct(i) = sum(abs(w(sel)));
          Ct(i) = sum(w(sel).^2);
        end
        % Ct is the criterion value from each group or feature
        % gid is the corresponding group number

		if ~isempty(force)
			pval = force(gid);
			pind = find(pval > 0);
			Ct(pind) = max(Ct) + pval(pind); % positive 'force' values mean that these features/groups will be eliminated last (in order of increasing force value) 
			pind = find(pval < 0);
			Ct(pind) = min(Ct) + pval(pind); % negative 'force' values mean that these features/groups will be eliminated first (in order of increasing force value)
		end
			        
        [val,index]= sort(Ct); % least influential first
        ordered_gid = gid(index);        

        if vDim>a.speed
            nToRemove = floor(vDim/2);    %% choose number of groups to take 
        else
            nToRemove = 1;             %% slow mode- remove one group 
        end
        if vDim-nToRemove < feat, nToRemove = vDim-feat; end
        
        keep = ~ismember(gi, ordered_gid(1:nToRemove)); % which features belong to the nToRemove lowest-ranked groups
        rank = [ordered_gid(1:nToRemove) rank]; % group labels come out ordered according to rank
        dat=get(dat,[],find(keep));     %% cut features
		gi(~keep) = [];
        vDim=vDim-nToRemove;
        
		if a.algorithm.verbosity
			if nToRemove == 1
				str = sprintf('removed %s %d', fstr, ordered_gid(1));
			elseif nToRemove < 8
				str = sprintf('removed %ss [%s]', fstr, num2str(ordered_gid(1:nToRemove)));
			else
				str = sprintf('removed %d %ss', nToRemove, fstr);
			end
       		if ~isempty(force)
				nforced = sum(force(ordered_gid(1:nToRemove))~=0);
				if nforced
					if nToRemove == nforced
						str = sprintf('%s (forced)', str);
					else
						str = sprintf('%s (%d forced)', str, nforced);
					end
				end
			end
			str = sprintf('%s => feat = %d%s', str, vDim, testerrstr);
			fprintf('%s\n', str)
		end
    end 
	
    rank=[unique(gi) rank];
	a.feat=feat;
    a.rank=rank;
  end; 
  
  if a.output_rank==0  %% output selected features, not label estimates
  	dat = res;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [a, str] = evaluate(a, ngroups, keep)
  
  str = '';
  if ~isfield(struct(a), 'test_on_the_fly'), return, end
  if isempty(a.test_on_the_fly) | ~isstruct(a.test_on_the_fly), return, end
  if isempty(a.test_on_the_fly.data), return, end

  if ~isfield(a.test_on_the_fly, 'result'), a.test_on_the_fly.result = []; end
  if ~isfield(a.test_on_the_fly, 'nfeatures'), a.test_on_the_fly.nfeatures = []; end
  if ~isfield(a.test_on_the_fly, 'ngroups'), a.test_on_the_fly.ngroups = []; end
  if ~isfield(a.test_on_the_fly, 'loss'), a.test_on_the_fly.loss = []; end

  if ~isa(a.test_on_the_fly.result, 'group'), a.test_on_the_fly.result = group; end
  if isempty(a.test_on_the_fly.loss), a.test_on_the_fly.loss = loss('class_loss'); end

  if nargin < 3, keep = []; end
  if ~isempty(keep)
	a.test_on_the_fly.data = get(a.test_on_the_fly.data, [], find(keep));
  end
  [ans, nfeatures] = get_dim(a.test_on_the_fly.data);  
  n = 1 + length(a.test_on_the_fly.result.child);
  res = test(a.child, a.test_on_the_fly.data, a.test_on_the_fly.loss);
  res.X = ngroups;
  str = sprintf(' (cv error: %.3g)', res.Y);
  a.test_on_the_fly.result.child{n} = res;
  a.test_on_the_fly.nfeatures(n) = nfeatures;
  a.test_on_the_fly.ngroups(n) = ngroups;
  
