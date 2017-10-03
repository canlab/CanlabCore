function OUT = correl_compare_permute(meth,dat,nperms,condition)
% General function for comparing correlation matrices on two different
% sets of observations
%
% :Usage:
% ::
%     correl_compare_permute(meth,dat,nperms,condition)
%
%     dat is obs x variables, and an n x n association matrix will be computed
%     on the n(n - 1) pairs of columns using one of the methods in
%     correlation.m and specified by meth
%
% :Inputs:
%
%   **meth:**
%        String indicating the method for computing the correlation
%        coefficient (e.g., 'taub') see correlation.m for details
%
%   **dat:**
%        Matrix of observations (n instances by p variables)
%
%   **nperms:**
%        Number of permutations for inferential tests 
%
%   **condition:**
%        Vector indicating condition membership for each instance
%        (currently only works for two conditions)
%
% :Outputs:
%
%   **OUT:**
%        Output stats structure
%
%
%
% :See also:
%   - correl_compare_indep_inputr, correl_compare_indep, correlation*
% ..
%    this is a work in progress, and documentation is incomplete.
%    tor wager, Jan 15, 2007
% ..

  if isempty(meth), meth = 'taub'; end

				% set up permutations
  fprintf(1,'Setting up permutations. \n');
  [n,m] = size(dat);

  permindx = permute_setupperms(n,nperms);


				% correct permutation
  fprintf(1,'Correct perm. \n');
  [cdiff, cdiffmax,c1,c2] = get_corrmat(meth,dat,condition);


				% init perms
  cdiffprm = zeros(m,m,nperms);
  cdiffmaxp = zeros(nperms,1);
  cdiff_matrix = zeros(nperms,m*(m-1)/2);
  fprintf(1,'Permuting:\n');
  parfor i = 1:nperms
    % permute rows to shuffle condition assignment labels
    datp = dat(permindx(i,:),:);
    [cdiffprm(:,:,i), cdiffmaxp(i),null,null,cdiff_matrix(i,:)] = get_corrmat(meth,datp,condition);
  end

  fprintf(1,' Done.\n');

  OUT = struct('meth',meth,'c1',c1,'c2',c2,'cdiff',cdiff,'cdiffprm',cdiffprm,'cdiff_matrix',cdiff_matrix,'cdiffmaxp',cdiffmaxp);

				% Results and thresholding
  OUT.diff_corrected05thr = prctile(OUT.cdiffmaxp,95);

  [rows,cols,ncorr] = corrcoef_indices(m);
  for cc = 1:ncorr
    
				% do threshold on each pair separately
    this_cor = squeeze(OUT.cdiffprm(rows(cc),cols(cc),:)); % null-hypothesis
				% one-tailed
    pv = sum(this_cor <= OUT.cdiff(rows(cc),cols(cc))) ./ nperms;
    pval(cc) = min(pv,1-pv);
    
  end

  OUT.pdiff = reconstruct(pval,m,ncorr,rows,cols);
  OUT.pdiff = OUT.pdiff + eye(m);

end




function [cdiff, cdiffmax,c1,c2,cdiff_matrix] = get_corrmat(meth,dat,condition)

  c1 = correlation('taub',dat(condition == 1,:));
  c2 = correlation('taub',dat(condition == 2,:));
  cdiff = c1 - c2;
  % works for any diff between correls, two tailed
  cdiff_matrix = squareform(cdiff);
  cdiffmax = max(abs(cdiff_matrix));
  % works for any correlation with ones on the diagonal
  % two-tailed
  %cdiffmax = max(abs(squareform(cdiff - 1) + 1));
end


function [rows,cols,ncorr] = corrcoef_indices(npairs)
    % upper triangle only
    tmp = triu(ones(npairs));
    tmp = tmp - eye(npairs);
    [rows,cols] = find(tmp);
    ncorr = length(rows);
end

    

function valmat = reconstruct(vals,corrmatsize,ncorr,rows,cols)

    valmat = zeros(corrmatsize);
    for i = 1:ncorr
        valmat(rows(i),cols(i)) = vals(i);
    end
    valmat = valmat + valmat';

end

    
