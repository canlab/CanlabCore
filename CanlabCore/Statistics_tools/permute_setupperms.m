function permindx = permute_setupperms(n,nperms)
% Set up permutations x observations matrix of observation indices for
% permutation test
%
% :Usage:
% ::
%
%    % approximate test with nperms obs.
%    permindx = permute_setupperms(n,nperms)
%    permindx = permute_setupperms(n,[])
%
%    % exact test (all permutations) of n observations
%    permindx = permute_setupperms(n)
%
% ..
%    Uses code from SnPM3b by Tom Nichols
%    adapted for this function by Tor Wager
% ..


  if nargin < 2 || isempty(nperms) || nperms == 0
    OUT.meth = 'exact';
    nperms = gamma(n+1);
    use_approximate = 0;
  else
    nperms = min(gamma(n+1),nperms);
    use_approximate = 1;

    if nperms > 50000
      warning('More than 50,000 permutations! May cause time/memory problems');
    end
  end

				% interactive
				% % nperms = gamma(n+1);
% % use_approximate = input(sprintf('%d Perms. Enter return for exact test, or num. of perms in random subset: ',nperms));
% % if isempty(use_approximate)
% %     use_approximate = 0;
% % else
% %     nperms = min(use_approximate,nperms);
% %     use_approximate = 1;
% % end

  if use_approximate
      % taken from snpm3b code, by Tom Nichols
      %-Approximate test :
      % Build up random subset of all (within n) permutations
      %===============================================================
    rand('seed',sum(100*clock))	%-Initialise random number generator
    permindx      = zeros(nperms,n);
    permindx(1,:) = 1+rem(0:(n-1), n);

    for i = 2:nperms+1
      permindx(i,:) = randperm(n);
    end
    b = unique(permindx, 'rows');
    if(size(b,1) ~= size(permindx, 1))
      while(size(b,1) ~= size(permindx, 1))
        permindx = b;
        for i = (size(permindx,1)+1):nperms+1
          permindx(i,:) = randperm(n);
        end
        b = unique(permindx, 'rows');
      end
      permindx = [b(1,:); b(randperm(size(b,1)-1)+1,:)];
    end
    
		       % null permutation is not included in perm indx
    permindx = permindx(2:end,:);
  else
      % taken from snpm3b code, by Tom Nichols
      %
      %-Full permutation test :
      % Build up exhaustive matrix of permutations
      %===============================================================
      %-Compute permutations for a single exchangability block
      %---------------------------------------------------------------
      %-Initialise Xblkpermindx & remaining numbers
    Xblkpermindx = [];
    lef = [1:n]';
%-Loop through numbers left to add to permutations, accumulating permindx
    for i = n:-1:1
				%-Expand Xblkpermindx & lef
      tmp = round(exp(gammaln(n+1)-gammaln(i+1)));
      Exp = meshgrid(1:tmp,1:i); Exp = Exp(:)';
      if ~isempty(Xblkpermindx), Xblkpermindx = Xblkpermindx(:,Exp); end
      lef = lef(:,Exp);
				%-Work out sampling for lef
      tmp1 = round(exp(gammaln(n+1)-gammaln(i+1)));
      tmp2 = round(exp(gammaln(n+1)-gammaln(i)));
      sam = 1+rem(0:i*tmp1-1,i) + ([1:tmp2]-1)*i;
			      %-Add samplings from lef to Xblkpermindx
      Xblkpermindx   = [Xblkpermindx; lef(sam)];
		      %-Delete sampled items from lef & condition size
      lef(sam) = [];
      tmp = round(exp(gammaln(n+1)-gammaln((i-1)+1)));
      lef = reshape(lef,(i-1),tmp);
		%NB:gamma(n+1)/gamma((i-1)+1) == size(Xblkpermindx,2);
    end
    clear lef Exp sam i
				%-Reorient so permutations are in rows
    permindx = Xblkpermindx';
  end

			     %-Check, condition and randomise permindx
%-----------------------------------------------------------------------
%-Check permindxs sum within Xblks to sum to 1
  if ~all(all(sum(permindx,2) == (n+1)*n/2 ))
    error('Invalid permindx computed!'), end

	%-Convert to full permutations from permutations within blocks
  nperms = size(permindx,1);

%-Randomise order of permindxs (except first) to allow interim analysis
  rand('seed',sum(100*clock))	%-Initialise random number generator
  permindx=[permindx(1,:);permindx(randperm(nperms-1)+1,:)];

			 %-Check first permutation is null permutation
  if ~use_approximate
    if ~all(permindx(1,:)==[1:n])
      error('permindx(1,:)~=[1:n]');
    end
  end

end
