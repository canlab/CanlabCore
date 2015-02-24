function sd=distosim(ds)
% function sd=distosim(ds)
% Convert a distance or dissimilarity matrix to a similarity matrix


% only 3 dimensions allowed
% --------------------------------------

if ndims(ds)>3
    nd=size(ds);
    ds=ds(:,:,:);
end


% Check Matrix
% --------------------------------------

% dimensionality, and cutoff value to assess symmetry

n = size(ds); m = n(2); n = n(1);
del = 10*eps(class(ds));

if n ~= m, 
    error('Matrix is not square: You must have subjects on the third dimension');
end

% check symmetry:  must be symmetrical
% --------------------------------------

for i = 1:size(ds,3)
    D = ds(:,:,i);
    symm(i) = all(all(D >= 0 & abs(D - D') <= del*max(max(D))));
end

if any(~symm), error('Matrices are not symmetrical or there are negative values: Not valid similarity or dissimilarity matrix.'),end


% whether it's a sim, dissim, or neither
% --------------------------------------

for i = 1:size(ds,3)
    D = ds(:,:,i);
        
    % it's a dissimilarity matrix, zeros on diagonal
    if all(diag(D) < del)
        dissim(i) = 1;
        
        % do the conversion
        D=(1-D.^2);
        sd(:,:,i) = D;
        
    else
        dissim(i) = 0;
    end
   
    % it's a similarity matrix, ones on diagonal, no entries > 1
    if all(abs(diag(D) - 1) < del) & all(all(D < 1+del))
        sim(i) = 1;
        
        % do the conversion (but don't do cause this is distosim)
        %D=(1-D) .^ .5;
        %sd(:,:,i) = D;

    else
        sim(i) = 0;
    end
end

if sum(dissim) == size(ds,3)
    simdis = 'dis';
    disp('converting dis-similarity to similarity')
elseif sum(sim) == size(ds,3)
    simdis = 'sim';
    disp('already a similarity matrix');
    sd = ds;
    %disp('converting similarity to dis-similarity');
else
    dissim
    sim
    error('Not all matrices are consistently either similarity or dissimilarity!');
end

return

