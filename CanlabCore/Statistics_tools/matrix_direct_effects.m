function [direct_mtx, mediated_mtx, mediators] = matrix_direct_effects(sig, data)
% :Usage:
% ::
%
%     [direct_mtx, mediated_mtx, mediators] = matrix_direct_effects(sig, data);
%
% Take an n x n matrix of significant correlations and the data for each
% variable. Test whether each link is completely mediated by another
% variable. Prune the correlation significance matrix, also noting which
% variables are mediators.
%
% Example:
% Take output from cluster_nmds (c structure) and prune to return direct
% links only:
% ::
%
%    [c.direct_mtx, c.mediated_mtx, c.mediators] = matrix_direct_effects(c.STATS.sigmat, c.dat);
%    nmdsfig_tools('removelines');
%    nmdsfig_tools('drawlines',c.GroupSpace, c.direct_mtx);
%
% ..
%    Tor Wager, Feb 2007
% ..

[rows,cols] = find(tril(sig));

n = size(sig,1);
indx = 1:n;

n_sig_pairs = length(rows);

% outputs
direct = ones(1, n_sig_pairs);
mediated = zeros(1, n_sig_pairs);

mediators = zeros(1, n);

t1 = clock;
fprintf(1, '\nTesting for mediators of %d significant bivariate relationships: %6d', n_sig_pairs, 1); 

for i = 1:n_sig_pairs
    
    fprintf(1, '\b\b\b\b\b\b%6d', i);
    
    this_pair = [rows(i) cols(i)];
    
    % x and y, the relationship potentially mediated by another
    x = data(:, rows(i));
    y = data(:, cols(i));
    
    % potential mediators
    indx_other_regions = indx(~(indx == rows(i) | indx == cols(i)));
    m_candidates = data(:, indx_other_regions);
    
    % search for mediators; see mediation.m
    med_results = mediation_search('M', x, y, m_candidates, 'noverbose');
    
    sigfx = med_results.pvals < .05;
    sig_mediation = sigfx(:,1) & sigfx(:,2) & sigfx(:,5);

    % preserve signs for later use; store as 1 or -1 in matrix
    mediation_sign = sign(med_results.paths(:,5));
    direct_sign = sign(med_results.paths(1,4));     % should all be same!
   
    complete_mediation = sig_mediation & ~sigfx(:,3);
    
    % collect final output
    if any(complete_mediation)
        direct(i) = 0;
        
        mediators(complete_mediation) = mediation_sign(complete_mediation);
    else
        direct(i) = direct_sign(1);
    end
    
    if any(sig_mediation)
        mediated(i) = 1;
    end
    
    
end

fprintf(1, '\nDone in %d s\n', etime(clock,t1))

% collect output matrices
mediated_mtx = zeros(n);
direct_mtx = zeros(n);

for i = 1:n_sig_pairs
    mediated_mtx(rows(i), cols(i)) = mediated(i);
    direct_mtx(rows(i), cols(i)) = direct(i);
    
    mediated_mtx(cols(i), rows(i)) = mediated(i);
    direct_mtx(cols(i), rows(i)) = direct(i);
end


end

