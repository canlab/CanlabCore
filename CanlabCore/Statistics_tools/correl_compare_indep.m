function [r1, r2, rdiff, rdiff_p, rdiff_zscore] = correl_compare_indep(y_g1, y_g2, varargin)
% Compute statistics on the difference between correlation coefficients
% for indepndent samples g1 and g2.
%
% :Usage:
% ::
%
%     [r1, r2, rdiff, rdiff_p, rdiff_zscore] = correl_compare_indep(y_g1, y_g2, varargin)
%
% Correlations are computed over y (n observations x k variables) for each
% group.
%
% Source, Hubert Blalock, Social Statistics, NY: McGraw-Hill, 1972:  406-407.
% From notes on Garson stats website.
%
%
% :Inputs:
%
%   **y_g1:**
%        Matrix of variables for group 1 (n observations x k variables)
%
%   **y_g2:**
%        Matrix of variables for group 2 (n observations x k variables)
%     
% :Optional Inputs:
%
%   **varargin:**
%        String 'noverbose' can be used to prevent output from being
%        printed to command line
%
% :Outputs:
%
%   **r1:**
%        Matrix of correlation coefficients for group 1
%
%   **r2:**
%        Matrix of correlation coefficients for group 2
%
%   **rdiff:**
%        Difference in correlation coefficient matrices (r1 - r2)
%
%   **rdiff_p:**
%        p values testing against no difference between r1 and r2. 
%        Hubert Blalock, Social Statistics, NY: McGraw-Hill, 1972:  406-407.
%
%   **rdiff_zscore:**
%        z scores for the difference between r1 and r2. 
%
% :Examples:
% ::
%
%    %generate two random input matrices (20 subjects x 5 variables)
%    y_g1 = rand(20,5);
%    y_g2 = rand(20,5);
%    [r1, r2, rdiff, rdiff_p, rdiff_zscore] = correl_compare_indep(y_g1,y_g2);
%
%    %supress output
%    [r1, r2, rdiff, rdiff_p, rdiff_zscore] = correl_compare_indep(y_g1,y_g2,'noverbose');
%
%
% :See also:
%   - correl_compare_indep_inputr, correl_compare_permute*
%
% ..
%    Tor Wager, Aug. 2008
%
%    Edit: March 2010-add no verbose option.  input 'noverbose'
%    Bug fix in Fisher's Z transform
% ..

    verbose = 1;
    if length(varargin) > 0
        for i = 1:length(varargin)
            if ischar(varargin{i})
                switch varargin{i}
                    case 'noverbose', verbose = 0;

                    otherwise, warning(['Unknown input string option:' varargin{i}]);
                end
            end
        end
    end

    k = size(y_g1, 2);

    if k ~= size(y_g2, 2), error('Num. variables for each group must match.'); end

    r1 = corrcoef(y_g1);
    r2 = corrcoef(y_g2);

    n1 = size(y_g1, 1);
    n2 = size(y_g1, 1);

    % Fisher's Z scores and stats on each group
    [tmp, r1_sig] = r2z(r1(:), n1);
    r1_z = fisherz(r1(:));
    
    r1_sig = reshape(r1_sig, k, k);
    r1_z = reshape(r1_z, k, k);

    [tmp, r2_sig] = r2z(r2(:), n2);
    r2_z = fisherz(r2(:));
    
    r2_sig = reshape(r2_sig, k, k);
    r2_z = reshape(r2_z, k, k);

    if verbose
        disp('Group 1'); correlation_to_text( r1, r1_sig);
        disp('Group 2'); correlation_to_text( r2, r2_sig);
    end

    % difference and stats

    rdiff = r1 - r2;

    rdiff_z = r1_z - r2_z;
    rdiff_se =  sqrt(1./(n1 - 3) + 1./(n2 - 3));

    rdiff_zscore = rdiff_z ./ rdiff_se;

    rdiff_p = 2 * (1 - normcdf(abs(rdiff_zscore)));

    if verbose
        disp('Difference');
        print_matrix(rdiff);

        disp('Difference p-values');
        print_matrix(rdiff_p);
    end

end


