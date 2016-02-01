function out = correl_compare_dep(y1,y2,varargin)
% Compare dependent correlations between pairs of vectors in y1 and y2.
%
% :Usage:
% ::
%
%     out = correl_compare_dep(y1,y2,['alpha',myalpha],['rank'],['table'])
%
% Each of y1 and y2 would contain at least 2 columns, which would be
% correlated and saved in r1 and r2 matrices in output.
% Then, the r1 and r2 matrices are subtracted, and P-values are
% returned for the differences.
%   - In the simplest case, y1 would contain vectors [a b] and y2 would
%     contain vectors [a c].  tests are provided on the a-b vs. a-c
%     difference in correlations.
%
% Repeats dep. correl. analysis for each pair of columns
%
% Returns results in correlation matrix form, where number of rows and
% cols. are the number of pairs [y1(:,i) y2(:,i)]
%
% myalpha is 2-tailed alpha value; p-values are 2-tailed
% FDR correction is at .01, 2-tailed
%
% Based on Steiger, 1980, tests for comparing dependent correlations.
%
% :Examples:
% ::
%
%    for i = 1:length(cl), y1(:,i) = cl.CONTRAST.data(:,2); y2(:,i) = cl.CONTRAST.data(:,1); end
%    for i = 1:length(cl), y1(:,i) = cl(i).CONTRAST.data(:,2); y2(:,i) = cl(i).CONTRAST.data(:,1); end
%
%    % y1 is matrix of obs x data vectors for condition 1
%    % y2 is matrix of obs x data vectors for condition 2
%    out = correl_compare_dep(y1,y2)
%
%    figure('Color','w');nmdsfig(c.GroupSpace,c.ClusterSolution.classes, ...
%    c.names,out.sig,1,{'Pos' 'Neg'});
%    nmdsfig_legend(c.ClusterSolution.X,c.r)
%
% :Examples:
% ::
%
%    c_compare = correl_compare_dep(y1avg,y2avg,'alpha',.05,'rank','table','names',c.APPLY_CLUSTER.names);
%
%    out = correl_compare_dep([ypred pain],[ypred temp], 'alpha', .06, 'table', 'names', {'biomarker resp' 'pain or temp'});


    myalpha = .05;
    dorankdata = 0;
    dotable = 0;
    names = [];

    for i = 1:length(varargin)
        if isstr(varargin{i})
            switch varargin{i}
                % functional commands
                case 'alpha', myalpha = varargin{i+1};
                case 'rank', dorankdata = 1;
                case 'table', dotable = 1;
                case 'names', names = varargin{i+1};
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end


    [N,npairs] = size(y1);


    [rows,cols,ncorr] = corrcoef_indices(npairs);

    if nargin < 4, dorankdata = 0; end
    if dorankdata
        str = sprintf('Ranking data: Nonparametric correlations'); fprintf(1,str);

        for i = 1:npairs
            y1(:,i) = rankdata(y1(:,i));
            y2(:,i) = rankdata(y2(:,i));
        end
    else
        str = sprintf('Assuming continuous data (no ranks).'); fprintf(1,str);
    end

    erase_string(str);

    i = 1; j = 2; k = 3; h = 4;

    str = sprintf('Computing differences among correlations %04d',0); fprintf(1,str);

    diffr = zeros(ncorr,1);
    Zstar2 = zeros(ncorr,1);
    rr = zeros(ncorr,2);

    for cc = 1:ncorr

        fprintf(1,'\b\b\b\b%04d',cc);

        % get full correlation matrix for the 4 variables involved
        dat = [y1(:,[rows(cc) cols(cc)]) y2(:,[rows(cc) cols(cc)])];

        % get differences between correlation z-values (estimate)
        [diffr(cc),r,rr(cc,:),diffrz,z] = correlation_diffs(dat,i,j,k,h);

        s = covcorr(r,i,j,k,h);         % covariance of corr coeffs

        Zstar2(cc) = diffrz * sqrt( (N-3) ./ (2-(2*s)) );

        % bootstrap
        %vals = bootstrp(5000,@correlation_diffs,dat,i,j,k,h);
        %Zboot(cc) = mean(vals) ./ std(vals);

        %if Zstar2(cc) > 2, keyboard, end
    end

    erase_string(str);


    pvec = 2 * (1 - normcdf(abs(Zstar2)));

    r1 = reconstruct(rr(:,1),npairs,ncorr,rows,cols);
    r2 = reconstruct(rr(:,2),npairs,ncorr,rows,cols);

    Z = reconstruct(Zstar2,npairs,ncorr,rows,cols);
    p = reconstruct(pvec,npairs,ncorr,rows,cols);

    dat = [rr diffr Zstar2 pvec];
    dat = [rows cols dat];
    dat = dat(pvec <= myalpha,:);

    diffr = reconstruct(diffr,npairs,ncorr,rows,cols);

    sig = (p <= myalpha - eye(size(p))) .* sign(diffr);

    % FDR corrected
    pthr = FDR(pvec,.05);
    if isempty(pthr), pthr = 0; end

    sigfdr = (p <= pthr) .* sign(diffr);

    out = struct('alpha',myalpha,'r1',r1,'r2',r2,'diffr',diffr, ...
        'Z',Z,'p',p,'sig',sig,'pthr',pthr,'sigfdr',sigfdr,'sigstats',dat);

    % output table...
    if dotable
        if isempty(names)
            disp(['No names entered; try ''names'' keyword to add them.']);
            for i = 1:npairs, names{i} = ['R' num2str(i)]; end
        end

        out.names = names;
        
        disp(['Uncorrected, p < ' num2str(myalpha)])
        maketable(out,'sig',names);
        disp('')

        disp(['FDR corrected, p < ' num2str(myalpha)])
        maketable(out,'sigfdr',names);
        disp('')
    end



    return


function [diffr,r,rr,diffrz,z] = correlation_diffs(dat,i,j,k,h)
    % get differences between correlation z-values
    r = corrcoef(dat);

    rr = [r(i,j) r(k,h)];                % correls to be compared
    diffr = -diff(rr);              % negative sign means we get (1) - (2)

    if nargout > 3
        z = .5 .* log( (1+rr) ./ (1-rr) );
        diffrz = -diff(z);               % diff btwn correl z values
    end

    return



function [rows,cols,ncorr] = corrcoef_indices(npairs)
    % upper triangle only
    tmp = triu(ones(npairs));
    tmp = tmp - eye(npairs);
    [rows,cols] = find(tmp);
    ncorr = length(rows);
    return




function s = covcorr(R,i,j,k,h)

    % pool correl. coeff for more stable var est.
    % as they are equal under Ho.  Steiger, 1980, eq. 14 for Z*
    pooledr = mean([R(i,j) R(k,h)]);
    s = pearsonf(R,i,j,k,h) ./  ( corvar(pooledr) );

    return


function cv = corvar(r)
    % variance of a correlation coefficient
    % simplified form (special case) of pearsonf
    cv = (1 - r^2) ^ 2;
    return


    % % function cv = co(R,x,y)
    % % % x and y are 2-vector indices of matrix R to compute covariance for
    % % % cv is covariance
    % % cv = pearsonf(R,x(1),x(2),y(1),y(2));
    % % return


function pf = pearsonf(R,i,j,k,h)

    % Pearson-Filon: covariance (or var) of element i,j with element k,h
    % depending on correlation values

    % for variance of 1 correl, works as well, and reduces to:
    % (1 - rr(1)^2)^2

    pf = (1/2) .* R(i,j) .* R(k,h) .* ...
        ( R(i,k).^2 + R(i,h).^2 + R(j,k).^2 + R(j,h).^2 ) + ...
        R(i,k) .* R(j,h) + R(i,h) .* R(j,k) - ...
        R(i,j) .* ( R(j,k) .* R(j,h) + R(i,k) .* R(i,h) ) - ...
        R(k,h) .* ( R(j,k) .* R(i,k) + R(j,h) .* R(i,h) );


    return




    % function pf = pearsonf2(r,j,k,h,m)
    %
    % pf = .5 * ( (r(j,h) - r(j,k)*r(k,h)) * (r(k,m) - r(k,h)*r(h,m))  + ...
    %     (r(j,m) - r(j,h)*r(h,m)) * (r(k,h) - r(k,j)*r(j,h)) + ...
    %     (r(j,h) - r(j,m)*r(m,h)) * (r(k,m) - r(k,j)*r(j,m)) + ...
    %     (r(j,m) - r(j,k)*r(k,m)) * (r(k,h) - r(k,m)*r(m,h))    );
    %
    % return

function valmat = reconstruct(vals,npairs,ncorr,rows,cols)

    valmat = zeros(npairs);
    for i = 1:ncorr
        valmat(rows(i),cols(i)) = vals(i);
    end
    valmat = valmat + valmat';

    return



function erase_string(str1)
    fprintf(1,repmat('\b',1,length(str1))); % erase string
    return


function maketable(c_compare,whfield,names)
    
    str = {'-' '+'};
    [rows,cols] = find(triu(c_compare.(whfield)));
    if isempty(rows)
        disp('No significant results.')
    else
        fprintf(1,'Name1\tName2\trow\tcol.\t+ or -\tr1\tr2\tZ\tp\n');
        for i = 1:length(rows)
            fprintf(1,'%s\t%s\t%3.0f\t%3.0f\t%s\t%3.3f\t%3.3f\t%3.2f\t%3.4f\n', ...
                names{rows(i)}, names{cols(i)},rows(i),cols(i), ...                 
                str{1.5 + (.5.*c_compare.(whfield)(rows(i),cols(i)))}, ...          % + or - sign
                c_compare.r1(rows(i),cols(i)),c_compare.r2(rows(i),cols(i)), ...   % correlations
                c_compare.Z(rows(i),cols(i)),c_compare.p(rows(i),cols(i)));        % Z and p
        end
    end
    fprintf(1,'\n');
    return
