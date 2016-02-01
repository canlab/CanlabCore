function out = correl_compare_dep_permtest(y1,y2,varargin)
% Compare dependent correlations between pairs of vectors in y1 and y2
%
% :Usage:
% ::
%
%     out = correl_compare_dep_permtest(y1,y2,['alpha',myalpha],['rank'],['table'])
%
% :PERMUTATION TEST: for correl_compare_dep
%
% Repeats dep. correl. analysis for each pair of columns
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
%    y1 is matrix of obs x data vectors for condition 1
%    y2 is matrix of obs x data vectors for condition 2
%    out = correl_compare_dep(y1,y2)
%
%    figure('Color','w');nmdsfig(c.GroupSpace,c.ClusterSolution.classes, ...
%    c.names,out.sig,1,{'Pos' 'Neg'});
%    nmdsfig_legend(c.ClusterSolution.X,c.r)
%
%    % compare correlations on cluster averages
%    c_compare = correl_compare_dep(y1avg,y2avg,'alpha',.05,'rank','table','names',c.APPLY_CLUSTER.names);


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

    % ----------------------------------------
    %%% Correct Permutation %%%
    % ----------------------------------------
    str = sprintf('Computing differences among correlations %04d',0); fprintf(1,str);
    
    [diffr,Zstar2,rr,pvec] = compare_dep_subfcn(ncorr,y1,y2,rows,cols,i,j,k,h,N);

    erase_string(str);
    
    % ----------------------------------------
    %%% All output stuff %%%
    % ----------------------------------------
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
            for ii = 1:npairs, names{ii} = ['R' num2str(ii)]; end
        end

        disp(['Uncorrected, p < ' num2str(myalpha)])
        maketable(out,'sig',names);
        disp('')

        disp(['FDR corrected, p < ' num2str(myalpha)])
        maketable(out,'sigfdr',names);
        disp('')
    end

    % ----------------------------------------
    %%% Permutations %%%
    % ----------------------------------------
    niter = 5000;

    [out.nalpha,out.n01,out.nalpha_posvsneg,out.n01_posvsneg] = do_perms;
    out.crit_n_at_alpha = prctile(out.nalpha,95);
    out.crit_n_at_01 = prctile(out.n01,95);
    out.crit_sumsigned_at_alpha = prctile(out.nalpha_posvsneg,97.5);
    out.crit_sumsigned_at_01 = prctile(out.n01_posvsneg,97.5);
    
    out.numsig_at_alpha = sum(squareform(out.p) <= out.alpha);
    out.numsig_at_01 = sum(squareform(out.p) <= .01);

    out.sumsigned_at_alpha = sum(  (squareform(out.p) <= out.alpha) .* sign(squareform(out.diffr))  );
    out.sumsigned_at_01 = sum(  (squareform(out.p) <= .01) .* sign(squareform(out.diffr))  );
    
    
    if dotable
        fprintf(1,'Alpha: %3.4f, Num. sig: %3.0f, Critical num sig: %3.0f\n',out.alpha,out.numsig_at_alpha,out.crit_n_at_alpha);
        fprintf(1,'Alpha: %3.4f, Num. sig: %3.0f, Critical num sig: %3.0f\n',.01,out.numsig_at_01,out.crit_n_at_01);
            
        fprintf(1,'Alpha: %3.4f, Sum of signed sig.: %3.0f, Critical num 2-tailed: %3.0f\n',out.alpha,out.sumsigned_at_alpha, out.crit_sumsigned_at_alpha);
        fprintf(1,'Alpha: %3.4f, Sum of signed sig.: %3.0f, Critical num 2-tailed: %3.0f\n',.01,out.sumsigned_at_01, out.crit_sumsigned_at_01);
    
    end
    
    
    % nested function
    % get number of significant effects under null hypothesis at alpha and at
    % .01
    function [nalpha,n01,nalpha_posvsneg,n01_posvsneg] = do_perms
        
        nalpha = zeros(1,niter);
        n01 = nalpha;
        nalpha_posvsneg = nalpha;
        n01_posvsneg = nalpha;
        
        fprintf(1,'Getting distribution for number of significant correlations.\n  Running %3.0f permutations: %04d',niter,0);
        
        for ii = 1:niter
            fprintf(1,'\b\b\b\b%04d',ii);
     
            % permute data
            for jj = 2:npairs
                wh = randperm(N);
                y1(:,jj) = y1(wh,jj);
                y2(:,jj) = y2(wh,jj);                
            end
            
            % test correlations
            [diffr,Zstar2,rr,pvec] = compare_dep_subfcn(ncorr,y1,y2,rows,cols,i,j,k,h,N,0);
           
            nalpha(ii) = sum(pvec <= myalpha);
            n01(ii) = sum(pvec <= .01);
            
            nalpha_posvsneg(ii) = sum(  (pvec <= myalpha) .* sign(diffr)  );
            n01_posvsneg(ii) = sum(  (pvec <= .01) .* sign(diffr)  );
            
        end
        fprintf(1,'\n');
        
    end
        
    end   %%% end main function

    
    
    
    
    
    
    function [diffr,Zstar2,rr,pvec] = compare_dep_subfcn(ncorr,y1,y2,rows,cols,i,j,k,h,N,dotext)
    
    diffr = zeros(ncorr,1);     % correlation difference
    Zstar2 = zeros(ncorr,1);    % from Steiger, the Z-score
    rr = zeros(ncorr,2);        % the correl values for y1 and y2

    if ~exist('dotext','var'), dotext = 1; end
    
    for cc = 1:ncorr

        if dotext, fprintf(1,'\b\b\b\b%04d',cc); end

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

    pvec = 2 * (1 - normcdf(abs(Zstar2)));

    end
    
    
    

function [diffr,r,rr,diffrz,z] = correlation_diffs(dat,i,j,k,h)
    % get differences between correlation z-values
    r = corrcoef(dat);

    rr = [r(i,j) r(k,h)];                % correls to be compared
    diffr = -diff(rr);              % negative sign means we get (1) - (2)

    if nargout > 3
        z = .5 .* log( (1+rr) ./ (1-rr) );
        diffrz = -diff(z);               % diff btwn correl z values
    end

    end



function [rows,cols,ncorr] = corrcoef_indices(npairs)
    % upper triangle only
    tmp = triu(ones(npairs));
    tmp = tmp - eye(npairs);
    [rows,cols] = find(tmp);
    ncorr = length(rows);
    end




function s = covcorr(R,i,j,k,h)

    % pool correl. coeff for more stable var est.
    % as they are equal under Ho.  Steiger, 1980, eq. 14 for Z*
    pooledr = mean([R(i,j) R(k,h)]);
    s = pearsonf(R,i,j,k,h) ./  ( corvar(pooledr) );

    end


function cv = corvar(r)
    % variance of a correlation coefficient
    % simplified form (special case) of pearsonf
    cv = (1 - r^2) ^ 2;
    end


    % % function cv = co(R,x,y)
    % % % x and y are 2-vector indices of matrix R to compute covariance for
    % % % cv is covariance
    % % cv = pearsonf(R,x(1),x(2),y(1),y(2));
    % % end


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


    end




    % function pf = pearsonf2(r,j,k,h,m)
    %
    % pf = .5 * ( (r(j,h) - r(j,k)*r(k,h)) * (r(k,m) - r(k,h)*r(h,m))  + ...
    %     (r(j,m) - r(j,h)*r(h,m)) * (r(k,h) - r(k,j)*r(j,h)) + ...
    %     (r(j,h) - r(j,m)*r(m,h)) * (r(k,m) - r(k,j)*r(j,m)) + ...
    %     (r(j,m) - r(j,k)*r(k,m)) * (r(k,h) - r(k,m)*r(m,h))    );
    %
    % end

function valmat = reconstruct(vals,npairs,ncorr,rows,cols)

    valmat = zeros(npairs);
    for i = 1:ncorr
        valmat(rows(i),cols(i)) = vals(i);
    end
    valmat = valmat + valmat';

    end



function erase_string(str1)
    fprintf(1,repmat('\b',1,length(str1))); % erase string
    end


% function maketable(c_compare,whfield,names)
%     [rows,cols] = find(triu(c_compare.(whfield)));
%     if isempty(rows)
%         disp('No significant results.')
%     else
%         fprintf(1,'Name1\tName2\trow\tcol.\t+ or -\tZ\tp\n');
%         for i = 1:length(rows)
%             fprintf(1,'%s\t%s\t%3.0f\t%3.0f\t%3.0f\t%3.2f\t%3.4f\n',names{rows(i)}, ...
%                 names{cols(i)},rows(i),cols(i),c_compare.(whfield)(rows(i),cols(i)),c_compare.Z(rows(i),cols(i)),c_compare.p(rows(i),cols(i)));
%         end
%     end
%     fprintf(1,'\n');
%     end


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
end
    
    
