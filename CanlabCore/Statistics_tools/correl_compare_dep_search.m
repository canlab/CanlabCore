function correl_compare_dep_search(seed1,seed2,y1,y2,varargin)
% Compare dependent correlations between seed1<->y1 and seed2<->y2
%
% :Usage:
% ::
%
%     correl_compare_dep_search(seed1,seed2,y1,y2,['alpha',myalpha],['rank'],['mask',maskimage])
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
%
%    % get image names
%    EXPT.subjects = {'ambar_carvalho' 'andreas_nguyen' 'angela_valle' 'brad_wilson' 'dominic_ricci'};
%    EXPT = getfunctnames2(EXPT,'con_0004.img','tmp','SPM_analysis/physical_pain/full_model_gv_p_v_np')
%    self = str2mat(EXPT.tmp{:})
%
%    EXPT = getfunctnames2(EXPT,'con_0003.img','tmp','SPM_analysis/Videos/event_pain_only')
%    other = str2mat(EXPT.tmp{:})
%
%    maskimage = which('scalped_avg152T1_graymatter_smoothed.img')
%
%    % get seeds
%    cd Jamil/GROUP_ANALYSES/Overlaps/Overlap_29_Aug/
%    cl = mask2clusters('con_0002.img');
%    cl = cl(4)
%    cl = extract_contrast_data([{self} {other}],cl);
%    seedself = cl(1).CONTRAST.data(:,1);
%    seedother = cl(1).CONTRAST.data(:,2);
%
%    correl_compare_dep_search(seedself,seedother,self,other,'alpha',.005,'mask',mask);
%
%    % RESULTS:
%    cl = mask2clusters('Correl_seed1_sig.img'); cluster_orthviews(cl,'bivalent');
%
%    % Try the whole thing on ranks:
%    correl_compare_dep_search(seedself,seedother,self,other,'alpha',.005,'mask',mask,'rank');
%
% ..
%    Edited: Tor Wager, Feb 2010; cosmetic fixes only 
% ..

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
                case 'mask', maskimage = varargin{i+1}; varargin{i+1} = [];
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end

    % if image names instead of data matrix, then load data
    if isstr(y1)

        imgs1 = y1;
        imgs2 = y2;
        [y1,volInfo] = iimg_get_data(maskimage,imgs1);
        [y2,volInfo] = iimg_get_data(maskimage,imgs2);

    end




    [N,nvox] = size(y1);


    %[rows,cols,ncorr] = corrcoef_indices(nvox);

    if nargin < 4, dorankdata = 0; end
    if dorankdata
        str = sprintf('Ranking data: Nonparametric correlations'); fprintf(1,str);

        seed1 = rankdata(seed1);
        seed2 = rankdata(seed2);

        for i = 1:nvox
            y1(:,i) = rankdata(y1(:,i));
            y2(:,i) = rankdata(y2(:,i));
        end
    else
        str = sprintf('Assuming continuous data (no ranks).'); fprintf(1,str);
    end

    erase_string(str);

    i = 1; j = 2; k = 3; h = 4;

    str = sprintf('Computing differences among correlations %06d',0); fprintf(1,str);

    diffr = zeros(nvox,1);
    Zstar2 = zeros(nvox,1);
    rr = zeros(nvox,2);
    p1 = zeros(nvox,1);
    p2 = zeros(nvox,1);
    
    % -------------------------------------------------------
    % loop through all voxels
    % -------------------------------------------------------

    for cc = 1:nvox

        fprintf(1,'\b\b\b\b\b\b%06d',cc);

        % get full correlation matrix for the 4 variables involved
        dat = [seed1 y1(:,cc) seed2 y2(:,cc)];

        % get differences between correlation z-values (estimate)
        % diffr = difference between correls
        % rr = correls for condition A, then B (columns)
        % diffrz = Z-score of difference

        [diffr(cc),r,rr(cc,:),diffrz,z,pp] = correlation_diffs(dat,i,j,k,h);

        % p-values for seed analysis of cond1, cond2
        p1(cc) = pp(i,j);
        p2(cc) = pp(k,h);

        s = covcorr(r,i,j,k,h);         % covariance of corr coeffs

        Zstar2(cc) = diffrz * sqrt( (N-3) ./ (2-(2*s)) );

    end

    erase_string(str);

    % all p-values
    % -------------------------------------------------------

    pdiff = 2 * (1 - normcdf(abs(Zstar2)));

    r1 = rr(:,1);
    r2 = rr(:,2);


    % threshold: differences
    [sigdiff,pthrdiff,sigfdrdiff] = apply_threshold(myalpha,pdiff,diffr);

    [sig1,pthr1,sigfdr1] = apply_threshold(myalpha,p1,rr(:,1));
    [sig2,pthr2,sigfdr2] = apply_threshold(myalpha,p2,rr(:,2));

    % write images: differences
    % ----------------------------
    painpag(diffr,volInfo,'outname','Correl_diff_r.img','descrip','Difference in correl with seed A - B');
    iimg_reconstruct_vols(pdiff,volInfo,'outname','Correl_diff_p.img','descrip','Difference in correl with seed A - B');
    iimg_reconstruct_vols(Zstar2,volInfo,'outname','Correl_diff_Z.img','descrip','Z*2 from Steiger, 1980');
    iimg_reconstruct_vols(pdiff,volInfo,'outname','Correl_diff_p.img','descrip','Difference in correl with seed A - B');

    if any(sigdiff)
        iimg_reconstruct_vols(sigdiff,volInfo,'outname','Correl_diff_sig.img','descrip',['Sig uncorrected: ' num2str(myalpha)]);
    else
        disp(['Difference: No significant results at ' num2str(myalpha)]);
    end

    if any(sigfdrdiff)
        iimg_reconstruct_vols(sigfdrdiff,volInfo,'outname','Correl_diff_sigfdr.img','descrip','Sig FDR, p<.05 ');
    else
        disp(['Difference: No significant results at FDR .05']);
    end
    
    % write images: seed 1
    % ----------------------------
    iimg_reconstruct_vols(r1,volInfo,'outname','Correl_seed1_r.img','descrip','seed1 correlation');
    iimg_reconstruct_vols(p1,volInfo,'outname','Correl_seed1_p.img','descrip','Seed1 correl with seed A - B');

    if any(sig1)
        iimg_reconstruct_vols(sig1,volInfo,'outname','Correl_seed1_sig.img','descrip',['Sig uncorrected: ' num2str(myalpha)]);
    else
        disp(['Seed1: No significant results at ' num2str(myalpha)]);
    end

    if any(sigfdr1)
        iimg_reconstruct_vols(sigfdr1,volInfo,'outname','Correl_seed1_sigfdr.img','descrip','Sig FDR, p<.05 ');
    else
        disp(['Seed1: No significant results at FDR .05']);
    end
    
        % write images: seed 2
    % ----------------------------
    iimg_reconstruct_vols(r2,volInfo,'outname','Correl_seed2_r.img','descrip','seed2 correlation');
    iimg_reconstruct_vols(p2,volInfo,'outname','Correl_seed2_p.img','descrip','Seed2 correl with seed A - B');

    if any(sig2)
        iimg_reconstruct_vols(sig2,volInfo,'outname','Correl_seed2_sig.img','descrip',['Sig uncorrected: ' num2str(myalpha)]);
    else
        disp(['Seed2: No significant results at ' num2str(myalpha)]);
    end

    if any(sigfdr2)
        iimg_reconstruct_vols(sigfdr2,volInfo,'outname','Correl_seed2_sigfdr.img','descrip','Sig FDR, p<.05 ');
    else
        disp(['Seed2: No significant results at FDR .05']);
    end
    
    
    
    return






function [sig,pthr,sigfdr] = apply_threshold(myalpha,p,r)

    % uncorrected sign. voxels
    sig = (p <= myalpha - eye(size(p))) .* sign(r);

    % FDR corrected
    pthr = FDR(p,.05);
    if isempty(pthr), pthr = 0; end

    sigfdr = (p <= pthr) .* sign(r);

    return



function [diffr,r,rr,diffrz,z,p] = correlation_diffs(dat,i,j,k,h)
    % get differences between correlation z-values
    [r,p] = corrcoef(dat);

    rr = [r(i,j) r(k,h)];                % correls to be compared
    diffr = -diff(rr);              % negative sign means we get (1) - (2)

    if nargout > 3
        z = .5 .* log( (1+rr) ./ (1-rr) );
        diffrz = -diff(z);               % diff btwn correl z values
    end

    return



function [rows,cols,ncorr] = corrcoef_indices(nvox)
    % upper triangle only
    tmp = triu(ones(nvox));
    tmp = tmp - eye(nvox);
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

function valmat = reconstruct(vals,nvox,ncorr,rows,cols)

    valmat = zeros(nvox);
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
