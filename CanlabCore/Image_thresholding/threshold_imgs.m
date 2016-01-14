function [P2,P,sigmat,sigmatneg] = threshold_imgs(dd,u,varargin)
% :Usage:
% ::
%
%     function [P2,P,sigmat,sigmatneg] = threshold_imgs(dd,u,[k],['pos' 'neg' 'both'])
%
% :Inputs:
%
%   **dd:**
%        list of filenames (str matrix)
%
%   **u:**
%        height threshold for images
%
%   **k:**
%        extent threshold for contiguous voxels
%
%   **[str]:**
%        'pos' 'neg' or 'both', to values above, below,
%         - and + threshold
%
% Output - threshold t - generic function
%
% also do: plot rob vs ols benefit by tissue class and ols-irls average
%
% :Examples:
% ::
%
%    [P2,P,s,sn] = threshold_imgs(p([5 8],:),tinv(1-.001,10),0,'both');
%    compare_filtered_t([],P2(1,:),P2(2,:))
%    [p2,p1] = threshold_imgs('irls-ols_z_0001.img',norminv(.9),0,'both');
%
%    % compare_filtered_t([],P2(1,:),P2(2,:),p2)
%    h = image_scatterplot(str2mat(P,p1),'avgvs3');
%    xlabel('Average OLS and Robust t-value'), ylabel('Z-score of Robust - OLS difference')
%
%    s = str2mat('rob_tmap_0002.img','rob_tmap_0003.img');
%    P = threshold_imgs(s,tinv(1-.05,36),[],'pos');P = threshold_imgs(s,tinv(1-.05,36),[],'neg');
%
% ..
%    tor wager
% ..

    t = u;      % t-threshold
    iind = 1;   % grandfathered
    jind = 1;   % index of images

    % build output string
    str = ['t_' sprintf('%3.2f',u)]; tmp = find(str=='.'); str(tmp)='-';, tmp = find(str==' ');str(tmp)=[];
    tmp = '0';
    if length(varargin) > 0, k = varargin{1}; if k > 0, tmp = num2str(k);,end, end
    str = [str '_k' tmp];
    if length(varargin) > 1, tmp = varargin{2};, else,tmp = 'pos';, end
    str = [str '_' tmp];

    for j = 1:size(dd,1)

        % go to directory so we can find image, save full path name of image

        cwd = pwd;
        [dr,ff,ee] = fileparts(deblank(dd(j,:)));
        if ~isempty(dr), cd(dr), end
        dr = pwd;
        name = which([ff ee]);
        outname = fullfile(dr,[ff '_filt_' str ee]);
        if j ==1, P = name; P2 = outname; else, P = str2mat(P,name);,P2 = str2mat(P2,outname);,end

        cd(cwd)

        % read the image

        V = spm_vol(name);
        v = spm_read_vols(V);

        coln{jind} = V.fname;

        %t = tinv(1-u,df);
        disp(['Threshold is ' num2str(t)])

        if length(varargin) > 1
            if strcmp(varargin{2},'pos')
                v(v < t) = NaN;
            elseif strcmp(varargin{2},'neg')
                v(-v < t) = NaN;
            elseif strcmp(varargin{2},'both')
                v(abs(v) < t) = NaN;
            else error('Unknown string input.')
            end
        else
            % save both pos and neg
            v(abs(v) < t) = NaN;
            %figure;hist(v(v ~= 0))
        end

        sigmat(iind,jind) = sum(~isnan(v(:)) & v(:) > 0);
        sigmatneg(iind,jind) = sum(~isnan(v(:)) & v(:) < 0);
        fprintf(1,'Sig + voxels: %3.0f, - vox: %3.0f', sigmat(iind,jind),sigmatneg(iind,jind))

        V.descrip = ['Image thresholded at ' num2str(u)];

        if length(varargin) > 0
            k = varargin{1};
            if k > 0
                fprintf(1,' Thr. k = %3.0f ',k)

                V.descrip = ['Image thresholded at ' num2str(u) ', k = ' num2str(k)];

                if any(v(:) > 0)
                    mask = clusterSizeMask(k,v);
                else
                    mask = zeros(size(v));
                end

                if any(v(:) < 0)
                    mask2 = clusterSizeMask(k,-v);	% to save negative results, filtered out earlier if not requested
                else
                    mask2 = zeros(size(v));
                end

                mask = mask | mask2;
                v = mask .* v;

                sigmat(iind,jind) = sum(~isnan(v(:)) & v(:) > 0);
                sigmatneg(iind,jind) = sum(~isnan(v(:)) & v(:) < 0);
                fprintf(1,'\tAfter size thresh: + %3.0f, - %3.0f',sigmat(iind,jind),sigmatneg(iind,jind))
            end

            V.descrip = [V.descrip ' and k = ' num2str(varargin{1})];
        end
        fprintf(1,'\n')

        jind = jind + 1;

        if length(varargin) > 1, V.descrip = [V.descrip varargin{2}];, end

        V.fname = outname;
        spm_write_vol(V,v);
        disp(['Written: ' V.fname])


    end        % img

    return

