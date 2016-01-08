function EXPT = scnlab_norm_check3(wt1,subjlabels,template,mask,varargin)
% WARNING:  scnlab_norm_check3 is deprecated! All improvements are being placed in scnlab_norm_check.
%
% :Usage:
% ::
%
%     EXPT = scnlab_norm_check3(wt1,subjlabels,template,mask,[print out MI])
%
% :Inputs:
%
%   **wt1:**
%        char array of wT1.img files, one per line
%
%   **subjlabels:**
%        cell array of subject labels
%
%   **template:**
%        template img that everything has been normalized to (usually the avg152T1.img file)
%
%   **mask:**
%        image to mask with
%
% :Optional Input:
%
%   print out mutual information table - flag for whether or not to print out the MI table; defaults to 0
%
% :Examples:
% ::
%
%    cd(studyroot); % wherever your study root is
%    wt1s = filenames('hr*/structural/wT1.img', 'char', 'absolute'); % assuming that hr is your study code
%    subjlabels = filenames('hr*');
%    template = which('avg152T1.nii');
%    mask = filenames('scalped_avg152T1_graymatter.img', 'char', 'absolute'); % set to wherever your mask is...
%
%    EXPT = scnlab_norm_check3(wt1, subjlabels, template, mask);
%
% THIS FUNCTION IS DEPRECATED; SCNLAB_NORM_CHECK IS PREFERRED

    warning('scnlab_norm_check3 is deprecated! All improvements are being placed in scnlab_norm_check.');
    return
        
    doprintMI = 0;

    if isempty(wt1)
        wt1 = spm_get(Inf,'w*img','Cannot find subjects'' normalized T1 imgs: please select');
    end
    if isempty(template)
        template = spm_get(1,'*img','Cannot find template for comparison: please select');
    end
    if isempty(mask)
        mask = spm_get(1,'*img','Cannot find mask: please select');
    end

    if length(varargin) > 0
        doprintMI = varargin{1};
    end

    
    
    [ds,g,mystd,d,d2,c,c2,mi,b,eigv,eigval] = compare_subjects(wt1,mask,0,'Norm T1 vs. template',1,subjlabels,template);

    if(doprintMI)
        nms = {'Avg. abs dev from mean' 'Correl with template' 'MI with template'};
        print_matrix([d' c' mi'], nms, subjlabels)
    end

    tmp = sortrows(ds,2);
    tmp = tmp(:,1) - tmp(:,3);

    EXPT.template = template;
    EXPT.SUBJECT.wt1 = wt1;
    EXPT.SUBJECT.meanfunct = [];
    EXPT.SUBJECT.norm_vs_template = [tmp mi' c'];
    EXPT.SUBJECT.global_t1 = g';
    EXPT.SUBJECT.std_t1 = mystd';
    EXPT.SUBJECT.names = ([{'Dist. from group, actual-expected chi2'}, {'Mutual info with template'},{'Correlation with template'}]);

    local_save_plots('Norm_vs_Temp');

    local_plot_results(EXPT, template, wt1, subjlabels)

end


function EXPT = local_get_images(EXPT, templ)


    % normalized structural vs. canonical
    % normalized structural vs. mean norm functional
    % mean norm functional vs. group

    % default
    P = 'avg152T1.img';  % canonical
    normhi = 'wT1.img';                      % normalized hi-res T1 - check against canonical
    %hires = 'het1spgr.img';                 % hi-res T1: check against functional
    funct = 'mean_swravol.img';             % mean functional in first scan directory
    wcard = 'swra*img';                     % wildcard to search for funct images if mean doesn't exist

    EXPT = get_expt_subdir(EXPT);
    EXPT.SUBJECT.wt1 = [];
    EXPT.SUBJECT.meanfunct = [];

    % canonical
    EXPT.template = which(P);
    if isempty(EXPT.template), disp(['Warning: cannot find canonical image: ' P]),end



    % get normalized images
    disp('Normalized Hi-res images:')

    for i = 1:length(EXPT.subjects)

        str = fullfile(pwd,EXPT.subjects{i},'anatomy',normhi);
        d = dir(str);

        if ~isempty(d)

            if isempty(EXPT.SUBJECT.wt1)
                EXPT.SUBJECT.wt1 = str;
            else
                EXPT.SUBJECT.wt1 = str2mat(EXPT.SUBJECT.wt1,str);
            end

        else
            disp(['Missing: ' str])
        end

    end


    % get functional images

    disp('Mean functional images:')

    for i = 1:length(EXPT.subjects)

        str = fullfile(pwd,EXPT.subjects{i},funct);
        d = dir(str);

        if isempty(d)
            disp(['Missing and attempting to create: ' str])

            str2 = fullfile(pwd,EXPT.subjects{i},'scan1',wcard);
            df = dir(str2);

            if isempty(df),
                disp(['Warning - No images: ' str2]);
            else
                str2 = repmat(fullfile(pwd,EXPT.subjects{i},['scan1' filesep]),length(df),1);
                str2 = [str2 str2mat(df.name)];

                tor_spm_mean_ui(str2,str);
            end
        else

            if isempty(EXPT.SUBJECT.meanfunct)
                EXPT.SUBJECT.meanfunct = str;
            else
                EXPT.SUBJECT.meanfunct = str2mat(EXPT.SUBJECT.meanfunct,str);
            end

        end

    end


end




function local_plot_results(EXPT, template, wt1, subjlabels)


    % plot results
    % ____________________________________________________________

    tmp = scale(EXPT.SUBJECT.norm_vs_template); tmp2 = tmp;

    % save best two and worst two overall
    overall = tmp;
    overall(:,2:3) = -1*overall(:,2:3);
    overall = mean(overall,2);   % higher is worse
    EXPT.SUBJECT.overall_t1 = overall;
    overall2 = overall;
    tmp = find(overall == min(overall)); wh(1) = tmp; overall(tmp) = NaN;
    tmp = find(overall == min(overall)); wh(2) = tmp; overall(tmp) = NaN;
    tmp = find(overall == max(overall)); wh(3) = tmp; overall(tmp) = NaN;
    tmp = find(overall == max(overall)); wh(4) = tmp; overall(tmp) = NaN;
    P = str2mat(template,wt1(wh,:));
    EXPT.SUBJECT.best_t1 = P(1:2,:);
    EXPT.SUBJECT.worst_t1 = P(3:4,:);

    spm_check_registration(P);

    % relies on the figure created by compare_subjects... too closely coupled
    h = sort(get(gcf,'Children'));
    axes(h(3)); title('Canonical template');
    axes(h(7)); title(['Best: ' subjlabels{wh(1)}]);
    axes(h(11)); title(['Best: ' subjlabels{wh(2)}]);
    axes(h(15)); title(['Worst: ' subjlabels{wh(3)}]);
    axes(h(19)); title(['Worst: ' subjlabels{wh(4)}]);

    h2 = axes('Position',[.6 .05 .3 .25]);
    plot(tmp2,1:length(tmp2)); %try legend(EXPT.subjects); catch end
    hold on; plot(overall2,1:length(overall),'ko-','LineWidth',2);drawnow
    set(gca,'YTickLabel',subjlabels,'YTick',1:length(overall))
    xlabel('Value (higher is worse for black and blue lines)')

    saveas(gcf,'Best_and_worst_t1','fig');
    saveas(gcf,'Best_and_worst_t1','png');

end


function local_save_plots(plot_name)
    try
        saveas(gcf,sprintf('subjplot_%s', plot_name),'fig');
        saveas(gcf,sprintf('subjplot_%s', plot_name),'png');
        saveas(gcf,sprintf('multdist_%s', plot_name),'fig');
        saveas(gcf,sprintf('multdist_%s', plot_name),'png');
    catch
        disp('Cannot save images');
    end
end
