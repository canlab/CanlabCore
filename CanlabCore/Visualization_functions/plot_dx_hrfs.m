function EXPT = plot_dx_hrfs(EXPT,clusters,varargin)
% :Usage:
% ::
%
%    EXPT = plot_dx_hrfs(EXPT,clusters,[dolegend],[dosave],[dosmooth],[doindiv])
%
% :Inputs:
%
%    uses EXPT.FIR and clusters
%
%    If not found, creates:
%
%    EXPT.FIR.regsofinterest = trial types
%    EXPT.FIR.mcol = colors
%
%    Seems to be a problem with showing the brain slice when it makes the
%    legend as well!  Weird bug.  Optional argument turns legend off.
%
% :Optional Inputs: (all defaults are zero)
%
%   **dolegend:**
%        1 on, 0 off
%
%   **dosave:**
%        1 saves tiff files in timecourse_plots subdir, 0 does not
%
%   **dosmooth:**
%        n smooths hrfs and re-calculates st. errors, 0 plots saved
%        values stored in clusters.HRF.HRF and .STE
%
%   **doindiff:**
%        plot low vs. high groups of individuals (indiv diffs)
%
% :See also: 
%   - extract_dxbeta_data.m, plot_dx_hrfs_indiffs.m
%
% :Examples:
% ::
%
%    % has some smoothing (0 weight @ 3 time pts)
%    plot_dx_hrfs(EXPT,cl(1),0,1,3);

if length(varargin) > 0, dolegend = varargin{1}; else, dolegend = 0;, end
if length(varargin) > 1, dosave = varargin{2}; else, dosave = 0;, end
if length(varargin) > 2, dosmooth = varargin{3}; else, dosmooth = 0;, end
if length(varargin) > 3, doindiv = varargin{4}; else, doindiv = 0;, end

if ~isfield(EXPT.FIR,'mcol')
    disp(EXPT.FIR.dxnames)
    wh = input(['Enter vector of conditions to plot: ']);
    mcol = input('Enter vector of colors, e.g. {''ro-'' ''gd'' etc.}: ');
    EXPT.FIR.regsofinterest = wh;
    EXPT.FIR.mcol = mcol;
else
    wh = EXPT.FIR.regsofinterest;
    mcol = EXPT.FIR.mcol;
    mcol = mcol(wh);
end

if ~isfield(EXPT.FIR,'TR')
    EXPT.FIR.TR = input('Enter TR in s: ');
end

% smoothing, if specified
if dosmooth
    fprintf(1,'Re-averaging cluster HRF and STE: Smoothing = %3.0f\n',dosmooth);
    clusters = dosmoothing(dosmooth,clusters,EXPT);
else
    disp('Using existing HRFs in clusters.HRF.HRF and .STE . No change in baseline subtraction will take effect.');
end

    for j = 1:length(clusters)
        
        
        if ~isempty(wh)
            
            if doindiv
                h = do_indiff_plot(clusters,j,wh,EXPT.FIR.TR,mcol);
            else
                h = do_group_plot(clusters,j,wh,EXPT.FIR.TR,mcol);
            end
            
            if dolegend
                legend(h,EXPT.FIR.dxnames(wh));
            end

            drawnow

            if dosave, savetiff(j), end
            
            
        end
    end
    
    
return
    
    
% ------------------------------------------------    
% MAIN FUNCTION TO MAKE GROUP PLOTS
% ------------------------------------------------
    
    
function h = do_group_plot(clusters,j,wh,TR,mcol)
    

        
    
    
            figure('Color','w'); 
            %h = axes('Position',[0.1100    0.1100    0.3347    0.8150]);
            subplot(1,2,1)
            montage_clusters_maxslice([],clusters(j),{'r'})
            
            %h = axes('Position',[0.6300    0.1100    0.3347    0.8150]); %
            subplot(1,2,2)
            hold on;set(gca,'FontSize',18),xlabel(['Time (s)'])
            
            for k = 1:length(wh)
                %eval(['tmp = clusters(j).HRF.HRF' wh{k} ';'])
                str = (['tmp = clusters(j).HRF.HRF{' num2str(wh(k)) '};']); 
                eval(str)

                x = 0:TR:length(tmp)*TR-1;
                
                h(k) = plot(x,tmp,mcol{k},'LineWidth',2);
                              
                % to plot peak error bar
                %ww = find(tmp==max(tmp));
                %tor_bar_steplot(tmp(ww),mean(clusters(j).HRF.STE{wh(k)}),mcol(wh(k)),x(ww)-1);
                
                % to plot transparent fill
                fill_around_line(tmp,clusters(j).HRF.STE{wh(k)},mcol{k},x);
                
 
            end
            %title(num2str(round(clusters(j).mm_center)))
            
            %subplot(1,3,3)
            %h = bar(clusters(j).HRF.BASE); set(h,'FaceColor',[.8 .8 .8]);
            %tor_bar_steplot(clusters(j).HRF.BASE,clusters(j).HRF.BASESTE,{'k'});
            %title('Baseline estimates')
            %set(gcf,'Position',[109         215        1099         501])
            
 return
 
 
% ------------------------------------------------    
% SAVE A TIFF FILE OF CLUSTER IN APPROPRIATE DIR
% ------------------------------------------------
    
  function savetiff(clnum)
        
        
        if ~(exist('timecourse_plots') == 2)
           mkdir timecourse_plots
        end
        
        clnum = num2str(clnum);
        
        % get unique name for this file
        n = 1;
        name = ['timecourse_cl' clnum '_' num2str(n) '.tiff'];
        while exist([pwd filesep 'timecourse_plots' filesep name]) == 2
            n = n + 1;
            name = ['timecourse_cl' clnum '_' num2str(n) '.tiff'];
        end
            
        saveas(gcf,['timecourse_plots' filesep name],'tiff');
        
    return
        
            
% ------------------------------------------------    
% SMOOTHING AND RE-AVERAGING ON BETAS
% ------------------------------------------------

   function clusters = dosmoothing(smlen,clusters,EXPT)
        % dosmooth(smoothing length)
        % with exponential function, as in smooth_timeseries
        
        st = cumsum([1 EXPT.FIR.numframes]);   
        en = st(2:end) - 1;         % ending values
        st = st(1:end-1);           % starting values

        % set up baseline
        if isfield(EXPT.FIR,'baseline')
            base = EXPT.FIR.baseline;
        else
            base = 0;
        end
        if isempty(base), base = 0;, end
        
        
        for j = 1:length(clusters)
            
            for i = 1:length(st) 
                
                % get individual hrfs for THIS condition
                dat = clusters(j).HRF.indiv(st(i):en(i),:);   
                
                [dat,V] = smooth_timeseries(dat,smlen);
                
                
                % subtract baseline, if specified
                if base
                    bvals = nanmean(dat(base,:));
                    bvals = repmat(bvals,size(dat,1),1);
                    dat = dat - bvals;
                end
                    
                clusters(j).HRF.HRF{i} = nanmean(dat');
                clusters(j).HRF.STE{i} = ste(dat');
                
                
                % hi/low individual diffs, if saved
                if isfield(clusters(j).HRF,'hihrf')
                    % high group
                    dat = clusters(j).HRF.hihrf{i};
                    [dat] = smooth_timeseries(dat',smlen)';
                    clusters(j).HRF.hiHRF{i} = nanmean(dat);
                    clusters(j).HRF.hiSTE{i} = ste(dat);                   
                end
                if isfield(clusters(j).HRF,'lowhrf')
                    % high group
                    dat = clusters(j).HRF.lowhrf{i};
                    [dat] = smooth_timeseries(dat',smlen)';
                    clusters(j).HRF.lowHRF{i} = nanmean(dat);
                    clusters(j).HRF.lowSTE{i} = ste(dat);                   
                end
                
            end
        end
        
    return
     
    
% ------------------------------------------------    
% MAIN FUNCTION TO MAKE INDIV DIFF PLOTS
% ------------------------------------------------
    
    
function h = do_indiff_plot(clusters,j,wh,TR,mcol) 
    
     % LOW
        
        if ~isempty(wh)
            figure('Color','w'); subplot(1,3,1)
            title(num2str(round(clusters(j).mm_center)))
            montage_clusters_maxslice([],clusters(j),{'r'})
            
            ah = subplot(1,3,2);
            hold on;set(gca,'FontSize',18),xlabel(['Time (s)'])
            
            for k = 1:length(wh)
                str = (['tmp = clusters(j).HRF.lowHRF{' num2str(wh(k)) '};']); 
                eval(str)
                
                x = 0:TR:length(tmp)*TR-1;
                
                h(k) = plot(x,tmp,mcol{k},'LineWidth',2);
                              
                % to plot peak error bar
                %ww = find(tmp==max(tmp));
                %tor_bar_steplot(tmp(ww),mean(clusters(j).HRF.STE{wh(k)}),mcol(wh(k)),x(ww)-1);
                
                % to plot transparent fill
                fill_around_line(tmp,clusters(j).HRF.lowSTE{wh(k)},mcol{k},x);
                
 
            end
            title(['Low participants'])
            %legend(h,EXPT.FIR.dxnames(wh));
            
            
            % HIGH
            
            ah(2) = subplot(1,3,3);
            
            hold on;set(gca,'FontSize',18),xlabel(['Time (s)'])
            
            for k = 1:length(wh)
                str = (['tmp = clusters(j).HRF.hiHRF{' num2str(wh(k)) '};']); 
                eval(str)
                
                x = 0:TR:length(tmp)*TR-1;
                
                h(k) = plot(x,tmp,mcol{k},'LineWidth',2);
                              
                % to plot peak error bar
                %ww = find(tmp==max(tmp));
                %tor_bar_steplot(tmp(ww),mean(clusters(j).HRF.STE{wh(k)}),mcol(wh(k)),x(ww)-1);
                
                % to plot transparent fill
                fill_around_line(tmp,clusters(j).HRF.hiSTE{wh(k)},mcol{k},x);
                
 
            end
            title(['High participants'])
            %legend(h,EXPT.FIR.dxnames(wh));
            
            
            equalize_axes(ah);
            set(gcf,'Position',[109         215        1099         501])
            
            drawnow
        end
        
   return
        
    
