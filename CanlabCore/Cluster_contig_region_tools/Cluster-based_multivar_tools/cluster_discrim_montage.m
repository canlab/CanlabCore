function cluster_discrim_montage(cl,thr,varargin)
% cluster_discrim_montage(cl,thr,[scores],[behavior])
%
% For viewing clusters with nonzero entries in each column of thr
% with color assignments based on values in thr
% i.e.,
% thr is thresholded component loadings, zero if not significant
%
% optional: scores and behavior arguments
% If entered, partial correlation plot will be shown
% using discrim_plot.m
% 
% see also:
% cluster_discrim
% discrim_plot
%
% examples:
% to run stand-alone from output of cluster_discrim:
% cluster_discrim_montage(cl,D.PCR.threshcomps,D.PCR.sscore,D.beh);

    % ---------------------------------------------     
    % save clusters and get figures for each
    % ---------------------------------------------
    
    
    %mycols = {'r' 'b' 'g' 'y' 'c' 'm' 'k' 'k' 'k' 'k'}; str = [];
    for i = 1:size(thr,2)
        
        fprintf(1,'\n--------------------------------\n');
        fprintf(1,'Sig. Component %3.0f',i);
        fprintf(1,'\n--------------------------------\n');
        
        wh = find(thr(:,i));    % positive and negative loading clusters
        
        fprintf(1,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n', ...
            'Cluster','Name','x','y','z','V','Comp.r','Comp.p', 'Beh.r', 'Beh.p','BehRobr', 'BehRobp','Zavg','Inc/Dec');
                
        for j = 1:length(cl)
            % save loading in Z-values
            cl(j).Z = thr(j,i) .* ones(1,size(cl(j).XYZ,2));
            
            if thr(j,i)
                % print output
                fprintf(1,'Cluster %3.0f\t',j);
                if isfield(cl,'shorttitle'), fprintf(1,'%s\t',cl(j).shorttitle),end
                if isfield(cl,'mm_center'),fprintf(1,'%3.0f\t%3.0f\t%3.0f\t',cl(j).mm_center);,end
                if isfield(cl,'numVox'),fprintf(1,'%3.0f\t',cl(j).numVox);,end            
                fprintf(1,'%3.2f\t',thr(j,i))
                
                if length(varargin) > 0
                    scores = varargin{1};
                    % correlations with components; we already have correls
                    % in thr input; get p-values
                    [dummy,dummy,r,p,rrob,prob]  = partialcor([scores],cl(j).timeseries,1);
                    fprintf(1,'%3.4f\t',[p])
                end
                
                if length(varargin) > 0
                    beh = varargin{2};
                    % correlations with behavior
                    [dummy,dummy,r,p,rrob,prob]  = partialcor(cl(j).timeseries,beh,1);
                    fprintf(1,'%3.2f\t%3.4f\t%3.2f\t%3.4f\t',[r p rrob prob])
                end
                
                % t-test
                [hh,p,ci,stats] = ttest(cl(j).timeseries);
                strs = {'Decrease' '--' 'Increase'};
                if ~hh, strs = strs{2};,elseif stats.tstat >0,strs=strs{3};,else,strs = strs{1};end
                fprintf(1,'%3.2f\t%s\t',stats.tstat,strs)
                
                fprintf(1,'\n');
            end
        end
        
        % orthviews
        whpos = find(thr(:,i) > 0);
        whneg = find(thr(:,i) < 0);
        try
            
        if any(thr(wh,i) > 0),
            cluster_orthviews(cl(whpos),{[1 0 0]});
            if any(thr(wh,i) < 0),
                cluster_orthviews(cl(whneg),{[1 0 0]},'add');   % should have Z-scores reversed, so be neg
            end
        elseif any(thr(wh,i) < 0),
            cluster_orthviews(cl(whneg),{[1 0 0]});
        end
        
        catch
            disp('cannot make cluster_orthviews'),
        end
        
        
        % scatterplot, if we have necessaries
        
        if length(varargin) > 0
            try, scores = varargin{1}; beh = varargin{2};, catch, error('Enter scores and beh as 3rd and 4th args or stick w/two args.'),end
            if ~exist('h'),h=[];,end
            if ishandle(h)
                axes(h); cla;
            else
                h = axes('Position',[.52 .08 .45 .38]);
            end
            y2 = discrim_plot(scores,beh,0,i);    % partial residual plot of ith column of scores vs. beh
            xlabel(['Component ' num2str(i)])
        end
            
        % montage
        
        montage_clusters([],cl(wh),{'r'},[2 2]);    % should plot z-scores in red/blue
        drawnow
        
        ans = input('Press return to continue.');
        
    end