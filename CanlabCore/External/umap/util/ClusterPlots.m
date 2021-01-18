%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

classdef ClusterPlots < Plots
    
    properties(SetAccess=private)
        clusterIds;
        numClues;
        numClusters;
        l;
        data;
        clr;
        clues;
        clue;
        ranks;
        clrMap;
        labelMap;
        step;
    end
    
    methods
        function this=ClusterPlots(data, clusterIds, clrMap, labelMap)
            this.clrMap=clrMap;
            this.data=data;
            clues=unique(clusterIds);
            this.clusterIds=clusterIds;
            this.clues=clues;
            this.numClusters=sum(clues>=0);
            numClues=length(clues);
            this.numClues=numClues;
            this.N=numClues;
            this.Hs=zeros(1, numClues);
            this.numClues=numClues;
            cnts=zeros(1,numClues);
            for i=1:numClues
                cnts(i)=sum(clusterIds==clues(i));
            end
            this.cnts=cnts;
            this.setMinMax;
            [~,this.ranks]=sort(cnts);
            nClrs=size(clrMap,1);
            this.step=floor(nClrs/numClues);
            this.labelMap=labelMap;
        end
        
        function H=plot3D(this, ax, i)
            H=plot3(ax, this.data(this.l,1), this.data(this.l,2), this.data(this.l,3), '.', ...
                'markerSize', 5, 'lineStyle', 'none', ...
                'markerEdgeColor', this.clr, ...
                'markerFaceColor', this.clr);
            this.Hs(i)=H;
            this.ax=ax;
        end
        
        function init2D(this)
            this.otherHs=zeros(1, this.numClues);
        end
        
        function H=plot2D(this, ax, i)
            H=plot(ax, this.data(this.l, 1), ...
                this.data(this.l, 2), '.', ...
                'visible', 'off',...
                'markerSize', 5, ...
                'lineStyle', 'none', ...
                'markerEdgeColor', this.clr, ...
                'markerFaceColor', this.clr);
            this.otherHs(i)=H;
        end
        
        function plot2DHere(this, ax, i)
            this.Hs(i)=plot(ax, this.data(this.l, 1), ...
                this.data(this.l, 2), '.', ...
                'markerSize', 5, ...
                'lineStyle', 'none', ...
                'markerEdgeColor', this.clr, ...
                'markerFaceColor', this.clr);
        end
        
        
        function setCluster(this, i)
            this.clue=this.clues(i);
            if this.clue<1
                this.clr=[.2 .2 .2];
            elseif ~isempty(this.labelMap)
                clr_=this.labelMap.get(...
                    [num2str(this.clues(i)) '.color']);
                if isempty(clr_)
                    this.clr=[.95 .9 .99];
                else
                    this.clr=str2num(clr_)/256;
                end
            else
                ranking=find(i==this.ranks,1);
                this.clr=this.clrMap(ranking * this.step, :);
            end
            if any(this.clr>1)
                this.clr(this.clr>1)=1;
            end
            this.l=this.clusterIds==this.clue;
        end
        
        function refresh(this, ax2D)
            for i=1:this.numClues
                this.setCluster(i);
                this.Hs(i)=this.plot2D(ax2D);
            end
        end
        
        function names=getNames(this, extractTraining)
            if isempty(this.names)
                if nargin<2
                    extractTraining=false;
                    word='background';
                else
                    word='untrained test subset';
                end
                lblMap=this.labelMap;
                doLabels=~isempty(lblMap);
                [addTrainingHtml, sup1, sup2, trStart, trEnd]=...
                    LabelBasics.AddTrainingHtml(lblMap, extractTraining);
                this.names=cell(1, this.numClues);
                for i=1:this.numClues
                    if this.clues(i)==0
                        this.names{i}=word;
                    elseif this.clues(i)<0
                        this.names{i}=[word ' ID=' num2str(0-this.clues(i))];
                    elseif doLabels
                        [key, ~, keyTraining]=LabelBasics.Keys(this.clues(i));
                        name=lblMap.get(key);
                        if isempty(name)
                            name=['Cluster ' num2str(this.clue)];
                        end
                        if addTrainingHtml
                            name=strrep(name, '^{', sup1);
                            name=strrep(name, '}', sup2);
                            this.names{i}=[name trStart lblMap.get(...
                                keyTraining) trEnd];
                        else
                            this.names{i}=name;
                        end
                    else
                        this.names{i}=['Cluster ' num2str(this.clues(i))];
                    end
                end
            end
            names=this.names;
        end
        
    end
    
    methods(Static)
        
        function plots=Go(ax3D, data, clusterIds, labelMap, ...
                xLabel, yLabel, zLabel, doLegend, ax2D, gray, doDensity,...
                tickOff, doJavaLegend)
            if nargin<13
                doJavaLegend=true;
                if nargin<12
                    tickOff=true;
                    if nargin<11
                        doDensity=true;
                        if nargin<10
                            gray=false;
                            if nargin<9
                                ax2D=[];
                                if nargin<8
                                    doLegend=true;
                                end
                            end
                        end
                    end
                end
            end
            if ~gray
                clrMap=jet(256);
            else
                clrMap=bone(256);
            end
            clrMap=clrMap(1:240, :);
            if isempty(clusterIds) %everything is background?
                clusterIds=zeros(size(data,1), 1);
            end
            plots=ClusterPlots(data, clusterIds, clrMap, labelMap);
            if isempty(ax3D)
                fig2=Gui.NewFigure(true);
                set(fig2, 'name', ...
                    [num2str(plots.numClusters) ' clusters found`...']);
                op=get(fig2, 'OuterPosition');
                w=op(3);
                h=op(4);
                set(fig2, 'OuterPosition', [op(1)+.1*w op(2)-.1*h, ...
                    w*.8, h*.8]);
                ax3D=gca;
            end
            cla(ax3D, 'reset');
            hold(ax3D, 'on');
            if tickOff
                set(ax3D, 'xtick', [], 'ytick', [], 'zTick', [])
            end
            mns=min(data);
            mxs=max(data);
            xlim(ax3D, [mns(1) mxs(1)]);
            ylim(ax3D, [mns(2) mxs(2)])
            n_components=size(data,2);
            if n_components>2
                zlim(ax3D, [mns(3) mxs(3)])
                view(ax3D, [1 1 1]);
            end
            if ~isempty(ax2D)
                plots.init2D;
            end
            if n_components>2 && (~doDensity || doLegend)
                names=plots.getNames(doLegend && doJavaLegend);
                if isempty(ax2D)
                    for i=1:plots.numClues
                        try
                        plots.setCluster(i);
                        plots.plot3D(ax3D, i);
                        catch ex
                            disp(ex);
                        end
                    end
                else
                    for i=1:plots.numClues
                        plots.setCluster(i);
                        plots.plot3D(ax3D, i);
                        plots.plot2D(ax2D, i);
                    end
                end
                if doLegend
                    plots.javaLegend=Plots.Legend(plots, names, [], -0.01, 0.061, ...
                        true, sum(plots.cnts), [], [], doJavaLegend);
                end
            elseif n_components==2
                if doLegend || ~doDensity
                    names=plots.getNames;
                    for i=1:plots.numClues
                        plots.setCluster(i);
                        plots.plot2DHere(ax3D, i);
                    end
                    plots.ax=ax3D;
                    if doLegend
                        Plots.Legend(plots, names, [], -0.01, 0.061, ...
                            true, sum(plots.cnts), [],[], doJavaLegend);
                    end
                else
                    for i=1:plots.numClues
                        plots.setCluster(i);
                        plot(ax3D, data(plots.l,1), data(plots.l,2), '.', ...
                            'markerSize', 1, 'lineStyle', 'none', ...
                            'markerEdgeColor', plots.clr, ...
                            'markerFaceColor', plots.clr);
                    end
                end
            else
                [D,~,I]=Density.Get3D(data);
                maxD=max(D(D(:)>0));
                minD=min(D(D(:)>0));
                rangeD=maxD-minD;
                for i=1:plots.numClues
                    plots.setCluster(i);
                    if n_components>2
                        if plots.clue<1 
                            plots.plot3D(ax3D, i);
                            continue;
                        end
                        l2=false(1, length(plots.l));
                        di=I(plots.l,:);
                        d=zeros(1, plots.cnts(i));
                        for j=1:plots.cnts(i)
                            d(j)=D(di(j,1), di(j,2), di(j,3));
                        end
                        ratios=(d-minD)./rangeD;
                        denominator=10;
                        for j=1:denominator
                            ratio=j/denominator;
                            ll=ratios<ratio & ratios>=(j-1)/denominator;
                            l2(plots.l)=ll;
                            clrRatio=1-(j-1)/(denominator+2);
                            clr2=clrRatio*plots.clr;
                            
                            plot3(ax3D, data(l2,1), data(l2,2), data(l2,3), '.', ...
                                'markerSize', 2, 'lineStyle', 'none', ...
                                'markerEdgeColor', clr2, ...
                                'markerFaceColor', clr2);
                            fprintf('%d of %d', sum(l2), plots.cnts(i));
                            fprintf('\n');
                        end
                        plots.ax=ax3D;
                    else
                        
                    end
                end
            end
            if nargin>2
                xlabel(ax3D, xLabel);
                if nargin>3
                    ylabel(ax3D, yLabel)
                    if nargin>4
                        if n_components>2
                            zlabel(ax3D, zLabel);
                        end
                    end
                end
            end
            grid(ax3D, 'on')
            set(ax3D, 'plotboxaspectratio', [1 1 1])            
        end
    end
end