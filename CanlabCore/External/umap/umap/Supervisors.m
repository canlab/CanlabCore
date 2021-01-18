%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef Supervisors < handle
    properties(Constant)
        DEBUG_QF=false;
        DEV_UNIT_LIMIT=4;
        BORDER_DEV_UNIT_LIMIT=.2;
        VERBOSE=false;
        COLOR_NEW_SUBSET=[220 220 223];
        CLUSTER_LABELS_IF_EQUAL=false;
        PROP_STOP_WARNING='ustStopNoNameWarn';
        MIN_FREQUENCY=10;
        CLUSTER_UNSUPERVISED=true;
    end
    
    properties
        density;
        embedding;
        inputData;
        description;
        context;
        knnIndices;
        supervisingMean;
        supervisedMean
        verbose='graphic';
    end
    
    properties(SetAccess=private)
        uuid;
        sourceUuid; % universal id of source of supervising labels.  This is a
                    % immutable universally unique identifier
        sourceSubId;
        sourceDataId;
        sourceDescription;
        ids;
        labelMap;
        mdns;
        mads;
        means;
        stds;
        cnts;
        N;
        xLimit;
        yLimit;
        zLimit;
        labels;
        nnLabels;% nearest neighbor match
        nnUnreducedLabels;
        matchType=1;
        warnNoName;
        warnings=0;
        fromD;
        contourPercent=10;
        clusterMethod2D;
        clusterDetail;
        epsilon;
        minopts;
        dbscanDistance;
        trainingSetPlots;
        unsupervisedClues;
    end
    
    methods
        function overrideColors(this, colorFile, beQuiet)
            if nargin<3
                beQuiet=false;
            end
            ColorsByName.Override(this.labelMap, colorFile, beQuiet);
        end
        
        function printSupervisors(this)
            nSupervisors=length(this.ids);
            for i=1:nSupervisors
                name=this.labelMap.get(java.lang.String(num2str(this.ids(i))));
                fprintf('#%d.  %d=%s %d events\n', i, this.ids(i), ...
                    name, this.cnts(i));
            end
        end
        
        function N=size(this)
            N=length(this.ids);
        end
        function this=Supervisors(labels, labelMap, embedding, ax, fromD)
            this.fromD=fromD;
            D=size(embedding,2);
            this.labelMap=labelMap;
            this.ids=unique(labels);
            this.labels=labels;
            nUniqueLabels=length(this.ids);
            this.cnts=zeros(nUniqueLabels,1);
            this.mdns=zeros(nUniqueLabels,D);
            this.means=zeros(nUniqueLabels,D);
            this.mads=zeros(nUniqueLabels,D);
            this.stds=zeros(nUniqueLabels,D);
            BasicMap.Global.set(Supervisors.PROP_STOP_WARNING, 'false');
            for i=1:nUniqueLabels
                id=this.ids(i);
                l=labels==id;
                this.cnts(i)=sum(l);
                this.mdns(i,:)=median(embedding(l,:));
                this.mads(i,:)=mad(embedding(l,:),1);
                this.means(i,:)=mean(embedding(l,:));
                this.stds(i,:)=std(embedding(l,:),1);
            end
            this.N=nUniqueLabels;
            if nargin<4 || isempty(ax) || ~ishandle(ax)
                mx=max(embedding);
                mn=min(embedding);
                this.xLimit=[mn(1) mx(1)];
                this.yLimit=[mn(2) mx(2)];
                if D>2
                    this.zLimit=[mn(3) mx(3)];
                end
            else
                this.xLimit=xlim(ax);
                this.yLimit=ylim(ax);
                if D>2
                    this.zLimit=zlim(ax);
                end
            end
            this.embedding=embedding;
        end
        
        function setUniversalIdentifier(this, uuid, subId, dataId, description)
            assert(isempty(this.sourceUuid), ...
                'universal id is not empty (immutable means set once)');
            this.sourceUuid=uuid;
            if nargin>1 
                this.sourceSubId=subId;
                if nargin>2
                    this.sourceDataId=dataId;
                    if nargin>3
                        this.sourceDescription=description;
                    end
                end
            end
            this.uuid=java.util.UUID.randomUUID;
        end
        

        function [name, color, label]=getNameByMedian(this, mdn)
            [~, iMeHearty]=pdist2(this.density.clusterMdns, mdn, ...
                'euclidean', 'smallest', 1);
            name=this.density.clusterNames{iMeHearty};
            color=this.density.clusterColors{iMeHearty};
            label=this.density.clusterLabels{iMeHearty};
        end        
        
    end
    
    methods(Static)        
        function clr=NewColor(id)
            clr=Supervisors.COLOR_NEW_SUBSET+id;            
            clr(clr<0)=0;
            clr(clr>252)=252;
        end
        
        function [mins, maxs]=GetMinsMaxs(data)
            [mins, maxs]=MatBasics.GetMinsMaxs(data, .15);
        end
        
        function plots=Plot(data, lbls, lblMap, fromD, ...
                umap, ax, doJavaLegend, tickOff, contourPercent)
            if nargin<9
                contourPercent=10;
                if nargin<8
                    tickOff=false;
                    if nargin<7
                        doJavaLegend=[];
                        if nargin<6
                            ax=[];
                            if nargin<5
                                umap=[];
                                if nargin<4
                                    fromD=[];
                                end
                            end
                        end
                    end
                end
            end
            if isempty(ax)
                fig=Gui.NewFigure(true);
                ax=gca;
                op=get(fig, 'OuterPosition');
                w=op(3);
                h=op(4);
                set(fig, 'OuterPosition', [op(1)+.1*w op(2)-.1*h, ...
                    w*.8, h*.8]);
                Gui.FitFigToScreen(fig);
            end
            toD=size(data, 2);
            
            dimInfo=sprintf('  %dD\\rightarrow%dD', fromD, toD);
            
            xLabel=['UMAP-X' dimInfo];
            yLabel=['UMAP-Y' dimInfo];
            if toD>2
                zLabel=['UMAP-Z' dimInfo];
                plots=ClusterPlots.Go(ax, data, lbls, lblMap, ...
                    xLabel, yLabel, zLabel, ~isempty(doJavaLegend), [], ...
                    false, false, tickOff, doJavaLegend);
            else
                [~,javaLegend,~,~,plots]=ProbabilityDensity2.DrawLabeled(ax,...
                    data, lbls, lblMap, contourPercent, true, [], [], ...
                    -0.01, 0.061, doJavaLegend);
                plots.javaLegend=javaLegend;
                if ~tickOff
                    grid(ax, 'on');
                end
                xlabel(ax, xLabel)
                ylabel(ax, yLabel);
            end
            if ~isempty(umap)
                umap.adjustLims(ax, data);
            end
            set(ax, 'plotboxaspectratio', [1 1 1])
        end

        
    end
    
    methods
        function [names, lbls, clrs]=getQfTraining(this)
            names={};
            clrs=[];
            lbls=this.labels;
            N_=length(this.ids);
            for i=1:N_
                id=this.ids(i); 
                if id>0      
                    if this.cnts(i)>=Supervisors.MIN_FREQUENCY
                        key=num2str(id);
                        names{end+1}=this.labelMap.get(java.lang.String(key));
                        if isempty(names{end})
                            names{end}=key;
                        end
                        clr_=this.labelMap.get([key '.color']);
                        if isempty(clr_)
                            clrs(end+1,:)=[.95 .9 .99];
                        else
                            clrs(end+1,:)=str2num(clr_)/256;
                        end
                    else
                        lbls(this.labels==id)=0;
                    end
                end
            end
        end
        
        function statement=computeNearestNeighborsUnreduced(...
                this, testData, pu)
            [R,D]=size(testData);
            if length(this.nnUnreducedLabels)~=R
                this.nnUnreducedLabels=[];
            end
            if ~isempty(this.nnUnreducedLabels)
                return;
            end 
            txt=['Finding nearest neighbors in ' num2str(D) 'D space'];
            if nargin<3
                pu=PopUp(txt, 'center', 'Supervising...', false);
            elseif ~isempty(pu)
                old=pu.label.getText;
                pu.label.setText(txt);
            end            
            if size(this.knnIndices, 1) ~= size(testData, 1) ...
                || mean(this.inputData(:)) ~= this.supervisingMean ...
                || mean(testData(:)) ~= this.supervisedMean ...
                [~,II]=pdist2(this.inputData, testData, 'euclidean', 'Smallest', 1);
            else
                II=this.knnIndices(:,1);
            end
            nn=this.labels(II);
            this.nnUnreducedLabels=nn;
            this.matchType=4;
            if nargin<3
                pu.close;
            elseif ~isempty(pu)
                pu.label.setText(old);
            end
                
            if nargout==0
                return;
            end
            statement=[];
            
            try
                if length(this.nnLabels)==R
                    lbls=this.nnLabels;
                elseif length(this.density.labels)==R
                    lbls=this.density.labels;
                else
                    % all done
                    return;
                end
                disagreements=sum(nn~=lbls);
                disagreementPercent=String.encodePercent(...
                    disagreements, R,0);
                statement=sprintf(['<html>%s of the %s events have different'...
                    '<br>supervisor labels when comparing raw data<br>'...
                    'that when comparing embedded data with pdist2!'...
                    '<hr></html>'], disagreementPercent, ...
                    String.encodeInteger(R));
                msg(statement, 0, 'south west', 'Matching raw data');
            catch ex
                disp(ex);
            end
        end

        function [names, lbls, clrs]=getQfTrained(this, data)
            names={};
            clrs=[];
            [lbls, lblMap]=this.supervise(data);
            lbls(lbls<0)=0;
            ids_=unique(lbls);
            cnts_ = UmapUtil.discreteCount(lbls, ids_);
            N_=length(ids_);
            for i=1:N_
                id=ids_(i);
                if id>0
                    key=num2str(id);
                    name=lblMap.get(java.lang.String(key));
                    if isempty(name)
                        name=['Subset #' key];
                    end
                    if cnts_(i)>=Supervisors.MIN_FREQUENCY
                        key=num2str(id);
                        names{end+1}=name;
                        try
                            clrs(end+1,:)=str2num(...
                                lblMap.get([key '.color']))/256;
                        catch 
                            clrs(end+1,:)=[.1 .1 .1];
                        end
                    else
                        if ~strcmpi(this.verbose, 'none')
                            warning('Ignoring "%s" since it has %d events...', ...
                                name, cnts_(i));
                        end
                        lbls(lbls==id)=0;
                    end
                else
                    lbls(lbls==id)=0;
                end
            end
        end
        
        
        function [names, lbls, clrs]=getOtherTrained(this, lbls)
            names={};
            clrs=[];
            lblMap=this.labelMap;
            lbls(lbls<0)=0;
            ids_=unique(lbls);
            cnts_ = UmapUtil.discreteCount(lbls, ids_);
            N_=length(ids_);
            for i=1:N_
                id=ids_(i);
                if id>0
                    key=num2str(id);
                    name=lblMap.get(java.lang.String(key));
                    if isempty(name)
                        name=['Subset #' key];
                    end
                    if cnts_(i)>=Supervisors.MIN_FREQUENCY
                        names{end+1}=name;
                        value=lblMap.get([key '.color']);
                        try
                            clrs(end+1,:)=str2num(value)/256;
                        catch
                            clrs(end+1,:)=[.1 .1 .1];
                        end
                    else
                        lbls(lbls==id)=0;
                        if ~strcmpi(this.verbose, 'none')
                            warning('Ignoring "%s" since it has %d events...', ...
                                name, cnts_(i));
                        end
                    end
                else
                    lbls(lbls==id)=0;
                end
            end
        end
    end
    
    methods(Access=private)
        
        function qf=qfMatchWithClusters(this, data, dns, numClusters, ...
                clusterIds, pu)
            D=size(data,2);
            if nargin<5
                pu=[];
                if nargin<3
                    [numClusters,clusterIds,dns]=this.findClusters(data, pu);
                end
            end
            qf=[];
            if isempty(this.embedding)
                return;
            end
            this.density=dns;
            [R,C]=size(clusterIds);
            if C>1 && R==1
                clusterIds=clusterIds';
            end
            [tNames, lbls, clrs]=this.getQfTraining;
            if Supervisors.CLUSTER_UNSUPERVISED
                [lbls, tNames, clrs, unsprvStartLabel]=...
                    this.addUnsupervisedClusters(lbls, tNames, clrs, pu);
            end
            if isequal(data, this.embedding(:,1:D))
                matchStrategy=2;
            else
                matchStrategy=3;
            end
            if isempty(pu)
                pu2='none';
            else
                pu2=pu;
            end
            qf=run_HiD_match(this.embedding(:,1:D), ...
                lbls, data, clusterIds, 'trainingNames', tNames, ...
                'matchStrategy', matchStrategy, 'log10', false, 'pu', pu2);
            qf.tClrs=clrs;
            [tCnt, sCnt, suprBestIdx4Clue, clusterMatch]=qf.getMatches;
            if Supervisors.CLUSTER_UNSUPERVISED
                unsupervisedClusterIdxs=find(clusterMatch>=unsprvStartLabel);
                if ~isempty(unsupervisedClusterIdxs)
                    suprBestIdx4Clue(unsupervisedClusterIdxs)=0;
                    clusterMatch(unsupervisedClusterIdxs)=0;
                    sCnt(unsupervisedClusterIdxs)=0;
                end
            end
            if Supervisors.DEBUG_QF
                [~, ~, tMi, sMi, supr1stIdx4Clue2, supr1stId4Clue, supervisorsUnmatched, ...
                    clustersUnmatched, tNoMatchIds, sNoMatchIds]=qf.getMatches2;
                qf.tNames(supervisorsUnmatched)
                assert(isequal(tNames(ismember(qf.tIds, tNoMatchIds)),...
                    qf.tNames(supervisorsUnmatched)));
                tNames(suprBestIdx4Clue(suprBestIdx4Clue>0))
                %assert(isequal(this.ids(supr1stIdx4Clue+1), supr1stId4Clue'))
                assert(isequal(suprBestIdx4Clue,supr1stIdx4Clue2))
                assert(isequal(clustersUnmatched, sCnt==0))
                assert(isequal(supervisorsUnmatched, tCnt==0))
                [tQ, sQ, tF, sF]=qf.getScores;
                [d2, unmatched]=qf.getTableData(clrs);
                tN=length(tQ);
                sN=length(sQ);
                NN=min([tN sN]);
                for i=1:NN
                    assert(isequal(num2str(d2{i,4}), num2str(tQ(i))))
                    assert(isequal(num2str(d2{tN+i,4}), num2str(sQ(i))))
                    assert(isequal(num2str(d2{i,5}), num2str(tF(i))))
                    assert(isequal(num2str(d2{tN+i,5}), num2str(sF(i))))
                    assert(d2{i, 8}==qf.tSizes(i))
                    assert(d2{tN+i, 8}==qf.sSizes(i))
                end
            end
            cluMdns=zeros(numClusters, D);
            clusterLabels=cell(1,numClusters);
            clusterNames=cell(1,numClusters);
            clusterColors=cell(1,numClusters);
            newSubsets=0;
            labels_=zeros(size(data, 1), 1);
            for i=1:numClusters
                l=clusterIds==i;
                if sCnt(i)==0
                    newSubsets=newSubsets+1;
                    clusterLabel=0-i;
                    clusterNames{i}=['New subset #' num2str(newSubsets) ];
                    clr=num2str(Supervisors.NewColor(newSubsets));
                else
                    clusterLabel=clusterMatch(i);
                    clusterNames{i}=tNames{suprBestIdx4Clue(i)};
                    clr=this.labelMap.get([num2str(clusterLabel) '.color']);
                    if Supervisors.DEBUG_QF
                        assert(isequal(clusterNames{i}, ...
                            this.labelMap.get(java.lang.String(num2str(clusterLabel)))))
                    end
                end
                clusterColors{i}=clr;
                clusterLabels{i}=clusterLabel;
                labels_(l)=clusterLabel;
                if Supervisors.VERBOSE
                    sum(l)
                end
                cluMdns(i,:)=median(data(l,:));
            end
            this.density.setLabels(labels_, clusterNames, ...
                clusterColors, cluMdns, clusterLabels);
            

            function d=normalize(d)
                mn=min(d);
                N_=length(mn);
                for j=1:N_
                    if mn(j)<=0
                        add_=1-mn(j);
                        d(:,j)=d(:,j)+add_;
                    end
                end
                for j=1:N_
                    mx_=max(d(:,j));
                    mn_=min(d(:,j));
                    r=mx_-mn_;
                    d(:,j)=(d(:,j)-mn_)/r;
                end
            end
        end
        
        function nnMatchWithClusters(this, data, density, ...
                numClusters, clusterIds, pu)
            if nargin<6
                pu=[];
            end
            this.density=density;
            if this.matchType==4 && isempty(this.nnUnreducedLabels)
                warning('matchType==4 without nnUnreducedLabes');
            else
                this.computeNearestNeighbors(data, pu);
            end
            clusterLabels=cell(1,numClusters);
            clusterNames=cell(1,numClusters);
            clusterColors=cell(1,numClusters);
            labels_=zeros(size(data, 1), 1);
            D=size(data,2);
            cluMdns=zeros(numClusters, D);
            newSubsets=0;
            for i=1:numClusters
                l=clusterIds==i;
                cluMdns(i,:)=median(data(l,:));
                clueLabels=this.nnLabels(l);
                u=unique(clueLabels);
                clusterLabelCnts=UmapUtil.discreteCount(clueLabels, u);
                [mxCnt, mxI]=max(clusterLabelCnts);
                lbl=u(mxI);
                if lbl==0
                    newSubsets=newSubsets+1;
                    lbl=0-i;
                    clusterNames{i}=['New subset #' num2str(newSubsets) ];
                    clusterColors{i}=num2str(...
                        Supervisors.NewColor(newSubsets));
                else
                    key=num2str(lbl);
                    clusterNames{i}=this.labelMap.get(java.lang.String(key));
                    clusterColors{i}=this.labelMap.get([key '.color']);
                end
                labels_(l)=lbl;
                clusterLabels{i}=lbl;
                if Supervisors.DEBUG_QF
                    fprintf(['%s (id=%d) has %d/%d events in cluster '...
                        '#%d''s %d events\n'], clusterNames{i}, ...
                        lbl, mxCnt, ...
                        sum(this.nnLabels==lbl), i, sum(l));
                end
            end
            this.density.setLabels(labels_, clusterNames, ...
                clusterColors, cluMdns, clusterLabels);
        end
                
        function matchWithClusters(this, data, density, ...
                numClusters, clusterIds)
            this.density=density;
            D=size(data, 2);
            if D>size(this.mdns,2)
                D=size(this.mdns,2);
                data=data(:, 1:D);
            end
            clusterLabels=cell(1,numClusters);
            clusterNames=cell(1,numClusters);
            clusterColors=cell(1,numClusters);
            avgs=this.mdns(:,1:D);
            if isprop(this, 'mads')
                devs=this.mads(:,1:D);
            else
                devs=[];
            end
            if any(this.cnts<Supervisors.MIN_FREQUENCY)
                avgs(this.cnts<Supervisors.MIN_FREQUENCY, 1)=this.xLimit(2)*5;
                avgs(this.cnts<Supervisors.MIN_FREQUENCY, 2)=this.yLimit(2)*5;
                if D>2
                    avgs(this.cnts<Supervisors.MIN_FREQUENCY, 3)=this.zLimit(2)*5;
                end
            end
            cluMdns=zeros(numClusters, D);
            cluDevs=zeros(numClusters, D);
            for i=1:numClusters
                l=clusterIds==i;
                if Supervisors.VERBOSE
                    sum(l)
                end
                cluMdns(i,:)=median(data(l,:));
                cluDevs(i,:)=mad(data(l,:), 1);
            end
            hasDevs=~isempty(devs);
            [D, I]=pdist2(avgs, cluMdns, 'euclidean', 'smallest', 1);
            labels_=zeros(size(data, 1), 1);
            reChecks={};
            for i=1:numClusters
                labelIdx=I(i);
                label=this.ids(labelIdx);
                if label==0
                    label=0-i;
                else
                    key=num2str(label);
                    clusterLabels{i}=label;
                    clusterNames{i}=this.labelMap.get(java.lang.String(key));
                    clusterColors{i}=this.labelMap.get([key '.color']);
                end
                l=clusterIds==i;
                if Supervisors.VERBOSE
                    sum(l)
                    this.labelMap.get(java.lang.String(num2str(label)))
                end
                if hasDevs
                    devDist=MatBasics.DevDist(cluMdns(i,:), cluDevs(i,:));
                    if any(D(i)>=devDist*Supervisors.DEV_UNIT_LIMIT)
                        reChecks{end+1}=struct('clustId', i, 'count',...
                            sum(l), 'label', label, 'labelIdx', labelIdx);
                        label=0-i;
                    end
                end
                labels_(l)=label;
            end
            if hasDevs
                N_=length(reChecks);
                while N_>0 
                    changes=[];
                    for i=1:N_
                        clustId=reChecks{i}.clustId;
                        label=reChecks{i}.label;
                        labelIdx=reChecks{i}.labelIdx;
                        closestLabelIdxs=labels_==label;
                        if Supervisors.VERBOSE
                            sum(closestLabelIdxs)
                            this.labelMap.get(java.lang.String(num2str(label)))
                            reChecks{i}
                        end
                        if any(closestLabelIdxs)
                            %Does this cluster with no label match
                            %   sit on the border of one with the 
                            %   closest label from the supervisor?
                            unlabeledClusterIdxs=clusterIds==clustId;
                            borderDistance=min(pdist2(...
                                data(unlabeledClusterIdxs, :), ...
                                data(closestLabelIdxs,:), 'euclidean', ...
                                'smallest', 1));
                            supervisorDevDistance=MatBasics.DevDist(...
                                avgs(labelIdx,:), devs(labelIdx,:));
                            limit=supervisorDevDistance*...
                                Supervisors.BORDER_DEV_UNIT_LIMIT;
                            if borderDistance<=limit
                                changes(end+1)=i;
                                labels_(unlabeledClusterIdxs)=label;
                            end
                        end
                    end
                    if isempty(changes)
                        break;
                    else
                        reChecks(changes)=[];
                    end
                    N_=length(reChecks);                    
                end
                newSubsetIds=unique(labels_(labels_<0));
                newSubsets=length(newSubsetIds);
                for i=1:newSubsets
                    clustId=0-newSubsetIds(i);
                    clusterLabels{clustId}=newSubsetIds(i);
                    clusterNames{clustId}=['New subset #' num2str(i)];
                    color_=Supervisors.NewColor(clustId);
                    if any(color_<0)
                        color_(color_<0)=0;
                    end
                    clusterColors{clustId}=num2str(color_);
                end
                N_=length(this.ids);
                remainder=data(labels_==0,:); 
                labels2=zeros(size(remainder,1),1);
                for i=1:N_
                    if ~any(find(labels_==this.ids(i),1))
                        label=this.ids(i);
                        if Supervisors.VERBOSE
                            this.labelMap.get(java.lang.String(num2str(label)))
                        end
                        [D2, I2]=pdist2(avgs(i,:), remainder, ...
                            'euclidean', 'smallest', 1);
                        pt=avgs(i,:)+devs(i,:);
                        devDist=pdist2(avgs(i,:), pt);
                        if any(D2<devDist*Supervisors.DEV_UNIT_LIMIT)
                            labels2(D2<devDist)=label;
                        end
                    end
                end
                if any(labels2)
                    labels_(labels_==0)=labels2;
                end
            end
            this.density.setLabels(labels_, clusterNames, ...
                clusterColors, cluMdns, clusterLabels);
        end
        
        function [lbls, names, clrs, next]=addUnsupervisedClusters(...
                this, lbls, names, clrs, pu)
            next=max(lbls)*2;
            un=lbls==0;
            if any(un)
                if isempty(this.unsupervisedClues)
                    if nargin<5
                        pu=[];
                    end
                    [~, clues]=Density.FindClusters(...
                        this.embedding(un,:), 'most high', 'dbm', pu);
                    clues(clues<0)=0;
                    u=unique(clues);
                    cnts=UmapUtil.discreteCount(clues, u)';
                    nClues=length(cnts);
                    for i=1:nClues
                        clue=u(i);
                        if clue>0
                            if cnts(i)<Supervisors.MIN_FREQUENCY
                                clues(clues==clue)=0;
                            end
                        end
                    end
                    this.unsupervisedClues=clues;
                else
                    clues=this.unsupervisedClues;
                end
                lbls(un)=next+clues;
                u=unique(clues);
                nClues=length(u);
                for i=1:nClues
                    clue=u(i);
                    if clue>0
                        names{end+1}=['Unsupervised cluster ID=' ...
                            num2str(next+clue)];
                        clrs(end+1,:)=[.25 .25 .25];
                    end
                end
                lbls(lbls==next)=0;
                
            end
        end
    end
    
    methods
        function changeMatchType(this, data, matchType, pu)
            if nargin<4
                pu=[];
            end
            this.matchType=matchType;
            if matchType==4 
                assert(size(data,2)==size(this.inputData,2), ...
                    sprintf('Training set unreduced is %d dimensions BUT test set is %d!!!',...
                    size(this.inputData,2), size(data,2)));
                this.computeNearestNeighborsUnreduced(data, pu);
            elseif matchType==3
                this.computeNearestNeighbors(data, pu);
            else
                this.computeAndMatchClusters(data, matchType, pu);
            end
        end
        
        function [numClusters, clusterIds, qf]=computeAndMatchClusters(...
                this, data, matchType, pu)
            if nargin<4
                pu=[];
                if nargin<3 
                    matchType=[];
                end
            end
            qf=[];
            if isempty(matchType)
                matchType=this.matchType;
            end
            [numClusters, clusterIds, dns]=this.findClusters(data, pu);
            if numClusters>0
                qf=this.matchClusters(data, dns, numClusters,  ...
                    clusterIds, matchType, pu);
            end
        end
        
        function setMatchType(this, matchType)
            this.matchType=matchType;
        end
        
        function qf=matchClusters(this, data, dns, numClusters, clusterIds, ...
                matchType, pu)
            if nargin<7
                pu=[];
                if nargin<6
                    matchType=this.matchType;
                else
                    this.matchType=matchType;
                end
            end
            qf=[];
            if matchType==0
                this.matchWithClusters(data, dns, numClusters, clusterIds);
            elseif matchType==1
                qf=this.qfMatchWithClusters(data, dns, numClusters, clusterIds, pu);
            elseif matchType>=2
                this.nnMatchWithClusters(data, dns, numClusters, clusterIds, pu);
            end
        end
        
        function nnLbls=prepareForTemplate(this)
            nnLbls=this.nnLabels;
            this.nnLabels=[];
        end
        
        function setNearestNeighborLabels(this, nnLabels)
            this.nnLabels=nnLabels;
        end
        
        function resetNearestNeighbors(this, data)
            this.nnLabels=[];
            this.computeNearestNeighbors(data);
        end
        
        function computeNearestNeighbors(this, data, pu)
            [R,D]=size(data);
            if length(this.nnLabels)~=R
                this.nnLabels=[];
            end
            if isempty(this.nnLabels)
                if R==size(this.embedding,1)
                    if isequal(this.embedding(:,1:D), data)
                        this.nnLabels=this.labels;
                        return;
                    end
                end
                txt=['Finding nearest neighbors in ' num2str(D) 'D space'];
                if nargin<3
                    pu=PopUp(txt, 'center', 'Supervising...', false);
                elseif ~isempty(pu)
                    old=pu.label.getText;                    
                    pu.label.setText(txt);
                end
                [~,II]=pdist2(this.embedding(:,1:D), data, 'euclidean', 'Smallest', 1);
                this.nnLabels=this.labels(II);
                if nargin<3
                    pu.close;
                elseif ~isempty(pu)
                    pu.label.setText(old);
                end
            end
        end
        
        function [labels, labelMap, reSupervised]=supervise(...
                this, data, doHtml, matchType)
            reSupervised=false;
            if nargin<4
                matchType=this.matchType;
                if nargin<3
                    doHtml=false;
                end
            else
                this.matchType=matchType;
            end
            if matchType==4 && isempty(this.nnUnreducedLabels)
                warning('matchType==4 without nnUnreducedLabes, setting matchType to 1');
                matchType=1;
            end
            if matchType==3 % nearest neighbor match alone (no cluster matching)
                reSupervised=length(this.nnLabels)~=size(data,1);
                this.computeNearestNeighbors(data);
                labels=this.nnLabels;
                labelMap=this.labelMap;
            elseif matchType==4
                labels=this.nnUnreducedLabels;
                labelMap=this.labelMap;
            else
                if isempty(this.density)
                    labels=[];
                    labelMap=[];
                    return;
                elseif length(this.density.labels) ~= size(data, 1)
                    reSupervised=true;
                    this.computeAndMatchClusters(data, matchType)
                end %0=median OR 2=nearest neighbor match to clusters
                labelMap=java.util.Properties;
                labels=this.density.labels;
                doMap(this.density.clusterColors, this.density.clusterNames);
            end
            
            function doMap(clusterColors, clusterNames)                
                ids_=unique(labels);
                N_=length(ids_);
                for i=1:N_
                    putInMap(ids_(i), clusterColors, clusterNames);
                end
            end
            
            function putInMap(id, clusterColors, clusterNames) 
                key=num2str(id);
                keyColor=[key '.color'];
                if id==0
                    if doHtml
                        name='<font color="62A162">unsupervised</font>';
                    else
                        name='\color[rgb]{0.4 0.65 0.4}\bf\itunsupervised';
                    end
                    color='92 92 128';
                elseif id<0
                    clustId=0-id;
                    nm=clusterNames{clustId};
                    if doHtml
                        name=['<font color="#4242BA"><i>' nm ' ?</i></font>'];
                    else
                        name=['\color[rgb]{0. 0.4 0.65}\bf\it' nm  ' ?'];
                    end
                    color=clusterColors(clustId);
                else
                    name=this.labelMap.get(java.lang.String(key));
                    if doHtml
                        if String.Contains(name, '^{')
                            name=strrep(name, '^{', '<sup>');
                            name=strrep(name, '}', '</sup>');
                        end
                    end
                    color=this.labelMap.get(keyColor);
                    if isempty(color)
                        color='215 205 225';
                    end
                end
                if isempty(name)
                    name=key;
                    if isempty(this.warnNoName)
                        this.warnNoName=java.util.HashSet;
                    end
                    if ~this.warnNoName.contains(java.lang.String(key))
                        this.warnings=this.warnings+1;
                        this.warnNoName.add(java.lang.String(key));
                        str=['There is no name for the data label "<b>' ...
                            key '"</b> ! <br>Hence we will use the '...
                            'label for the name...'];
                        warning(str);
                        if ~BasicMap.Global.is(Supervisors.PROP_STOP_WARNING)
                            bp=Gui.BorderPanel;
                            whine=['<html><center>' str '<br><br>'...
                                BasicMap.Global.smallStart '(<b><i>' ...
                                String.Pluralize2('name-less</i> label', ...
                                this.warnings) ' so far)' ...
                                BasicMap.Global.smallStart '<hr><br>'...
                                '</center></html>'];
                            bp.add(Gui.Label(whine), 'Center');
                            bp.add(Gui.CheckBox('STOP warning me ...', ...
                                false, BasicMap.Global, ...
                                Supervisors.PROP_STOP_WARNING), 'South');
                            msg(bp, 8, 'south east+', 'Name-less label...',...
                                'warning.png');
                        end
                    end
                end
                labelMap.put(java.lang.String(key), name);
                labelMap.put(keyColor, color);
            end
            
        end
        
        function [qf, qft]=qfTreeSupervisors(this, visible, pu, ttl)
            if nargin<4
                ttl='UMAP template''s training set';
                if nargin<3
                    pu=[];
                    if nargin<2
                        visible=true;
                    end
                end
            end
            if isempty(this.inputData)
                msg(Html.WrapHr(['Parameter reduction is out of date!<br>'...
                    '...redo with <u>this version</u> of the software!']));
                qf=[];
                qft=[];
                return;
            end
            [tNames, lbls, clrs]=this.getQfTraining;
            if isempty(pu)
                pu='none';
            end
            [qf, qft]=run_QfTree(this.inputData, lbls, {ttl}, ...
                'trainingNames', tNames, 'log10', true, 'colors', clrs, ...
                'pu', pu);
            if visible
                set(qft.fig, 'visible', 'on');
            end
        end
        
        function [ax, trainingSet, testSet, lbls]=plotTrainingAndTestSets(...
                this, data, ax, umap, pu, matchType, doJavaLegend)            
            if nargin<7
                doJavaLegend=false;
                if nargin<6
                    matchType=[];
                    if nargin<5
                        pu=[];
                        if nargin<4
                            ax=[];
                        end
                    end
                end
            end
            trainingSet=this.plotTrainingSet(umap, ax);
            if isempty(this.trainingSetPlots) || isempty(ax)
                fig=get(trainingSet.ax, 'Parent');
                p=get(fig, 'OuterPosition');
                if doJavaLegend || size(this.embedding, 2)>2
                    wFactor=.66;
                else
                    wFactor=.8;
                end
                if ispc
                    set(fig, 'OuterPosition', [p(1) p(2)-.29*p(4) ...
                        floor(p(3)*wFactor) p(4)+floor(p(4)*.43)]);
                else
                    set(fig, 'OuterPosition', [p(1) p(2) ...
                        floor(p(3)*wFactor) p(4)+floor(p(4)*.43)]);
                end
                Gui.FitFigToScreen(fig);
            end
            this.trainingSetPlots=trainingSet;
            
            if doJavaLegend
                subplot(2, 1, 1, trainingSet.ax)
            else
                subplot(4, 3, [1 2 4 5], trainingSet.ax)
            end
            if isempty(matchType)
                matchType=this.matchType;
            end
            title(trainingSet.ax, 'Training set');
            if doJavaLegend
                ax=subplot(2, 1, 2, ...
                    'Parent', get(trainingSet.ax, 'Parent'));
            else
                ax=subplot(4, 3, [7  8  10 11], ...
                    'Parent', get(trainingSet.ax, 'Parent'));
            end
            [testSet,lbls]=this.plotTestSet(umap, ax, data, ...
                pu, matchType, doJavaLegend);
            if isa(testSet.legendH, 'matlab.graphics.illustration.Legend')
                p=get(testSet.legendH, 'position');
                set(testSet.legendH, 'position',[1-p(3) .05 p(3) p(4)]);
            end
            ids_=unique(lbls);
            nSupervisors=sum(ids_>0);
            fig2=get(trainingSet.ax, 'Parent');
            set(fig2, 'name', ...
                [num2str(nSupervisors) ' supervisors found...']);
        end
        
        function [testSet, lbls,qfForClusterMatch]=plotTestSet(this, ...
                umap, ax, data, pu, matchType, doJavaLegend, finished)
            if doJavaLegend
                disp('javaLegend building`')
            end

            dns=this.density;
            mt=this.matchType;
            qfForClusterMatch=[];
            toD=size(data,2);
            if matchType==4 && isempty(this.nnUnreducedLabels)
                warning('matchType==4 without nnUnreducedLabes, setting matchType to 1');
                matchType=1;
            end
            nFound=0;
            found=[];
            if matchType<0
                good=false;
            else
                if matchType<3
                    [~,~,qfForClusterMatch]=...
                        this.computeAndMatchClusters(data, matchType, pu);
                    found=unique(this.density.labels);
                elseif matchType==3
                    this.computeNearestNeighbors(data, pu);
                    found=unique(this.nnLabels);
                else
                    this.computeNearestNeighborsUnreduced(data, pu);
                    found=unique(this.nnUnreducedLabels);
                end
                nFound=sum(found>0);
                good=nFound>0;
            end
            if good
                [lbls, lblMap]=this.supervise(data, ...
                    ~isempty(doJavaLegend) && doJavaLegend, matchType);
                if doJavaLegend 
                    LabelBasics.Frequency(this.labels, lblMap, true);
                    LabelBasics.Frequency(lbls, lblMap,  false);
                end
                testSet=Supervisors.Plot(data, lbls, lblMap, ...
                    this.fromD, umap, ax, doJavaLegend, false, ...
                    this.contourPercent);
                if ~isempty(this.trainingSetPlots)
                    testSet.addOtherPlots(this.trainingSetPlots);
                    title(ax, 'Test set');
                    if nargin>7 && finished
                        if doJavaLegend
                            subplot(2, 1, 1, this.trainingSetPlots.ax)
                            subplot(2, 1, 2, ax)
                        end
                    end
                end
            else
                lbls=[];
                if toD>2
                    Gui.PlotDensity3D(ax, data, 64, 'iso');
                else
                    ProbabilityDensity2.Draw(ax, data(:,1:2));
                end
                testSet.ax=ax;
                testSet.legendH=[];
            end
            if nargin<8 || ~finished
                this.density=dns;
                this.matchType=mt;
            end
            if isempty(this.trainingSetPlots)
                dimInfo=sprintf('  %dD\\rightarrow%dD', this.fromD, toD);
                xlabel(ax, ['UMAP-X' dimInfo]);
                ylabel(ax, ['UMAP-Y' dimInfo]);
                if toD>2
                    zlabel(ax, ['UMAP-Z' dimInfo]);
                end
            end
            nIds=length(this.ids);
            ttl=[String.encodeInteger(nFound) '/'  ...
                String.encodeInteger(nIds-1) ' subsets trained...'];
            set(get(ax, 'Parent'), 'name', ttl);
            if doJavaLegend
                if nFound<nIds-1
                    unfound={};
                    for i=1:nIds
                        label=this.ids(i);
                        if isempty(find(found==label, 1))
                            unfound{end+1}=char(this.labelMap.get(...
                                java.lang.String(num2str(label))));
                        end
                    end
                    html=['<html><b>Training subsets not found</b>:' ...
                        Html.ToList(unfound, 'ol') '<hr></html>'];
                    try
                        rr=testSet.javaLegend.getComponent(0);
                        rr.setToolTipText(html)
                    catch ex
                        ex.getReport
                    end
                end
                testSet.javaLegend.setTitle(ttl);
            end
        end
        
        function plots=plotTrainingSet(this, umap, ax)
            plots=Supervisors.Plot(this.embedding, this.labels,...
                 this.labelMap, this.fromD, umap, ax, [], false, ...
                    this.contourPercent);
        end
        
        function [qf, qft]=qfTreeSupervisees(this, embedding, rawData, ...
                visible, pu, ttl)
            if nargin<6
                ttl='"test set" found by UMAP supervised template';
                if nargin<5
                    pu=[];
                    if nargin<4
                        visible=true;
                    end
                end
            end
            [tNames, lbls, clrs]=this.getQfTrained(embedding);
            if ~isempty(lbls)
                if isempty(pu)
                    pu='none';
                end
                [qf, qft]=run_QfTree(rawData, lbls, {ttl}, ...
                    'pu', pu, 'trainingNames', tNames, ...
                    'log10', true, 'colors', clrs);
                if visible
                    set(qft.fig, 'visible', 'on');
                end
            else
                qf=[];
                qft=[];
            end
        end

        function [qf, qft]=qfDissimilarity(this, reducedData, ...
                unreducedData, visible, pu, priorQf, matchStrategy)
            if nargin<7
                matchStrategy=3;
                if nargin<6
                    priorQf=[];
                    if nargin<5
                        pu=[];
                        if nargin<4
                            visible=true;
                        end
                    end
                end
            end
            if isempty(this.inputData)
                msg(Html.WrapHr(['Parameter reduction is out of date!<br>'...
                    '...redo with <u>this version</u> of the software!']));
                return;
            end
            [tNames, tLbls, clrs]=this.getQfTraining;
            [sNames, sLbls]=this.getQfTrained(reducedData);
            if isempty(sLbls) %likely no clusters found by dbm or dbscan
                qft=[];
                qf=[];
                return;
            end
            qft=this.qfMatch(priorQf, visible, pu, 2, this.inputData,...
                unreducedData, tLbls, sLbls, tNames, sNames, clrs, ...
                false, matchStrategy);
            if ~isempty(qft)
                qf=qft.qf;
            else
                qf=[];
            end
        end
        
        function qft=qfDissimilarityTestSetPrior(this, reducedData, ...
                unreducedData, priorTestSetLbls, withTrainingSet,...
                visible, pu, file, matchStrategy)
            if nargin<9
                matchStrategy=1;
                if nargin<8
                    file=[];
                    if nargin<7
                        pu=[];
                        if nargin<6
                            visible=true;
                        end
                    end
                end
            end
            if withTrainingSet
                [tNames, tLbls, clrs]=this.getQfTraining;
                if isempty(this.inputData)
                    msg(Html.WrapHr(['Parameter reduction is out of date!<br>'...
                        '...redo with <u>this version</u> of the software!']));
                    return;
                end
                otherData=this.inputData;
                [sNames, sLbls]=this.getOtherTrained(priorTestSetLbls);
            else
                [sNames, sLbls]=this.getQfTrained(reducedData);
                otherData=unreducedData;
                [tNames, tLbls, clrs]=this.getOtherTrained(priorTestSetLbls);
            end
            if withTrainingSet
                scenario=1;
            else
                scenario=3;
            end
            qft=this.qfMatch(file, visible, pu, scenario, otherData, ...
                unreducedData, tLbls, sLbls, tNames, sNames, clrs, ...
                ~withTrainingSet, matchStrategy);
        end
        
        function qft=qfMatch(this, priorQf, visible, pu, scenario, ...
                tUnreducedData, sUnreducedData, tLbls, sLbls, ...
                tNames, sNames, clrs, fHist, matchStrategy)
            if isa(priorQf, 'QfHiDM')
                qf=priorQf;
                clrs=qf.tClrs;
            elseif exist(priorQf, 'file')
                qf=QfTable.Load(priorQf);                
            else
                if isempty(pu)
                    pu2='none';
                else
                    pu2=pu;
                end
                qf=run_HiD_match(tUnreducedData, tLbls,...
                    sUnreducedData, sLbls, 'trainingNames', tNames, ...
                    'matchStrategy', matchStrategy, 'log10', true, ...
                    'testNames', sNames, 'pu', pu2);
            end
            qft=QfTable(qf, clrs, [], get(0, 'currentFig'), visible);
            if ~qft.doHistQF(visible)
                qft=[];
                return;
            end
            if fHist || matchStrategy==2
                qft.doHistF(visible);
            end
            if ischar(priorQf) && ~exist(priorQf, 'file')
                qft.save(qf, priorQf);
            end
            rt=UMAP.REDUCTION_SUPERVISED_TEMPLATE;
            mt=UmapUtil.GetMatchTypeText(this.matchType, rt,...
                size(this.embedding,2), size(sUnreducedData, 2));
            sc=UmapUtil.GetMatchScenarioText(scenario, rt);
            qft.addSuffixToFigs([mt ': ' sc ' ' this.description]);
            qft.contextDescription=[mt ': ' sc];
            cntxt=this.context;
            cntxt.matchType=this.matchType;
            cntxt.matchScenario=scenario;
            cntxt.matchStrategy=matchStrategy;
            cntxt.reductionType=rt;
            cntxt.clusterDetail=this.clusterDetail;
            qft.context=cntxt;
        end
        
        function drawClusterBorders(this, ax)
            if isempty(this.density) % no clusters found?
                return;
            end
            wasHeld=ishold(ax);
            if ~wasHeld
                hold(ax, 'on');
            end
            N_=length(this.density.clusterColors);
            for i=1:N_
                if ~isempty(this.density.clusterColors{i})
                    clr=(str2num(this.density.clusterColors{i})/256)*.85;
                else
                    clr = [0 0 0];
                end
                gridEdge(this.density, true, i, clr, ax, .8, '.', '-', .5);
                if Supervisors.VERBOSE
                    disp(this.density.clusterColors{i});
                    disp(clr);
                    disp('ok');
                end
            end
            if ~wasHeld
                hold(ax, 'off');
            end
        end
        
        function initPlots(this, contourPercent)
            this.contourPercent=contourPercent;
        end
        
        function initClustering(this, detail, method2D, minopts, ...
                epsilon, distance)
            this.clusterDetail=detail;
            this.clusterMethod2D=method2D;
            this.minopts=minopts;
            this.dbscanDistance=distance;
            this.epsilon=epsilon;
        end
        
        function setClusterDetail(this, detail)
            this.clusterDetail=detail;
        end
        
        function [numClusts, clusterIds, dns]=findClusters(this, data, pu)
            if nargin<3
                pu=[];
            end
            if isempty(this.clusterDetail)
                this.clusterDetail='most high';
                this.clusterMethod2D='dbm';
                this.minopts=5;
                this.dbscanDistance='euclidean';
                this.epsilon=.6;
            end
            [mins, maxs]=Supervisors.GetMinsMaxs(data);
            [numClusts, clusterIds, dns]=Density.FindClusters(data, ...
                this.clusterDetail, this.clusterMethod2D, pu, ...
                this.epsilon, this.minopts, this.dbscanDistance, ...
                mins, maxs);
        end

    end
end