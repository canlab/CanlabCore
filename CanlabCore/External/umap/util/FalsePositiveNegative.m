%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Math Lead & Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University
%   License: BSD 3 clause
%

classdef FalsePositiveNegative < handle
    
    properties(Constant)
        FILE_EXTENSION='.fpnSubsets.mat';
        FILE_LABELS='sbOrLabels';
        FILE_NAMES='al';
    end
    
    properties
        fcnImportGates;
        fcnMoreHtml;
    end
    
    properties(SetAccess=private)
        map4PredictedClasses;
        map4Reclassifications;
        nClasses;
        nPredictedClasses;
        nPredicted;
        nNotPredicted;
        nReClassifications;
        nTestSets;
        t;
        tPredicted;
        tNotPredicted;
        file;
    end
    
    methods
        function this=FalsePositiveNegative(file)
            if nargin<1
                file='~/Documents/run_umap/examples/ustTest.pair/samusikImported_29D/matches_s2-3-4-5-6-7-8-9-10_t1_30nn_2D.txt';
            end
            this.file;
            file=strrep(file, '~', File.Home);
            if ~exist(file, 'file')
                [fldr,fn,ext]=fileparts(file);
                File.mkDir(fldr);
                [cancelled, bad]=WebDownload.Get({WebDownload.Url([fn ext], ...
                    UmapUtil.PATH)}, {file}, false, true);
                if cancelled || bad
                    msg('Could not get example');
                end
            end
            this.t=readtable(file, 'Delimiter', '\t');
            classes=unique(table2array(this.t(:, {'trainingClass'})));
            this.nClasses=length(classes);
            emptyRows=isnan(this.t.loD);
            unmatched=this.t.testSize==0 & ~emptyRows;
            this.nNotPredicted=sum(unmatched);
            this.tNotPredicted=this.t(unmatched,:);
            predicted=this.t.testSize>0&~emptyRows;
            this.nPredicted=sum(predicted);
            this.tPredicted=this.t(predicted,:);
            this.t=this.t(~emptyRows,:);
            this.map4Reclassifications=Map;
            predictedClasses=unique(table2array(this.tPredicted(:, {'trainingClass'})));
            allPredictedClasses=StringArray.List(predictedClasses);
            
            testSets=unique(table2array(this.tPredicted(:, {'testSet'})));
            this.nTestSets=length(testSets);
            [m, ~, rIdxs]=unique(table2array(this.tPredicted(1:end, ...
                {'reduction', 'matchType', 'clusterDetail', 'hiD', 'loD'})), 'rows');
            this.nReClassifications=size(m,1);
            for i=1:this.nReClassifications
                key=m(i,:);
                subt=this.tPredicted(rIdxs==i, :);
                reclassification=struct;
                reclassification.idx=i;
                reclassification.key=key;
                [matchType, clusterDetail]=...
                    FalsePositiveNegative.TranslateMatch(key(2), ...
                    key(3), key(4), key(5));
                labelDone=false;
                if this.nTestSets>1 
                    try
                        if isnumeric(subt.testSet) && ~any(isnan(subt.testSet))
                            reclassification.label=sprintf('%d; #%d %s %s', ...
                                subt.testSet(1), ...
                                key(1), matchType, clusterDetail);
                            labelDone=true;
                        elseif iscell(subt.testSet)
                            reclassification.label=sprintf('%s; #%d %s %s', ...
                                subt.testSet{1}, ...
                                key(1), matchType, clusterDetail);
                            labelDone=true;
                        end
                    catch ex
                    end
                end
                if ~labelDone
                    reclassification.label=sprintf('#%d %s %s', ...
                        key(1), matchType, clusterDetail);
                end
                x=table2array(subt(:, 'falsePosRatio'));
                y=table2array(subt(:, 'falseNegRatio'));
                z=table2array(subt(:, 'fMeasure'));
                predictionIdxs=StringArray.FirstListIndexes(...
                    allPredictedClasses,...
                    table2array(subt(:, 'trainingClass')));
                reclassification.xyz=[x y z predictionIdxs'];
                this.map4Reclassifications.set(num2str(i), reclassification);
            end
            cbn=ColorsByName;
            this.map4PredictedClasses=Map;
            app=BasicMap.Global;
            subStart=app.subStart;
            subEnd=app.subEnd;
            supStart=app.supStart;
            supEnd=app.supEnd;
            this.nPredictedClasses=length(predictedClasses);
            for i=1:this.nPredictedClasses
                match=struct;
                trainer=predictedClasses{i};
                if contains(trainer, supStart)
                    colorKeyNoTex=strrep(trainer, supStart, '');
                    colorKeyNoTex=strrep(colorKeyNoTex, supEnd, '');
                    match.color=cbn.get(colorKeyNoTex);
                else
                    match.color=cbn.get(trainer);
                end
                if isempty(match.color)
                    fprintf('No global color for "%s"\n', trainer);
                    match.color=Gui.HslColor(i,this.nPredictedClasses);
                end
                l=strcmp(this.tPredicted.trainingClass(:), trainer);
                match.reclassifications=rIdxs(l);
                assert(isequal(m(match.reclassifications,:),[this.tPredicted.reduction(find(l)) this.tPredicted.matchType(find(l)) this.tPredicted.clusterDetail(find(l)) this.tPredicted.hiD(find(l)) this.tPredicted.loD(find(l)) ]));
                match.matches=table2array(this.tPredicted(l, 'trainingSize'))-table2array(this.tPredicted(l, 'falseNeg'));
                match.falsePosRatios=table2array(this.tPredicted(l, 'falsePosRatio'));
                match.falseNegRatios=table2array(this.tPredicted(l, 'falseNegRatio'));
                match.fMeasures=table2array(this.tPredicted(l, 'fMeasure'));
                avgX=median(match.falsePosRatios);
                avgY=median(match.falseNegRatios);
                avgZ=median(match.fMeasures);
                
                txt=['(' num2str(length(match.fMeasures)) ', ' ...
                    String.encodePercent(avgX, 1, 0) '/' ...
                    String.encodePercent(avgY, 1, 0) '/' ...
                    String.encodePercent(avgZ, 1, 0) ')'];
                if contains(trainer, supStart)
                    lbl=strrep(trainer, supStart, '^{');
                    lbl=strrep(lbl, supEnd, '');
                    match.label=[lbl ' _{' txt '}'];
                else
                    match.label=[trainer ' _{' txt '}'];
                end
                match.html=[trainer '&nbsp;' subStart ...
                    '<b>' txt '</b>' subEnd];
                match.idx=i;
                this.map4PredictedClasses.set(trainer, match)
            end
        end
        
        function html=htmlDetails(this, fullDetails)
            N=size(this.t, 1);
            if N>0
                html='';
                td='<td align="right">';
                hdr=['<table cellpadding="2" border="1"><thead><tr>'...
                    '<th colspan="3">Predicted subsets</th>'...
                    '<th colspan="3">Trained subsets</th><th colspan="4">'...
                    'Prediction results</th></tr>'...
                    '<tr><th>Name</th><th>ID</th><th>Size</th>'...
                    '<th>Name</th><th>ID</th><th>Size</th>'...
                    '<th>F-measure</th><th>Precision<br>(false +)</th>'...
                    '<th>Recall<br>(false -)</th><th>QF Dis-'...
                    '<br>similarity</th></tr></thead>'];
                reduction=-1;
                groupMatchType=-1;
                groupClusterDetail=-1;
                breakHdr=['<table border="2" cellpadding="2"><thead><tr>'...
                    '<th>Reduction</th>'...
                    '<th>sampleSet</th>'...
                    '<th>trainingSet</th>'...
                    '<th>testSet</th>'...
                    '<th>neighbors</th>'...
                    '<th>hiD</th>'...
                    '<th>loD</th>'...
                    '<th>matchType</th>'...
                    '<th>clusterDetail</th>'...
                    '</tr></thead>'];
                if ~isempty(this.fcnMoreHtml)
                    moreHtml=feval(this.fcnMoreHtml);
                else
                    moreHtml='';
                end
                htmlIdx=1;
                for i=1:N
                    r=this.t(i, :);
                    if r.reduction ~= reduction ...
                            || r.matchType ~= groupMatchType ...
                            || r.clusterDetail ~= groupClusterDetail
                        reduction=r.reduction;
                        groupMatchType=r.matchType;
                        groupClusterDetail=r.clusterDetail;
                        if i==1
                            b='<hr>';
                        else
                            b='</table><hr>';
                        end
                        [matchType, clusterDetail]=...
                            FalsePositiveNegative.TranslateMatch(...
                            r.matchType, r.clusterDetail, r.hiD, r.loD);                        
                        breakBody=sprintf(['<tr>'...
                            '<td>%d</td>'...
                            '<td>%s</td>'...
                            '<td>%s</td>'...
                            '<td>%s</td>'...
                            '<td align="right">%d</td>'...
                            '<td align="right">%d</td>'...
                            '<td align="right">%d</td>'...
                            '<td>%s</td>'...
                            '<td>%s</td>'...
                            '</tr></table>'], reduction, r.sampleSet{1},...
                            r.trainingSet{1}, r.testSet{1}, r.neighbors, ...
                            r.hiD, r.loD, matchType, clusterDetail);
                        html=[html b breakHdr breakBody];
                        if iscell(moreHtml)
                            html=[html moreHtml{htmlIdx}];
                            htmlIdx=htmlIdx+1;
                        else
                            html=[html moreHtml];
                        end
                        if fullDetails
                            html=[html hdr];
                        end
                    end
                    if ~fullDetails
                        continue;
                    end
                    html=[html '<tr>'...
                        '<td>' r.trainingClass{1} '</td>'...
                        '<td>' num2str(r.trainingId) '</td>'...
                        td String.encodeInteger(r.trainingSize) '</td>'];
                    testIds=r.testIds{1};
                    if ~isempty(testIds)
                        testIds=r.testIds{1};
                        testIds=testIds(5:end);
                        testIdNums=str2num(testIds);
                        if isequal(r.trainingId, testIdNums)
                            testIds='';
                            testName='';
                        else
                            testIds=...
                                String.Num2Str(testIdNums, '<br>');
                            testName=['<small><b>'...
                                strrep(r.testClasses{1}, ', ', '<br>') ...
                                '</b></small>'];
                        end
                        html=[html ...
                            '<td>' testName '</td>'...
                            '<td>' testIds '</td>'...
                            td String.encodeInteger(r.testSize) '</td>'...
                            td  String.encodePercent(r.fMeasure, 1, 2) '</td>'...
                            td String.encodePercent(...
                            r.testSize-r.falsePos, r.testSize, 2) '<br>(' ...
                            String.encodePercent(r.falsePos, r.testSize, 2) ...
                            ')</td>'...
                            td String.encodePercent(...
                            r.trainingSize-r.falseNeg, r.trainingSize, 2) ...
                            '<br>(' ...
                            String.encodePercent(r.falseNeg, r.trainingSize,2) ...
                            ')</td>'...
                            td String.encodePercent(...
                            r.qfDissimilarity, 1, 2) '</td></tr>'];
                    else
                        html=[html '<td colspan="7">'];
                    end
                end
                html=[html '</table>'];
            end
        end
    end
    
    methods(Static)
        function [fig, this, legendH, classHtmls, ttlFound]=Plot(...
                plotType, fileNameOrObj, ax, legendType, className, ...
                testSetName, suppressTitle, suppressXTickLabel, ...
                drawLines)
            legendH=[];
            if nargin<9
                drawLines=true;
                if nargin<8
                    suppressXTickLabel=false;
                    if nargin<7
                        suppressTitle=false;
                        if nargin<6
                            testSetName='sample';
                            if nargin<5
                                className='subset';
                                if nargin<4
                                    legendType=2; % window with checkboxes
                                    if nargin<3
                                        ax=[];
                                        if nargin<2
                                            fileNameOrObj=[];
                                            if nargin<1
                                                plotType=[0 1]; %3D of +/- & fm
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if length(plotType)==1 && (plotType(1)<0||plotType(1)>4)
                plotType=3; %2D of +/- 
            end
            if isempty(fileNameOrObj)
                fileNameOrObj=fullfile(...
                    UmapUtil.LocalSamplesFolder,...
                        'matches_s2-3-4-5-6-7-8-9-10_t1_30nn_2D.txt');
            end
            newFigNeeded=isempty(ax);
            if newFigNeeded
                priorFig=get(0, 'CurrentFigure');
                [ax, fig]=Gui.GetOrCreateAxes(figure('visible', 'off'));
            else
                fig=get(ax, 'Parent');
                priorFig=[];
            end
            if length(plotType)>1
                plotType=unique(plotType);
                if ~all(plotType>=0 & plotType<5)
                    warning('Only D values of 0-3 are valid')
                    plotType=plotType(plotType>=0 & plotType<4);
                end
                secondPlotIs1D=ismember(plotType(2), [1 2]);
                N=length(plotType);
                if N>4
                    R=2;
                    C=3;
                    thirdPlotIs1D=ismember(plotType(3), [1 2]);
                elseif N>2
                    R=2;
                    C=2;
                else
                    R=2;
                    C=1;
                end
                if newFigNeeded
                    ip=get(fig, 'position');
                    if N>2
                        ip(3)=ip(3)*1.5;
                        ip(4)=ip(4)*1.4;
                    else
                        ip(3)=ip(3)*.9;
                        ip(4)=ip(4)*1.6;
                    end
                    set(fig, ...
                        'name', 'UST results computing!',...
                        'NumberTitle', 'off', 'position', ip);
                else
                    set(fig, 'visible', 'on');
                end
                if ~isempty(priorFig)
                    set(0, 'CurrentFigure', priorFig);
                end
                pu=PopUp('Building false positive/negative plots',  ...
                    'center', 'Computing...',...
                    true, false, 'genieSearch.png');
                if ~isempty(priorFig)
                    set(0, 'CurrentFigure', fig);
                end
                
                pu.initProgress(N*2);
                didLegend3=false;
                for i=1:N
                    if N>4
                        if i==2 && secondPlotIs1D
                            p=4;
                        elseif i==3 && thirdPlotIs1D
                            p=5;
                        elseif i==4
                            p=2;
                        elseif i==5
                            p=3;
                        else
                            p=i;
                        end
                    elseif N>2
                        if i==2 
                            if secondPlotIs1D
                                p=3;
                            else
                                p=2;
                            end
                        elseif i==3 
                            if secondPlotIs1D
                                p=2;
                            else
                                p=3;
                            end
                        else
                            p=i;
                        end
                    else
                        p=i;
                    end
                    ax2=subplot(R, C, p, 'Parent', fig);
                    drawnow;
                    pu.incrementProgress;
                    if plotType(i)>2
                        if didLegend3
                            lt=0;
                        else
                            if N>4 
                                if i==5
                                    lt=1;
                                    didLegend3=true;
                                else
                                    lt=0;
                                end
                            else
                                lt=1;
                                didLegend3=true;
                            end
                        end
                    elseif i==1
                        lt=1;
                    else
                        lt=0;
                    end
                    [~, this, legendH, classHtmls, ttlFound]=...
                        FalsePositiveNegative.Plot(plotType(i), ...
                        fileNameOrObj, ax2, lt, ...
                        className, testSetName, i>1,...
                        i==1 && secondPlotIs1D);
                    fileNameOrObj=this;
                    drawnow;
                    if i==1 && ismember(plotType(1), [0 1]) ...
                            && secondPlotIs1D
                        pos=get(ax2, 'position');
                        set(ax2, 'position', [pos(1) pos(2)-(.15*pos(4)) pos(3) pos(4)*1.15])
                    elseif i==2 && secondPlotIs1D
                        pos=get(ax2, 'position');
                        set(ax2, 'position', [pos(1) pos(2) pos(3) pos(4)*1.25])
                    end
                    pu.incrementProgress;
                end
                if newFigNeeded
                    set(fig, 'name', ['UST found ' ttlFound '!']);
                    doMenus;
                    Gui.SetFigVisible(fig);
                end
                pu.close;
                return;
            end
            if ischar(fileNameOrObj)
                this=FalsePositiveNegative(fileNameOrObj);
            else
                this=fileNameOrObj;
            end
            is1D=plotType<3;
            keys=this.map4Reclassifications.keys;
            N=length(keys);
            lineClrs=zeros(N,3);
            for i=1:N
                lineClrs(i,:)=Gui.HslColor(i, N);
            end
            keys=this.map4PredictedClasses.keys;
            N=length(keys);
            assert(N==this.nPredictedClasses);
            freqs=zeros(1,N);
            if plotType>2
                Hs=zeros(1,N);
                htmls=cell(1,N);
                labels=cell(1,N);
            else
                htmls={};
                labels={};
                Hs=[];
            end
            for i=1:N
                key=keys{i};
                d=this.map4PredictedClasses.get(key);
                freqs(i)=sum(d.matches);
                if plotType>2
                    htmls{i}=d.html;
                    labels{i}=d.label;
                end
            end
            cla(ax, 'reset');
            hold(ax, 'on');
            mx=max(freqs);
            if plotType>2
                classHtmls=containers.Map(keys, htmls);
                xlim(ax, [-.1 1.1]);
            else
                classHtmls=[];
                xlim(ax, [0 this.nPredictedClasses+1]);
            end
            ylim(ax, [-.1 1.1]);
            if plotType==4
                zlim(ax, [-.1 1.1]);
            end
            labels1D=cell(1,N);
            app=BasicMap.Global;
            supStart=app.supStart;
            supEnd=app.supEnd;
            
            classHsPerReclassification=cell(this.nReClassifications,1);
            for i=1:N
                key=keys{i};
                d=this.map4PredictedClasses.get(key);
                sz=20+floor(35*freqs(i)/mx);
                clr=d.color;
                N2=length(d.falsePosRatios);
                if plotType==4
                    Hs(i)=plot3(ax, d.falsePosRatios, d.falseNegRatios, ...
                        d.fMeasures, 'linestyle', 'none', 'marker', '.', ...
                        'markerSize', sz, 'MarkerFaceColor', clr, ...
                        'MarkerEdgeColor', clr);
                elseif plotType==3
                    Hs(i)=plot(ax, d.falsePosRatios, d.falseNegRatios, ...
                        'linestyle', 'none', 'marker', '.', ...
                        'markerSize', sz, 'MarkerFaceColor', clr, ...
                        'MarkerEdgeColor', clr);
                else
                    perc=String.encodePercent(freqs(i), mx,0);
                    if contains(key, supStart)
                        lbl=strrep(key, supStart, '^{');
                        lbl=strrep(lbl, supEnd, '}');
                        labels1D{d.idx}=[lbl ' ^{\color{blue}\bf' perc '}'];
                    else
                        labels1D{d.idx}=[key ' ^{\color{blue}\bf' perc '}'];
                    end
                end
                for j=1:N2
                    reclassificationIdx=d.reclassifications(j);
                    x=d.falsePosRatios(j);
                    y=d.falseNegRatios(j);
                    f=d.fMeasures(j);
                    if plotType==4 && legendType==2
                        H=plot3(ax, x, y, f, ...
                            'linestyle', 'none', 'marker', '.', ...
                            'markerSize', sz, 'MarkerFaceColor', clr, ...
                            'MarkerEdgeColor', clr, 'Visible', 'off');
                    elseif plotType==3 && legendType==2
                        H=plot(ax, x, y, ...
                            'linestyle', 'none', 'marker', '.', ...
                            'markerSize', sz, 'MarkerFaceColor', clr, ...
                            'MarkerEdgeColor', clr, 'Visible', 'off');
                    elseif plotType<3
                        clIdx=d.idx;
                        if ~drawLines
                            clr=lineClrs(reclassificationIdx,:);
                        end
                        if plotType==2
                            H=plot(ax, clIdx, f, ...
                                'linestyle', 'none', 'marker', '.', ...
                                'markerSize', sz, 'MarkerFaceColor', clr, ...
                                'MarkerEdgeColor', clr);
                        elseif plotType==1
                            H=plot(ax, clIdx, y, ...
                                'linestyle', 'none', 'marker', '.', ...
                                'markerSize', sz, 'MarkerFaceColor', clr, ...
                                'MarkerEdgeColor', clr);
                        else
                            H=plot(ax, clIdx, x, ...
                                'linestyle', 'none', 'marker', '.', ...
                                'markerSize', sz, 'MarkerFaceColor', clr, ...
                                'MarkerEdgeColor', clr);
                        end
                    else
                        continue;
                    end
                    Hs2=classHsPerReclassification{reclassificationIdx};
                    Hs2(end+1)=H;
                    classHsPerReclassification{reclassificationIdx}=Hs2;
                end
            end
            uistack(Hs, 'top');
            if plotType==4
                view(ax, [1 1 1]);
            end
            styles={'-', '--', '-.'};
            widths=[.3 .4 .51 .7]; 
            keys=this.map4Reclassifications.keys;
            N=length(keys);
            assert(this.nReClassifications==N);
            if legendType==2
                otherKeys=zeros(1,N);
                otherValues=cell(1,N);
            end
            if plotType<3
                freqs=[];
            end
            for i=1:N
                key=keys{i};
                d=this.map4Reclassifications.get(key);
                clr=lineClrs(i,:);
                subsets4Line=classHsPerReclassification{d.idx};
                ls=styles{MatBasics.Mod(i,length(styles))};
                lw=widths(MatBasics.Mod(i,length(widths)));
                xyz=d.xyz;
                if plotType>2
                    idxs=MatBasics.OrderByCloseness(xyz(:, 1:2), [1 1]);
                    xyz=xyz(idxs,:);                    
                    if plotType==4
                        H=plot3(ax, xyz(:,1), xyz(:,2), xyz(:,3),...
                            'LineWidth', lw, 'LineStyle', ls,  'visible', 'off', ...
                            'marker', '.', 'markerSize', 2, 'Color', clr, ...
                            'MarkerEdgeColor', clr);
                    elseif plotType==3
                        H=plot(ax, xyz(:,1), xyz(:,2),...
                            'LineWidth', lw, 'LineStyle', ls, 'visible', 'off', ...
                            'marker', '.', 'markerSize', 2, 'Color', clr,...
                            'MarkerEdgeColor', clr);
                    end
                elseif drawLines
                    clIdxs=xyz(:,4);
                    [~,idxs]=sort(clIdxs);
                    xyz=xyz(idxs,:);
                    clIdxs=xyz(:,4);
                    if plotType==2
                        H=plot(ax, clIdxs, xyz(:,3),...
                            'LineWidth', lw, 'LineStyle', ls, ...
                            'marker', '.', 'markerSize', 2, 'Color', clr,...
                            'MarkerEdgeColor', clr);
                    elseif plotType==1
                        H=plot(ax, clIdxs, xyz(:,2),...
                            'LineWidth', lw, 'LineStyle', ls, ...
                            'marker', '.', 'markerSize', 2, 'Color', clr,...
                            'MarkerEdgeColor', clr);
                    else
                        H=plot(ax, clIdxs, xyz(:,1),...
                            'LineWidth', lw, 'LineStyle', ls, ...
                            'marker', '.', 'markerSize', 2, 'Color', clr,...
                            'MarkerEdgeColor', clr);
                    end
                end
                if legendType>0 && legendType<3
                    if plotType>2
                        htmls{end+1}=[testSetName ' ' d.label];
                        labels{end+1}=[testSetName ' ' d.label];
                    else
                        htmls{end+1}=d.label;
                        labels{end+1}=d.label;
                    end
                    freqs(end+1)=nan;
                    if drawLines
                        Hs(end+1)=H;
                    else
                        Hs(end+1)=subsets4Line(1);
                    end
                end
                if legendType==2
                    otherKeys(i)=Hs(end);
                    otherValues{i}=subsets4Line;
                    uistack(otherValues{i}, 'top')
                end
            end
            if plotType==2
                tickLabels={'\bf\color{red}Epic fail', ...
                '\bf\color{magenta}0.25', '0.5', ...
                '\bf\color[rgb]{.21 .21 .7}0.75', ...
                '\bf\color{blue}Perfect'};
            else
                tickLabels={'\bf\color{blue}Perfect', ...
                '\bf\color[rgb]{.21 .21 .7}25%', '50%', ...
                '\bf\color{magenta}75%', ...
                '\bf\color{red}Epic fail'};
            end
            ticks=[0 .25 .5 .75 1];
            if plotType>2
                if plotType==4
                    xlabel(ax, 'X: False positive rate');
                    ylabel(ax, 'Y: False negative rate');
                    zlabel(ax, 'Z:  F-measure');
                    tickLabels2={'\bf\color{red}No overlap', ...
                        '\bf\color[rgb]{.21 .21 .7}25%', '50%', ...
                        '\bf\color{magenta}75%', ...
                        '\bf\color{blue}100% overlap'};
                    set(ax, 'xtick', ticks, 'xTickLabel', tickLabels, ...
                        'xTickLabelRotation', -25, ...
                        'ytick', ticks, 'yTickLabel', tickLabels, ...
                        'yTickLabelRotation', 25,...
                        'ztick', ticks, ...
                        'zTickLabel', tickLabels2, ...
                        'zTickLabelRotation', -25);
                elseif plotType==3
                    set(ax, 'xtick', ticks, 'xTickLabel', tickLabels, ...
                        'xTickLabelRotation', -25, ...
                        'ytick', ticks, 'yTickLabel', tickLabels, ...
                        'yTickLabelRotation', -25);
                    xlabel(ax, 'False positive rate');
                    ylabel(ax, 'False negative rate');
                end
                grid(ax, 'on')
            else
                if ~suppressXTickLabel
                    xlabel(ax, String.Pluralize2(className, this.nPredictedClasses));
                    set(ax, 'xtick', 1:this.nPredictedClasses, ...
                        'xTickLabel', labels1D,...
                        'xTickLabelRotation', 45);
                else
                    set(ax, 'xtick', 1:this.nPredictedClasses, ...
                        'xTickLabel', [])
                end
                set(ax, 'ytick', ticks, ...
                    'yTickLabel', tickLabels, ...
                    'yTickLabelRotation', -25);
                if plotType==2
                    ylabel(ax, 'F-measure ');
                elseif plotType==1
                    ylabel(ax, 'False negative rate');
                else
                    ylabel(ax, 'False positive rate');
                end           
                set(ax, 'xGrid', 'on')
                set(ax, 'yGrid', 'on')
            end
            if plotType==4            
                view(ax, [1 1 1]);
            end
            avgFp=String.encodePercent(median(this.tPredicted.falsePosRatio), 1, 0);
            avgFn=String.encodePercent(median(this.tPredicted.falseNegRatio), 1, 0);
            avgFm=String.encodeRounded(median(this.tPredicted.fMeasure), 2);
            ttlFound=[num2str(this.nPredictedClasses) ...
                '/' String.Pluralize2(className, this.nClasses) ' in '...
                String.Pluralize2(testSetName, this.nTestSets) ];
            ttlAvgs=[', avg f+/f-/fm is ' avgFp '/' avgFn '/' avgFm];
            ttlPrefix=[ttlFound ttlAvgs];
            ttlSuffix=[num2str(this.nPredicted) '/' ...
                String.Pluralize2('prediction', ...
                this.nReClassifications*this.nClasses) ...
                ' were made'];
            if this.nReClassifications>this.nTestSets
                ttlSuffix=[ ttlSuffix ' in '...
                    String.Pluralize2('reclassification', ...
                    this.nReClassifications)]; 
            end
            if ~suppressTitle
                title(ax, {['UST:  ' ttlPrefix], ttlSuffix})
            end
            %set(ax3D, 'plotboxaspectratio', [1 1 1])    
            if newFigNeeded
                set(fig, 'visible', 'on', ...
                    'name', ['UST found ' ttlFound '!'],...
                    'NumberTitle', 'off');
                menus=doMenus(ax);
            end
            if legendType>0 && length(Hs)>1
                if legendType==2
                    labels=htmls;
                end
                [legendH,~,~,~,plots]=Plots.Legend(Hs, labels, [], ...
                    [], [], true, freqs, [], [], legendType==2);
                if plotType>2
                    lTtl=Gui.UST_LEGEND1;
                else
                    lTtl=Gui.UST_LEGEND2;
                end
                if legendType==2 
                    legendH.setTitle(lTtl);
                    plots.otherPlotMap=containers.Map(otherKeys, otherValues);
                else
                    title(legendH, lTtl);
                end
            end
            function c=doMenus(newAx, um)
                if nargin<2
                    um=uimenu(fig, 'Label', 'More +/-');
                    if nargin<1
                        newAx=[];
                    end
                    FalsePositiveNegative.MoreBtn(fig, @(h,e)popupMenu(h, newAx));
                end
                c={};
                c{end+1}=uimenu(um, 'Label', 'False positives', ...
                     'Callback', @(h,e)changeD(0, newAx));
                c{end+1}=uimenu(um, 'Label', 'False negatives', ...
                     'Callback', @(h,e)changeD(1, newAx));
                 c{end+1}=uimenu(um, 'Label', 'F-measure (harmonic mean of pos/neg)', ...
                     'Callback', @(h,e)changeD(2, newAx));                
                c{end+1}=uimenu(um, 'Label', 'False pos x neg', ...
                     'Callback', @(h,e)changeD(3, newAx));
                c{end+1}=uimenu(um, 'Label', 'False pos x neg x f-measure',...
                     'Callback', @(h,e)changeD(4, newAx));
                if length(plotType)==1
                    c{end+1}=uimenu(um, 'Separator', 'on', 'Label', ...
                        'No legend', 'Callback', @(h,e)changeLegend(0));
                    c{end+1}=uimenu(um,  'Label', 'Legend inside', ...
                        'Callback', @(h,e)changeLegend(1));
                    c{end+1}=uimenu(um, 'Label', 'Legend outside', ...
                        'Callback', @(h,e)changeLegend(2));
                    set(c{plotType+1}, 'Enable', 'off');
                    set(c{6+legendType}, 'Enable', 'off');
                end
            end
            function popupMenu(h, newAx)
                jm=PopUp.Menu;
                app=BasicMap.Global;
                c={};
                c{end+1}=Gui.NewMenuItem(jm, 'False positives', @(h,e)changeD(0, newAx));
                c{end+1}=Gui.NewMenuItem(jm, 'False negatives', ...
                    @(h,e)changeD(1, newAx));
                c{end+1}=Gui.NewMenuItem(jm, ['<html>F-measure' app.supStart ...
                    ' (harmonic mean of pos/neg)' app.supEnd '</html>'], ...
                    @(h,e)changeD(2, newAx));
                c{end+1}=Gui.NewMenuItem(jm, 'False pos x neg', ...
                    @(h,e)changeD(3, newAx));
                c{end+1}=Gui.NewMenuItem(jm, 'False pos x neg x f-measure',...
                    @(h,e)changeD(4, newAx));
                if length(plotType)==1
                    jm.addSeparator;
                    c{end+1}=Gui.NewMenuItem(jm, 'No legend', ...
                        @(h,e)changeLegend(0));
                    c{end+1}=Gui.NewMenuItem(jm, 'Legend inside', ...
                        @(h,e)changeLegend(1));
                    c{end+1}=Gui.NewMenuItem(jm, 'Legend outside', ...
                        @(h,e)changeLegend(2));
                    mi=c{plotType+1};
                    mi.setEnabled(false);
                    mi=c{6+legendType};
                    mi.setEnabled(false);
                    jm.addSeparator;
                    doLines=~drawLines;
                    if doLines
                        word='Draw';
                    else
                        word='Remove';
                    end
                    c{end+1}=Gui.NewMenuItem(jm, ...
                        [word ' line for reclassification runs'],...
                        @(h,e)toggleLine());
                end
                jm.addSeparator;
                Gui.NewMenuItem(jm, 'See details in browser', ...
                    @(h,e)browseDetails(true, false));
                Gui.NewMenuItem(jm, 'See full details in browser', ...
                    @(h,e)browseDetails(true, true));
                if ~isempty(this.fcnImportGates)
                    jm.addSeparator;
                    Gui.NewMenuItem(jm, 'Import into AutoGate GatingTree', ...
                        @(h,e)doImport, 'tree.png');
                end
                jm.show(h, 15, 15);    
                
                function doImport
                    feval(this.fcnImportGates);
                end
            end
            
            function browseDetails(figToo, fullDetails)
                html=this.htmlDetails(fullDetails);
                if figToo
                    img=Html.TempImg(fig);
                    html=[img '<hr>' html '<hr>' img];
                end
                
                Html.BrowseString(Html.Wrap(html));
            end
            
            function toggleLine
                if isa(legendH, 'javax.swing.JDialog')
                    legendH.dispose;
                elseif ~isempty(legendH)
                    delete(legendH);
                end
                drawLines=~drawLines;
                [~,~,legendH]=FalsePositiveNegative.Plot(plotType, ...
                    this, ax, legendType, className, testSetName, ...
                    false, false, drawLines);
            end
            
            function changeLegend(lt)
                if isa(legendH, 'javax.swing.JDialog')
                    legendH.dispose;
                elseif ~isempty(legendH)
                    delete(legendH);
                end
                if lt==0
                    legendH=[];
                else
                    [~,~,legendH]=FalsePositiveNegative.Plot(plotType,...
                        this, ax, lt, className, testSetName, ...
                        false, false, drawLines);
                end
                for k=0:2
                    if lt==k
                        en='off';
                    else
                        en='on';
                    end
                    set(menus{6+k}, 'Enable', en);
                end
                legendType=lt;
            end
            
            function changeD(newPlotType, newAx)
                if isempty(newAx)
                    FalsePositiveNegative.Plot(newPlotType, this, [], 2, ...
                        className, testSetName);
                else
                    if isa(legendH, 'javax.swing.JDialog')
                        legendH.dispose;
                    end
                    if newPlotType>2
                        drawLines=true;
                    end
                    [~,~,legendH]=FalsePositiveNegative.Plot(...
                        newPlotType, this, newAx, legendType, ...
                        className, testSetName, false, false, drawLines);
                    plotType=newPlotType;
                    for k=0:4
                        if newPlotType==k
                            en='off';
                        else
                            en='on';
                        end
                        set(menus{k+1}, 'Enable', en);
                    end
                end
            end
        end
        
        function [H, J]=MoreBtn(fig, callback)
            if nargin<2
                callback=@(h,e)disp('No callback for this');
            end
            app=BasicMap.Global;
            if app.highDef
                tip=['<html><center>' app.smallStart ...
                    'See additional false positive <br>and'...
                    ' negative plots<hr>' Html.ImgXy(...
                    'wayneMoore1.png', [], 1.5)...
                    app.smallEnd '<hr>Wayne Moore</center></html>'];
                heightFactor=.75;
            else
                tip=['<html><center>See additional false positive '...
                    'and<br>negative plots<hr>' Html.ImgXy(...
                    'wayneMoore1.png', [], .91)...
                    '<hr>Wayne Moore</center></html>'];
                heightFactor=1;
            end
            [H,J]=Gui.ImageLabel('  Weighing more...', 'plusMinus.png', ...
                tip, callback, fig, 1, 1, true);
            J.setBackground(java.awt.Color(1, 1, .7))
            if BasicMap.Global.highDef
                p=get(H, 'position');
                set(H, 'position', [p(1) p(2) p(3)*.5 p(4) * heightFactor]);
            end
        end
        
        function ff = PlotPrecisionRecall(precisions, recalls, labels, varargin)
            %   AUTHORSHIP
            %   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
            %   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
            %   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
            %   Provided by the Herzenberg Lab at Stanford University
            %   License: BSD 3 clause
            %
            if nargin<1
                ff=FalsePositiveNegative.PlotPrecisionRecall(...
                    [.92 .5 .1], [.984 .8 .05], ...
                    {'Basophils', 'T-cells', 'B-cells'}', 'sizes', [42 24 83],...
                    'colors', [1 1 0;.5 .5 0;0 0 1], 'invert', false);
                msg(Html.WrapHr(['<h2>Incorrect arguments</h2>'...
                    'This is an example of a FalsePositiveNegative.PlotPrecisionRecall...<br><br>'...
                    '(See the top of the file function  in FalsePositiveNegative for this...)']),...
                    0,'south east++', 'Incorrect arguments!');
                return
            end
            if nargin < 3 || isempty(labels)
                labels = 1:length(precisions);
                warning('Labels were not given for the precision-recall scatter plot!');
            end
            
            p=parseArguments();
            parse(p,varargin{:});
            args=p.Results;
            invert = args.invert;
            if args.visible
                visibility = 'on';
            else
                visibility = 'off';
            end
            if max(args.sizes)>75 %normalize sizes for use with markers
                mx=max(args.sizes);
                args.sizes=25+(args.sizes/mx*50);
                nSizes=length(args.sizes);
                if nSizes>1 % sort biggest first for legend sake
                    [~,I]=sort(args.sizes, 'descend');
                    labels=labels(I);
                    precisions=precisions(I);
                    recalls=recalls(I);
                    if size(args.colors,1)==nSizes
                        args.colors=args.colors(I,:);
                    end
                    args.sizes=args.sizes(I);
                end
            end
            if invert
                xName = 'False positive rate';
                yName = 'False negative rate';
                precisions = 1 - precisions;
                recalls = 1 - recalls;
            else
                xName = 'Precision';
                yName = 'Recall';
            end
            
            ff=figure('Name', [yName ' vs. ' xName],'visible', visibility);
            ax=axes('Parent', ff);
            
            try
                if isempty(args.colors)
                    cbn=ColorsByName;
                    N=length(labels);
                    args.colors=zeros(N,3);
                    for i=1:N
                        clr=cbn.get(labels{i});
                        if isempty(clr)
                            clr=Gui.HslColor(i,N);
                        end
                        args.colors(i,:)=clr;
                    end
                end
                gscatter(ax, precisions, recalls, labels,args.colors,[],args.sizes);
            catch ex
                %r2019a and earlier do not support ax argument
                gscatter(precisions, recalls, labels,args.colors,[],args.sizes);
            end
            xlim(ax, [-.1 1.1]);
            ylim(ax, [-.1 1.1]);
            xlabel(ax, xName);
            ylabel(ax, yName);
            tickLabels={'\bf\color{blue}Perfect', ...
                '\bf\color[rgb]{.21 .21 .7}25%', '50%', ...
                '\bf\color{magenta}75%', ...
                '\bf\color{red}Failed'};
            ticks=[0 .25 .5 .75 1];
            if ~invert
                tickLabels=flip(tickLabels);
            end
            set(ax, 'xtick', ticks, ...
                'xTickLabel', tickLabels, 'xTickLabelRotation', -25, ...
                'ytick', ticks, ...
                'yTickLabel', tickLabels, 'yTickLabelRotation', -25);
            
            function p=parseArguments()
                p = inputParser;
                addParameter(p,'invert',false,@(x) islogical(x));
                addParameter(p,'visible',true,@(x) islogical(x));
                addParameter(p, 'sizes', 25, @(x) all(x>=0) );
                addParameter(p, 'colors', {}, @(x) isempty(x) || isnumeric(x) && size(x, 2)==3);
            end
        end
        
        function [matchType, clusterDetail]=TranslateMatch(...
                matchType, clusterDetail, hiD, loD)
            if matchType==4
                matchType=['nn' num2str(hiD) 'D'];
            elseif matchType==3
                matchType=['nn' num2str(loD) 'D'];
            elseif matchType==2
                matchType=['nnClust'];
            else
                matchType=['disClust'];
            end
            try
                if matchType<3
                    clusterDetail=Density.DETAILS{clusterDetail};
                else
                    clusterDetail='';
                end
            catch
                clusterDetail='?';
            end
        end
        
        function head=TabHead
            head=sprintf(['reduction\tsampleSet\t'...
                'trainingSet\ttestSet\t'...
                'neighbors\thiD\tloD\tmatchType\tclusterDetail'...
                '\ttrainingClass\ttrainingId\ttrainingSize\t'...
                'testClasses\ttestIds\ttestSize\tfMeasure\tfalsePos\t'...
                'falsePosRatio\tfalseNeg\tfalseNegRatio\t'...
                'qfDissimilarity\n']);
            
        end
        
        function rowHead=TabRowHead(reduction, sampleSet, ...
                trainingSet, testSet, n_neighbors, hiD, ...
                n_components, matchType, clusterDetail)
            rowHead=sprintf(...
                '%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%e\t', ...
                reduction, sampleSet, trainingSet,...
                testSet, n_neighbors, hiD, ...
                n_components, matchType, clusterDetail);
        end
        
        function [row, found]=TabRow(rowHead, rec)
            found=rec.testSize>0;
            if found
                ratioFalsePos=String.encodeRounded(...
                    rec.falsePos/rec.testSize, 3, false, ...
                    [], false);
                ratioFalseNeg=String.encodeRounded(...
                    rec.falseNeg/rec.trainingSize, 3, false, ...
                    [], false);
            else
                ratioFalsePos='0';
                ratioFalseNeg='0';
            end
            
            row=[sprintf(['%s%s\t%d\t%d\t%s\t%s\t%d\t'...
                '%s\t%d\t%s\t%d\t%s\t%s'], rowHead, ...
                rec.trainingClass, rec.trainingId, rec.trainingSize, ...
                StringArray.toString(rec.testClasses), ...
                ['ids ' num2str(rec.testIds)], rec.testSize, ...
                String.encodeRounded(rec.fMeasure,3, ...
                false, [], false), ...
                rec.falsePos, ratioFalsePos,...
                rec.falseNeg, ratioFalseNeg,...
                String.encodeRounded(rec.qfDissimilarity,3,...
                false,[],false)) newline];
        end
        
        function [body, notFound]=TabRows(recs, reduction, sampleSet, ...
                trainingSet, testSet, n_neighbors, hiD, ...
                n_components, matchType, clusterDetail)
            body='';
            notFound='';
            N2=length(recs);
            if N2==0
                return;
            end    
            notFound='';
            context=sprintf(...
                '%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t', ...
                reduction, sampleSet, trainingSet,...
                testSet, n_neighbors, hiD, ...
                n_components, matchType, clusterDetail);
            for j=1:N2
                [line, found]=...
                    FalsePositiveNegative.TabRow(context, recs{j});
                if found
                    body=[body line];
                else
                    notFound=[notFound line];
                end
            end
            if ~isempty(notFound)
                body=[body newline notFound];
            end
            body=[body newline];
        end
        
        function SaveAndPlot(file, head, body, plots)
            if nargin<4
                plots=[0 1];
            end
            if ~isempty(file)
            File.SaveTextFile(file, [head body])
            fig=FalsePositiveNegative.Plot([0 1], file);
            end
        end
        
       function [out, precisions, recalls]=CreateRecords(...
               tIds, tSizes, tNames, sIds, sSizes, sNames, ...
               trainingEvents, trainedTestEvents, trainingMap, ...
               dissimilarities, fMeasures)
            out={};
            N=length(tIds);
            precisions = zeros(1,N);
            recalls = zeros(1,N);
            C=size(trainingEvents, 2);
            for i=1:N
                trainingId=tIds(i);
                for c=1:C
                    if c~=C && ~any(trainingEvents(:,c)==trainingId)
                        continue;
                    end                    
                    r=struct;
                    r.trainingClass=String.ToHtmlSupFromTex(tNames{i});
                    r.trainingId=trainingId;
                    r.trainingSize=tSizes(i);
                    r.truth=sum(trainingEvents(:,c)==trainingId & trainedTestEvents==trainingId);
                    r.falsePos=sum(trainingEvents(:,c)~=trainingId & trainedTestEvents==trainingId);
                    r.falseNeg=sum(trainingEvents(:,c)==trainingId & trainedTestEvents~=trainingId);
                    precisions(i) = r.truth/(r.truth + r.falsePos);
                    recalls(i) = r.truth/(r.truth + r.falseNeg);
                    sIdxs=trainingMap.getAll(num2str(trainingId));
                    nMatches=length(sIdxs);
                    if nMatches==1
                        sIdx=sIdxs{1};
                        sName=String.ToHtml(sNames{sIdx});
                        r.testIds=sIds(sIdx);
                        r.testSize=sSizes(sIdx);
                        r.testClasses={String.ToHtmlSupFromTex(sName)};
                        r.fMeasure=fMeasures(sIdx);
                        r.qfDissimilarity=dissimilarities(sIdx);
                    elseif nMatches>1
                        testSize=0;
                        testClasses=cell(1, nMatches);
                        testIds=zeros(1, nMatches);
                        for j=1:nMatches
                            sIdx=sIdxs{j};
                            sName=String.ToHtml(sNames{sIdx});
                            testIds(j)=sIds(sIdx);
                            testSize=testSize+sSizes(sIdx);
                            testClasses{j}={String.ToHtmlSupFromTex(sName)};
                        end
                        r.testIds=testIds;
                        r.testSize=testSize;
                        r.testClasses=testClasses;
                        r.fMeasure=fMeasures(sIdx);
                        r.qfDissimilarity=dissimilarities(sIdx);
                    else
                        r.testSize=0;
                        r.testIds=0;
                        r.testClasses={};
                        r.fMeasure=-1;
                        r.qfDissimilarity=-1;
                    end
                    out{end+1}=r;
                    break;
                end
            end
       end
       
       function [ok, fldr, fl]=ImportGates(f)
           longTtl=['<html><center>Save gate (AKA subset) definitions '...
               '<br>for importing into AutoGate</center></html>'];
           [fldr, fl]=File.PutFile(fullfile(File.Home, 'Documents'), f, ...
               BasicMap.Global, 'fpnSubsets', ...
               'Save AutoGate import file', longTtl, ...
               FalsePositiveNegative.FILE_EXTENSION);
           ok=false;
           if ~isempty(fldr)
               try
                   [ok, msgTxt]=copyfile(f, fullfile(fldr, fl), 'f');
                   if ~ok
                       msgError(Html.H2Small2('Could not save!', msgTxt));
                   end
               catch ex
                   msgError(Html.H2Small2('Could not save!', ex.message));
               end
           end
       end
       
              
       function html=ToHtml(matrix, columnLabels, rowLabels, notes1, notes2)
           [rows,cols]=size(matrix);
           if nargin<4
               noteCols=0;
               nNotes=0;
               if nargin<3
                   rowLabels=[];
                   if nargin<2
                       columnLabels=[];
                   end
               end
           else
               nNotes=length(notes1);
               noteCols=2;
               assert(nNotes==length(notes2));
           end
           if isempty(columnLabels)
               columnLabels=cell(1,cols);
               for col=1:cols
                   columnLabels{col}=num2str(col);
               end
           end
           if isempty(rowLabels)
               rowLabels=columnLabels;
           end
           blue='<font color="blue">';
           legend=['<i>Rows are ' blue 'predicted classes</font> (e.g. manual gates), ' ...
               'columns are ' blue ' predicting classes </font> (e.g. automated gates from EPP, UST, X-shift etc).' ...
               '<br><b>Bold</b> indicates maximum row intersection for column unless it occurs in the row labeled "no subset"'...
               '<br><font color="red">Red</font> indicates a non maximum row intersection '...
               'for the column or a maximum in the "no subset" row.</b><i>'];
           html=['<h2>Classification predictions </h2>'...
               legend '<table cellpadding="3" border="1"><tr><th></th>'];
           app=BasicMap.Global;
           for col=1:cols+noteCols
               html=[html '<th width=40px><font color="blue">' ...
                   app.smallStart columnLabels{col} app.smallEnd ...
                   '</font></th>'];
           end
           html=[html '</tr>'];
           mx=max(matrix(2:end-2, 2:end-2));
           mx(end+1)=matrix(end-1, end-1);
           for row=1:rows
               html=[html '<tr><td><font color="blue">' rowLabels{row} '</font></td>'];
               for col=1:cols
                   v=matrix(row,col);
                   if v==0
                       html=[html '<td></td>'];
                       continue;
                   end
                   num=String.encodeInteger(v);
                   decorate=col>1 && row>1 && col<cols && row < rows;
                   if decorate
                       if v==mx(col-1)
                           if row<rows-1 && col==cols-1
                               html=[html '<td align="right"><font '...
                                   'color="red"><b><i>' num ...
                                   '</i></b></td>'];
                           else
                               html=[html '<td align="right"><b>' ...
                                   num '</b></td>'];
                           end
                       elseif v>0
                           html=[html '<td align="right"><font '...
                               'color="red"><b><i>' num '</i></b></td>'];
                       else
                           html=[html '<td align="right">' num '</td>'];
                       end
                   else
                       if col==1 || row==1
                           f='<font color="blue">';
                       else
                           f='';
                       end
                       html=[html '<td align="right">' f  num '</td>'];
                   end
               end
               if nargin>3
                   if row>nNotes || row==1
                       html=[html '<td colspan="2"></td>'];
                   else
                       html=[html '<td>' notes1{row-1} '</td><td>' notes2{row-1} '</td>'];
                   end
               end
               html=[html '</tr>'];
           end
           html=[html '</table>'];
       end
       
       function html=MatrixHtml(qf)
           html=FalsePositiveNegative.ToHtml(qf.falsePosNegCnts, ...
               ['size', qf.sNames, 'no subset', 'overlap','false -', 'false +'], ...
               ['size', qf.tNames, 'no subset', 'overlap'], ...
               qf.falseNegCulprits, qf.falsePosCulprits);
       end
       

    end
end