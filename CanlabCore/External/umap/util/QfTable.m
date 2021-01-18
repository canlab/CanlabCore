%   Class for Hi-D matching with merging of data subsets using 
%       QF match or F-measure or both 
%   This produces a visual table and histogram of dissimilarities and
%   overlaps between subset matches.
%
%   The function getAverages answers the question of overall goodness 
%   by returning the median/mean dissimilarity or overlap
%
%
%   QF Algorithm is described in 
%   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6586874/
%   and
%   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5818510/
%
%   Bioinformatics lead: Darya Orlova <dyorlova@gmail.com>
%   Software Developer: Stephen Meehan <swmeehan@stanford.edu> 
%
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef QfTable < handle
    properties
        contextDescription;
        context;
        fcnFpnSubsets;
    end
    
    properties(SetAccess=private)
        fig;
        data;
        unmatched;
        sortTable;
        tb;
        qf;
        tClrs;
        priorFig;
        fHistFig;
        qHistFig;
        falsePosNegFig;
        otherFigs={};
        R;
    end
    
    properties(Constant)
        PROP_OUTER_POS2='OuterPosition2_V2';
        DEBUG=false;
    end
    
    methods
        function this=QfTable(qf, tClrs, props, priorFig, visible)
            if nargin<5
                visible=true;
                if nargin<4
                    priorFig=get(0, 'currentFig');
                    if nargin<3
                        props=[];
                    end
                end
            end
            app=BasicMap.Global;
            if isempty(props)
                props=app;
            end
            this.tClrs=tClrs;
            this.qf=qf;
            PROP_COL_W='mdsColumnWidthsV1.33';
            PROP_COL_ORD='mdsColumnOrder_v2';
            PROP_SORT='mdsRowOrder_v2';
            path=BasicMap.Path;
            if visible
                pu=PopUp('Preparing QF dissimilarity table', 'north west+', ...
                    'Note...', false, true, fullfile(path, 'genieSearch.png'));
            else
                pu=[];
            end
            this.priorFig=priorFig;
            [this.fig, tb_]=Gui.Figure;
            this.tb=tb_;
            figName=['QF dissimilarity  ' ...
                num2str(length(qf.tIds)) ' X ' num2str(length(qf.sIds))...
                ' subsets'];
            set(this.fig, 'CloseRequestFcn', @(h, e)hush(h), ...
                'Name', figName);
            [this.data, labels, fmts, tips,  this.unmatched, groupIdx, ...
                freqIdx, rankIdx, symIdx]=QfTable.Contents(qf, tClrs, pu);
            [this.R, C]=size(this.data);
            J=edu.stanford.facs.swing.Basics;
            [sData, widths]=SortTable.Convert(this.data, fmts);
            this.sortTable=SortTable(this.fig, sData, ...
                labels, [], @(h,e)select(h,e), tips);
            st=this.sortTable;
            jt=st.jtable;
            N=length(widths);
            for i=1:N
                if app.highDef
                    factor=app.toolBarFactor;
                else
                    factor=1;
                end
                st.setColumnWidth(i, widths(i)*factor)
                st.setColumnWidth(i, widths(i)*factor)
            end
            preSorted=SortTable.SetRowOrder(jt, ...
                BasicMap.GetNumbers(props, PROP_SORT));
            if ~preSorted
                jt.sortColumn(rankIdx-1, true, true);
                jt.sortColumn(groupIdx-1, false, true);
                jt.sortColumn(freqIdx-1, false, false);
            end
            startingRowOrder=SortTable.GetRowOrder(jt);
            if ~isempty(qf.matrixHtml)
                ToolBarMethods.addButton(tb_, ...
                    fullfile(path, 'world_16.png'), 'See table in default `r', ...
                    @(h,e)browse(h));
            end
            bb=ToolBarMethods.addButton(tb_,fullfile(path, 'table.gif'), ...
                'Restore default row and column order', ...
                @(h,e)defaultOrder());
            ToolBarMethods.addButton(tb_,fullfile(path, 'leftArrowNarrow.png'), ...
                'Shift sorted columns to the left', ...
                @(h,e)shiftSortLeft());
            
            ToolBarMethods.addButton(tb_, fullfile(path, 'histQF.png'), ...
                'Open histogram of QF scores', ...
                @(h,e)doHistQF(this));            
            ToolBarMethods.addButton(tb_, fullfile(path, 'histF.png'), ...
                'Open histogram of F-measure scores', ...
                @(h,e)doHistF(this, true));
            savedOrder=BasicMap.GetNumbers(props, PROP_COL_ORD);
            if ~isempty(savedOrder)
                try
                    SortTable.SetColumnOrder(jt, savedOrder,...
                        BasicMap.GetNumbers(props, PROP_COL_W));
                catch
                end
            else
                defaultColumnOrder;
            end
            op=BasicMap.GetNumbers(props, QfTable.PROP_OUTER_POS2);
            if length(op)==4
                newFigPos=Gui.RepositionOnSameScreenIfRequired(op);
                set(this.fig, 'OuterPosition', newFigPos);
                Gui.FitFigToScreen(this.fig);
            elseif ~isempty(priorFig)
                Gui.SetToRight(this.priorFig, this.fig, false, 80);
            end
            if visible
                pu.close;
                Gui.SetFigVisible(this.fig);
            end
            if preSorted
                jt.unsort;
                SortTable.SetRowOrder(jt, ...
                    BasicMap.GetNumbers(props, PROP_SORT));
            elseif isempty(savedOrder)
                bb.doClick;
            end
            MatBasics.DoLater(@(h,e)heighten(), .31);
            this.qf.tClrs=this.tClrs;
            try
                if isempty(this.qf.falsePosNegs)
                    this.qf.getFalsePosNegRecords;
                end
                if ~isempty(this.qf.falsePosNegs)
                    try
                        this.qf.getFalsePosNegsTab(false);
                        this.qf.getFalsePosNegSubsets(true);
                    catch
                        %this.qf is a structure and not an instance
                        %so the 2 getFalse*() methods were called and 
                        %results are stored in structure
                    end
                    ToolBarMethods.addButton(tb_, ...
                        fullfile(path, 'plusMinus.png'), ...
                        'See false positive/negative results', ...
                        @(h,e)browseFalsePosNeg(this));
                    fpn=this.qf.falsePosNegSubsetsFile;
                    if ~isempty(fpn) && exist(fpn, 'file')
                        ToolBarMethods.addButton(tb_, ...
                            fullfile(path, 'tree.png'), ...
                            'Import false +/- into AutoGate GatingTree', ...
                            @(h,e)importFalsePosNegGates(this));
                    end
                end
            catch ex
                if ~isstruct(this.qf)
                    ex.getReport
                else
                    this.qf.falsePosNegs=[];
                end
            end
            function heighten
                drawnow;
                height=jt.getRowHeight;
                jt.setRowHeight(floor(height*1.33));
                try
                    if props.multiProps.is('mds_figHistQF', false)
                        this.doHistQF(visible);
                    end
                    if props.multiProps.is('mds_figHistF', false)
                        this.doHistF;
                    end
                catch ex
                    %disp('Using QF table outside of CytoGenie''s AutoGate')
                    %ex.getReport
                end
            end
            
            function shiftSortLeft(h)
                if ~SortTable.MoveSortColumnsLeft(jt)
                    if ismac
                        key='command';
                    else
                        key='control';
                    end
                    msg(Html.Wrap(['There is no sort order currently.'...
                        '<br>To sort the table you<ul>'...
                        '<li>Click on a column to sort it'...
                        '<li>Click again for descending order'...
                        '<li>Hold ' key ' key for multi column</ul>']), ...
                        8, 'north', 'No sort order...');
                else
                    r=jt.getSelectedRow;
                    if r<0
                        r=0;
                    end
                    rect=jt.getCellRect(r,0, true);
                    jt.scrollRectToVisible(rect);
                end
            end
            
            function defaultOrder()
                defaultRowOrder;
                defaultColumnOrder;
            end
            
            function defaultRowOrder
                SortTable.SetRowOrder(jt, [(rankIdx-1) true 0; ...
                    (groupIdx-1) false 0; (freqIdx-1) false 1]);
            end
            
            function defaultColumnOrder()
                idxs=1:jt.getColumnCount;
                idxs(rankIdx)=[];
                idxs(symIdx)=[];
                idxs=[rankIdx symIdx idxs]-1;
                SortTable.SetColumnOrder(jt, idxs, [], true);
            end
            
            function browse(h)
                Html.Browse(Html.Wrap([SortTable.ToHtml(jt) '<hr>' ...
                    Html.remove(qf.matrixHtml)]));
            end
            
            function hush(h)
                if ~isempty(priorFig) && ishandle(this.priorFig)
                    [columnOrder, widths]=SortTable.GetColumnOrder(jt);
                    props.set(PROP_COL_ORD, num2str(columnOrder));
                    props.set(PROP_COL_W, num2str(widths));
                    rowOrder=SortTable.GetRowOrder(jt);
                    if ~isequal(rowOrder, startingRowOrder)
                        props.set(PROP_SORT, MatBasics.Encode(...
                            rowOrder));
                    end
                    props.set(QfTable.PROP_OUTER_POS2, ...
                        num2str(get(h, 'OuterPosition')));
                    figure(this.priorFig);
                    Gui.CloseFigs(this.otherFigs);
                end
                delete(h);
            end
            
            function showSelected()
                selectedIdxs=SortTable.SelectedRows(jt);
                try
                    N_=length(selectedIdxs);
                    for ii=1:N_
                        
                    end
                catch ex
                    ex.getReport
                end
            end
            function select(h, e)
                if ~isempty(e.Indices)
                    rowIdx=e.Indices(1,1)-1;
                    if rowIdx>=0
                        try
                            colIdx=jt.convertColumnIndexToModel(...
                                e.Indices(1,2)-1);
                            st.showTip(colIdx+1);
                            rowIdx=jt.getActualRowAt(...
                                jt.convertRowIndexToModel(rowIdx))+1;
                            
                            showSelected;
                        catch ex
                            ex.getReport
                        end
                        return;
                    end
                end
                st.showTip(0);
            end
        end
        
        function [medianDissimilarity, meanDissimilarity, ...
                medianOverlap, meanOverlap]=getAverages(this)
            [~,medianDissimilarity, meanDissimilarity]=this.getData(true);
            [~,medianOverlap, meanOverlap]=this.getData(false);
        end
        
        function [qfData, mdn, mn]=getData(this, dissimilarity)
            if nargin<2 || dissimilarity
                idx=4;
            else
                idx=5;
            end
            qfData=cell2mat(this.data(:, idx))*100;
            qfData=qfData(~isnan(qfData));
            mdn=median(qfData);
            mn=mean(qfData);
        end
        
        function browseFalsePosNeg(this)
            if ~isempty(this.falsePosNegFig) ...
                    && ishandle(this.falsePosNegFig)
                figure(this.falsePosNegFig);
                return;
            end
            set(0, 'CurrentFigure', this.fig);
            [fig_,that]=FalsePositiveNegative.Plot([0 1], this.qf.falsePosNegFile);
            this.otherFigs{end+1}=fig_;
            figure(fig_);
            that.fcnImportGates=@()importFalsePosNegGates(this);
            that.fcnMoreHtml=@()getFalsePosNegMatrixHtml(this);
            Gui.Locate (fig_, this.fig, 'north east+');
            this.falsePosNegFig=fig_;
        end
        
        function html=getFalsePosNegMatrixHtml(this)
            html=FalsePositiveNegative.MatrixHtml(this.qf);
        end
        
        function importFalsePosNegGates(this)
            fpn=this.qf.falsePosNegSubsetsFile;
            if ~isempty(this.fcnFpnSubsets)
                feval(this.fcnFpnSubsets);
            else
                FalsePositiveNegative.ImportGates(fpn);
            end
        end

        function ok=doHistQF(this, visible)
            if nargin<2
                visible=true;
            end
            ok=true;
            if ishandle(this.qHistFig)
                figure(this.qHistFig);
            else
                fig_=Gui.NewFigure;
                if ~visible
                    set(fig_, 'visible', 'off');
                end
                this.qHistFig=fig_;
                Gui.Resize(fig_, .66);
                Gui.Locate(fig_, this.fig, 'south west+');
                ax=gca;
                [qfData, avgMdn, avgMn]=this.getData;
                if isempty(qfData)
                    ok=false;
                    return;
                end
                    
                histogram(ax, qfData, length(unique(qfData)));
                set(ax, 'units', 'normalized', 'Position', ...
                    [.1 .16 .8 .7]);
                Gui.StretchUpper(ax, @ylim, .1);
                Gui.StretchUpper(ax, @xlim, .1);
                xlabel(ax, '% dissimilarity', 'FontName', 'Arial')
                xlim(ax, [-5 100]);
                ylabel(ax, '# of subset matches', 'FontName', 'Arial')
                figName=['Dissimilarity ' ...
                    num2str(length(this.qf.tIds)) ' X ' ...
                    num2str(length(this.qf.sIds))];
                this.addTitleToHist(fig_, ax, ...
                    figName, avgMdn, avgMn);
                this.addFalsePosNegButton(fig_);
                %disp(num2str(get(fig_, 'OuterPosition')));
                set(fig_, 'CloseRequestFcn', @(h, e)hush(h));
            end
            
            function hush(h)
                if ishandle(this.fig)
                    figure(this.fig);
                end
                delete(h);
            end
        end
        
        function addTitleToHist(this, fig, ax, score, avgMdn, avgMn)
            ttl=[score ', median/mean=\color{blue}' ...
                String.encodeRounded(avgMdn, 1) ...
                '\color{black}/\color{blue}'...
                String.encodeRounded(avgMn, 1) '\color{black}% '];
            set(fig, 'name', String.RemoveTex(score), 'NumberTitle', 'off');
            if this.unmatched>0
                title(ax, {ttl, [ ...
                    num2str(this.R-this.unmatched) ' subsets matched, ' ...
                    '\color{red}' num2str(this.unmatched) ...
                    ' \color{black}NOT matched ']},'FontName', 'Arial');
            else
                title(ax, {ttl, [ ...
                    num2str(this.R) ' matches']}, ...
                    'FontName', 'Arial');
            end
        end
        
        function doHistF(this, visible)
            if nargin<3
                showFalsePosNeg=false;
                if nargin<2
                    visible=true;
                end
            end
            if ishandle(this.fHistFig)
                figure(this.fHistFig);
            else
                fig_=Gui.NewFigure;
                if ~visible
                    set(fig_, 'visible', 'off');
                end
                this.fHistFig=fig_;
                Gui.Resize(fig_, .66);
                Gui.Locate(fig_, this.fig, 'south east+');
                ax=axes('parent', fig_);
                [fData, avgMdn, avgMn]=this.getData(false);
                histogram(ax, fData, length(unique(fData)));
                
                xlabel(ax, '% overlap', 'FontName', 'Arial')
                xlim(ax, [-5 100]);
                ylabel(ax, '# of subset matches', 'FontName', 'Arial')
                set(ax, 'xlim', [0 100]);
                this.addTitleToHist(fig_, ax, ...
                    'Overlap ^{(F-measure)}', avgMdn, avgMn);
                set(ax, 'units', 'normalized', 'Position', ...
                    [.1 .16 .8 .7]);
                Gui.StretchUpper(ax, @ylim, .1);
                lbl='See false +/-';
                width=.2;
                this.addFalsePosNegButton(fig_);
                set(fig_, 'CloseRequestFcn', @(h, e)hush(h));
            end
            
            function hush(h)
                if ishandle(this.fig)
                    figure(this.fig);
                end
                delete(h);
            end
        end
        
        function addFalsePosNegButton(this, fig_, lbl, width)
            if nargin<4
                width=.2;
                if nargin<3
                    lbl='See false +/-';
                end
            end
            if isempty(this.qf.falsePosNegs)
                return;
            end
            uicontrol(fig_, 'style', 'pushbutton','String', lbl,...
                'Units', 'normalized', ...
                'FontWeight', 'bold', ...
                'ForegroundColor', 'blue',...
                'BackgroundColor', [1 1 .80],...
                'Position',[1-(width+.01), .015, width, .071],...
                'Callback', @(btn,event) browseFalsePosNeg(this));
        end
        
        function QF=save(this, qf, file) 
            QF.tIds=qf.tIds;
            QF.sIds=qf.sIds;
            %QF.matches=qf.matches;
            %QF.getStudNames=qf.getStudNames();
            QF.tNames=qf.tNames;
            QF.tSizes=qf.tSizes;
            QF.sSizes=qf.sSizes;
            QF.matrixHtml=qf.matrixHtml;
            QF.qfTableData=this.data;
            QF.unmatched=this.unmatched;
            QF.tClrs=this.tClrs;
            QF.falsePosNegs=qf.falsePosNegs;
            qf.relocateFalsePosNegFiles(file);
            QF.falsePosNegFile=qf.falsePosNegFile;
            QF.falsePosNegSubsetsFile=qf.falsePosNegSubsetsFile;
            QF.falsePosNegCnts=qf.falsePosNegCnts;
            QF.falsePosCulprits=qf.falsePosCulprits;
            QF.falseNegCulprits=qf.falseNegCulprits;
            QF.sNames=qf.sNames;
            if nargin>2 && ~isempty(file)
                save(file, 'QF');
            end
        end
        
        function addSuffixToFigs(this, suffix)
            renameFig(this.fig);
            renameFig(this.qHistFig);
            renameFig(this.fHistFig);
            function renameFig(fig)
                if ~isempty(fig)
                    set(fig, 'name', [ get(fig, 'name') ' ' suffix])
                end
            end
        end
        
        function closeFigs(this)
            Gui.CloseFig(this.fig);
            Gui.CloseFig(this.qHistFig);
            Gui.CloseFig(this.fHistFig);
        end

        function [tUnmatched, tN, sUnmatched, sN]=getMatchCounts(this)
            [tUnmatched, tN, sUnmatched, sN]=this.qf.getMatchCounts;
        end
    end
    
    methods(Static)
        
        
        function QF=Load(file)
            load(file, 'QF');
        end
        
        
        function [data, names, widths, tips, unmatched, plotIdx, freqIdx,...
                rankIdx, symIdx]=Contents(qf, clrs, pu)            
            unmatched=0;
            plotIdx=6;
            freqIdx=7;
            widths=[3 0; 17 nan; 5 0; 3 -2; 3 -2; ...
                    6 nan; 3 -1;  10 0; 5 nan; 8 0];
            names={'#', ...
                'Name', ...
                Html.WrapSm('Mat-<br>ches '), ...
                Html.WrapSm('QF dis-<br>similaritry'), ...
                Html.WrapSm('Overlap<br>F-measure'),...
                Html.WrapSm('<br>Trainer'), ...
                Html.WrapSm('<br>Freq.'), ...
                Html.WrapSm('# of <br>events'), ...
                '', ...
                Html.WrapSm('QF<br>Rank')};
            rankIdx=length(names);
            symIdx=rankIdx-1;
            tips={...
                'The subset''s item # ',...
                'The name of the subset including match name if applicable`',...
                'The number of subsets in this match set',...
                'The lower the % the more similar the match',...
                'The higher the % the more overlap (similarity) of cells (or bins)',...
                'yes for training set, no for test set ',...
                'The frequency of the subset within this selection of cells/events',...
                'The # of cells/events for the subset ',...
                'The subset''s frequency-sized symbol',...
                'The ranking of the match set'};
            try
                [data, unmatched]=qf.getTableData(clrs);
            catch ex
                if isstruct(qf)
                    data=qf.qfTableData;
                    unmatched=qf.unmatched;
                else
                    data=[];
                    unmatched=0;
                    ex.getReport
                end
            end
            if exist('pu', 'var') && ~isempty(pu)
                if ~isempty(pu.pb)
                    pu.pb.setValue(N);
                end
            end
        end

        function qft=RunWithGt(supervisors, gt, gid, fcs, fcsIdxs, ...
                data, file, embedding, visible, pu)
            qft=[];
            pFig=get(0, 'currentFig');
            if exist(file, 'file') && ~QfTable.DEBUG
                qf=QfTable.Load(file);    
                qft=QfTable(qf, qf.tClrs, [], pFig, visible);
                qft.doHistQF(visible);
                if ~isempty(qft.qf.falsePosNegs)
                    qft.doHistF(visible);
                end
                return;
            end
            if isempty(fcsIdxs)
                return;
            end
            newPu=nargin<10||isempty(pu);
            if newPu
                path=BasicMap.Path;
                pu=PopUp('Generating QF dissimilarity table', 'north', ...
                    'Note...', false, true, fullfile(path, 'genieSearch.png'));
                pu.setTimeSpentTic(tic);
            else
                pu.setText2('Generating QF dissimilarity table');
            end
            [gtData, ~, ~, gtIds]=QfHiDM.GetRequiredData2(fcs, fcsIdxs, ...
                gt, gid, pu, visible);
            if isempty(gtData)
                return;
            end
            if isequal(gtData, data)
                matchStrategy=2;
            else
                matchStrategy=1;
            end
            gids=MatBasics.GetUniqueIdsAndCntIfGt0(gtIds);
            hghl=gt.highlighter;
            N=length(gids);
            gtNames=cell(1,N);
            for i=1:N
                id=num2str(gids(i));
                gtNames{i}=gt.tp.getDescription(id);
            end
            [sNames, sLbls, clrs]=supervisors.getQfTrained(embedding);
            
            qf=run_HiD_match(gtData, gtIds, data, sLbls, ...
                'trainingNames', gtNames, ...
                'testNames', sNames, ...
                'matchStrategy', matchStrategy, 'log10', true, 'pu', pu);
            [~,~,tIdxForS, tIdForS]=qf.getMatches;
            gtClrs=zeros(N, 3);
            for i=1:N
                clrIdx=find(tIdxForS==i, 1, 'first');
                if clrIdx>0
                    gid=tIdForS(clrIdx);
                    gtClrs(i,:)=clrs(clrIdx,:);
                else
                    gid=gids(i);
                    gtClrs(i,:)=gt.highlighter.getColor(num2str(gid));
                end
            end
            if ~isempty(qf)
                qft=QfTable(qf, gtClrs, [], pFig, visible);
                qft.doHistQF(visible);
                if matchStrategy==2
                    qft.doHistF(visible);
                end
                if ~isempty(file) && ~exist(file, 'file')
                    qft.save(qf, file);
                end
            end
            if newPu
                pu.close;
            end
        end 
        
        function qft=Run(fcs, fcsIdxs, gt, gid1, gid2, file, visible, pu)
            pFig=get(0, 'currentFig');
            qft=[];
            if exist(file, 'file')
                qf=QfTable.Load(file);    
                qft=QfTable(qf, qf.tClrs, [], pFig, visible);
                qft.doHistQF(visible);
                return;
            end
            if isempty(fcsIdxs)
                return;
            end
            newPu=nargin<8|| isempty(pu);
            if newPu
                path=BasicMap.Path;
                pu=PopUp('Generating QF dissimilarity table', 'north', ...
                    'Note...', false, true, fullfile(path, 'genieSearch.png'));
                pu.setTimeSpentTic(tic);
            else
                pu.setText2('Generating QF dissimilarity table');
            end
            [data1, ids1, names1, clrs1]=go(gid1);
            if isempty(data1)
                return;
            end
            [data2, ids2, names2, clrs2]=go(gid2);
            if isempty(data2)
                return;
            end
            matchStrategy=1;
            qf=run_HiD_match(data1, ids1, data2, ids2, ...
                'trainingNames', names1, 'testNames', names2, ...
                'matchStrategy', matchStrategy, 'log10', true, ...
                'pu', pu);
            if ~isempty(qf)
                qft=QfTable(qf, clrs1, [], pFig, visible);
                qft.doHistQF(visible);
                if ~isempty(file) && ~exist(file, 'file')
                    qft.save(qf, file);
                end
            end
            if newPu
                pu.close;
            end
            
            function [gtData, gtIds, names, clrs]=go(gid)
                [gtData, ~, ~, gtIds]=QfHiDM.GetRequiredData2(fcs, ...
                    fcsIdxs, gt, gid, pu, visible);
                if isempty(gtData)
                    return;
                end
                gids=MatBasics.GetUniqueIdsAndCntIfGt0(gtIds);
                hghl=gt.highlighter;
                N=length(gids);
                clrs=zeros(N, 3);
                names=cell(1,N);
                for i=1:N
                    id=num2str(gids(i));
                    names{i}=gt.tp.getDescription(id);
                    clrs(i,:)=hghl.getColor(id, Gui.HslColor(i, N));
                end
            end
        end        
    end
end
