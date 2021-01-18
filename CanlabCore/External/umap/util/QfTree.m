%   Class for  building a dendrogram based on QF dissimilarity measurements
%
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
%   QF dissimilarity is described in 
%   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6586874/
%
%   QF Tree (phenogram) is described in 
%   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6586874/

classdef QfTree < handle
    properties(Constant)
        SHOW_MARKER_HTML=false;
        PROP_WND ='phenoWndV3';
        PROP_IDS='phenogramNodeIDS';
        PROP_D='phenogramNodeD';
        FONT_NAME='Arial';
        FONT_SIZE_NAME='\fontsize{14}';
        FONT_SIZE_IDS='\fontsize{12}';
        FONT_SIZE_FREQ='\fontsize{11}';
        FONT_SIZE_QF_SCORE='\fontsize{11}';
        PROP_FONT='QfTreeV4';
        PROP_FREQ_BASIS='QfTreeFreq';
        PROP_DISTANCE='phenogramDistance';
        DFLT_DISTANCE=9;
        DISTANCES={'QF', 'CityBlock', 'Chebychev',...
                        'Cosine', 'Euclidean', 'SquaredEuclidean', ...
                        'Earth mover''s (EMD)', 'QF + Euclidean', ...
                        'QF + CityBlock', 'Fast EMD'};
    end
    
    properties(SetAccess=private)
        heatMapArgs;
        props;
        tb;
        wndProp;
        cascadePos=1;
        eppImg=Html.Img('epp.png',[],1,true) ;
        idxById;
        doubleClicker;
        selected=[];
        lineWidths;
        edgeColors;
        sm1;
        sm2;
        fig;
        ax;
        nodes;
        xyNormalized;
        qf;
        gt;
        gid;
        app;
        fg;
        figs1D;
        leafNames;
        btnGenie;
        jBtnGenie;
        nodeNames;
        lastNear=0;
        mnuOpts;
        gtMap;
        xyIn;
        jdMrker;
        textHs;
        plotHs;
        showMeasurements;
        folder;
        fcsIdxs;
        fcsIdxsStr;
        lastTipFnc;
        nNodes=0;
        nLeafs=0;
        btnSvg;
        figures;
        tipIsNew=false;
        sup1;
        sup2;
        originalNames;
        isEpp;
        ttl;
    end
    
    methods
        function this=QfTree(qf, ttl, props, propPfx, colors, ...
                edgeColors, lineWidths, tNames)
            nT=length(tNames);
            this.originalNames=tNames;
            if isempty(propPfx)
                propPfx=[QfTree.PROP_WND '.' propPfx '.' num2str(nT)];
            else
                propPfx=[QfTree.PROP_WND  '.' num2str(nT)];
            end
            app=BasicMap.Global;
            for i=1:nT
                tNames{i}=app.getMatchName(tNames{i}, false, '^{', '}');
            end
            this.app=app;
            if app.highDef
                this.sup1='<b>(';
                this.sup2=')</b>';
            else
                this.sup1='<sup>(';
                this.sup2=')</sup>';
            end
            this.qf=qf;
            this.nNodes=length(this.qf.nodeSzs);
            this.props=props;
            [this.fig, this.ax, this.nodes, this.xyNormalized, ...
                this.figs1D, this.nodeNames, this.xyIn,...
                this.textHs, this.plotHs, this.wndProp]=...
                QfTree.New(qf, ttl, props, propPfx, colors, edgeColors, ...
                lineWidths, tNames, this);
            set(this.fig, 'UserData', this);
            N=length(this.plotHs);
            this.edgeColors=cell(1, N);
            this.lineWidths=zeros(1,N);
            for i=1:N 
                this.edgeColors{i}=get(this.plotHs(i), 'markerEdgeColor');
                this.lineWidths(i)=get(this.plotHs(i), 'lineWidth');
            end
            this.selected=java.util.LinkedHashSet;
            this.leafNames=tNames;  
            this.nLeafs=length(this.leafNames);
            try
                CytoGate.setHelp('AG_Phenogram');
            catch
            end
            this.ttl=ttl;
        end           
        
        function idx=getIdx(this, id)
            if isempty(this.idxById)
                this.idxById=java.util.TreeMap;
                N=length(this.leafNames);
                for i=1:N
                    id_=this.nodes{i};
                    id_=id_{1};
                    this.idxById.put(id_, i);
                end
            end
            if this.idxById.containsKey(id)
                idx=this.idxById.get(id);
            else
                if ~isempty(this.gtMap) && ~isempty(this.gtMap{end})
                    nLeaves=length(this.leafNames);
                    N=length(this.gtMap);
                    for i=1:N
                        this.idxById.put(this.gtMap{i}, i+nLeaves);
                    end
                end
                if this.idxById.containsKey(id)
                    idx=this.idxById.get(id);
                else
                    idx=0;
                end
            end
        end
        
        function setMouse(this)
            set(this.fig,'WindowButtonMotionFcn', @(h,e)motion(this));
            this.doubleClicker=DoubleClicker(this.fig, this.ax);
            this.doubleClicker.setFnc([], @(h, cp, handles)wbu(this, cp), [], ...
                @(h, cp, handles)wdcu(this, cp), []);
        end
        
        function [filePath, file]=getMainPngPath(this)
            file=['qft_v3_' this.gid '.png'];
            filePath=fullfile(this.folder, file);
        end
        
        function img=getMainPngImg(this, scale, forBrowser)
            if nargin<3
                forBrowser=false;
                if nargin<2
                    scale=.66;
                end
            end
            [~, file]=this.getMainPngPath;
            img=Html.ImgXy(file, this.folder, scale, forBrowser);
        end
        
        function saveMainPng(this)
            [newAx, newFig]=Gui.CopyAxes(this.fig, false);
            op=get(this.fig, 'outerPosition');
            nLeafs_=this.nLeafs;
            MX=25;
            if nLeafs_>MX
                hFactor=1.33+((nLeafs_-MX)*.025);
                if hFactor>2.5
                    hFactor=2.5;
                end
            else 
                hFactor=1.33;
            end
            set(newFig, 'visible', 'off', 'outerposition', ...
                [op(1) op(2) op(3)*hFactor op(4)*hFactor]);
            set(newAx, 'title', []);
            if nLeafs_<=MX
                Gui.UpdateFontSize(newFig, QfTree.PROP_FONT, ...
                    BasicMap.Global, 2);
            end
            saveas(newFig, this.getMainPngPath);
            if nLeafs_<=MX
                Gui.UpdateFontSize(newFig, QfTree.PROP_FONT, ...
                    BasicMap.Global,-2);
            end
        end
        
        function setGt(this, gt, gid)
            this.sm1=gt.app.smallStart;
            this.sm2=gt.app.smallEnd;
            this.gt=gt;
            this.isEpp=gt.tp.isEpp(gid);
            this.gid=gid;
            this.setMouse;
            tb_=ToolBar.New(this.fig, true, false, false, false);
            this.tb=tb_;
            this.btnSvg=ToolBarMethods.addButton(tb_, 'svg.png', 'Create SVG File',...
                @(h,e)svg(this));
            ToolBarMethods.addButton(tb_,'wrench.png', ...
                'See options', @(h,e)contextMenu(this, 0, h))                        
            ToolBarMethods.addButton(tb_,'heatMap.png', ...
                'See heat map of median marker epression', ...
                @(h,e)showMarkerHeatMap(this, false));
            ToolBarMethods.addButton(tb_,'table.gif', ...
                'See selected/all 1D PathFinders in browser', ...
                @(h,e)browseAll(this, false));
            ToolBarMethods.addButton(tb_,'phenogram.png', ...
                'See selected/all 1D PathFinders in browser', ...
                @(h,e)browseAll(this, true));
            ToolBarMethods.addButton(tb_,'pseudoBarHi.png', ...
                'See selected high res 1D PathFinders in browser', ...
                @(h,e)browseSelectedHiRes1DPF(this, [], true))
            ToolBarMethods.addButton(tb_,'pseudoBarLo.png', ...
                'See selected low res 1D PathFinders in browser', ...
                @(h,e)browseSelectedLoRes1DPF(this, [], true))
            ToolBarMethods.addButton(tb_,'comicStrip.png', ...
                'See gating sequences for selected subsets ', ...
                @(h,e)seeComicStrip(this))
            btnPin=ToolBarMethods.addButton(tb_,'pinFlashlight.png', ...
                'Highlight selected subsets ', ...
                @(h,e)highlightSubsetsNonModal(this));
            ToolBarMethods.addButton(tb_, 'windows.png',...
                Figures.TIP, @(h,e)doFigures(this, h));

            Gui.SetIcon1(btnPin, 'pinFlashlight.png',  ...
                [1 1 1], [.95 .95 .95]);
            pbeBtn=ToolBarMethods.addButton(tb_, 'pbe.png',...
                'See all events with probability bin embedding', ...
                @(h,e)seePbe(this));
            pbeBtn.setIconTextGap(2);
            pbeBtn.setText(['<html>' this.app.smallStart ...
                '<font color="blue"><u>See events</u></font>' ...
                this.app.smallEnd '</html>']);
            pbeBtn.setCursor(java.awt.Cursor.getPredefinedCursor(...
                java.awt.Cursor.HAND_CURSOR));
            this.gt.otherFigs{end+1}=this.fig;
            [this.btnGenie, this.jBtnGenie]=...
                Gui.ImageLabel([], 'genieSearch.png', ...
                'Click for help', @(h,e)genieSearch(this), ...
                this.fig, 4, 3);
            Gui.SetTransparent(this.jBtnGenie);
            this.findBranch;
            Gui.UpdateFontName(this.fig, QfTree.PROP_FONT);
            Gui.UpdateFontSize(this.fig, QfTree.PROP_FONT);
            fcs=gt.getFirstFcs(gid);
            this.fcsIdxs=fcs.findFcsIdxs(this.qf.columnNames);
            this.fcsIdxsStr=MatBasics.toString(this.fcsIdxs,'-');
            this.folder=gt.cytoGateFolder;    
            if ~this.props.has(this.wndProp)
                Gui.ShrinkNormalized(this.ax, .06, .1);
                if this.nNodes>20
                    xPerc=(this.nNodes-20)*1;
                    yPerc=(this.nNodes-20)*3;
                    if xPerc>35
                        xPerc=35;
                    end
                    if yPerc>45
                        yPerc=45;
                    end
                    Gui.Grow(this.fig, xPerc/100, yPerc/100);
                end
            end
            gt.phenogramMap.set(gid, this.fig);
            if ~isstruct(this.qf)
                this.qf.timing=toc(this.qf.timing);
            end
            try
                txt=String.MinutesSeconds(this.qf.timing);
                if ~isempty(txt)
                    P=Gui.GetPosition(gcf, 'characters');
                    nTxt=length(txt);
                    lbl=uicontrol('style', 'text', 'parent', ...
                        this.fig, 'String', txt, 'fontSize', 7, ...
                        'units', 'characters',...                    
                        'ToolTipString', ['The original matching ' ...
                        'computation cost ' txt]);
                    set(lbl, 'position', [P(3)-nTxt 0 nTxt 1.15]);
                end
            catch ex
                disp(ex);
            end
            this.figures=Figures.NewWindowMenu(this.fig, gt);
        end
        
        function doFigures(this, h)
            this.figures.doJavaMenu(h, this.figs1D.values);
        end
        
        function seePbe(this)
            this.gt.gtt.openPlotEditor(this.gt.tp, this.gid, ...
                num2str(this.fcsIdxs));
        end
        
        function fileName=getFileName(this, idx)
            [~,dsc, ~,~,key]=this.unpackNode(idx);
            key=strrep(key, ' OR ', '_');
            dsc=char(...
                edu.stanford.facs.swing.Basics.RemoveXml(...
                dsc));
            dsc=strrep(strtrim(dsc), ' ', '_');
            dsc=strrep(String.ToSystem(dsc), '"', '');
            fileName=[dsc '-' key];            
        end
        
        function file=getPngFile(this, id)
            file=[strrep(id, ' OR ', '_') '__' this.fcsIdxsStr ...
                '_' StatPlot.SFX ];
        end
        
        function file=getSvgFile(this, idx)
            file=[this.getFileName(idx) '__' this.fcsIdxsStr ...
                '-' StatPlot.SFX_SVG ];
        end
        
        function [path, file]=getPngPath(this, id)
            file=this.getPngFile(id);
            path=fullfile(this.folder, file);
            path=this.app.getCached(path);
        end
        
        function [img, preExists]=getPngImg(this, idOrInfix, ...
                prePendIfExists, scale, forBrowser)
            [fp, file]=this.getPngPath(idOrInfix);
            if ~exist(fp, 'file')
                preExists=false;
                img=Html.Img('blank.png');
            else
                preExists=true;
                if this.app.highDef
                    scale=scale*this.app.toolBarFactor;
                end
                img=Html.ImgXy(file, fileparts(fp), scale, forBrowser);
                if nargin>2 && ~isempty(prePendIfExists)
                    img=[prePendIfExists img ];
                end
            end
        end
        
        function save(this, filePath, gt_, gid_)
            QF=QfHiDM.TreeData(this);
            if nargin>2
                QF.ids=gt_.tp.getEppOrNonEppLeaves(gid_);
            end
            QF.timing=this.qf.timing;
            save(filePath, 'QF');
        end
        
        function svg(this)
            [ids, N]=this.getSelectedLeaves(true);            
            if ~askYesOrNo(Html.WrapHr(['Create vector graphics files for <br>' ...
                    'your ' String.Pluralize2('selected subset', N) '?']),...
                    'Vector graphics ...', 'center', true, 'phenogramSvg')
                return;
            end
            svgFolder=fullfile(this.gt.rootFolder, 'vectorGraphics');
            fldr=this.gt.tp.get(FcsTreeGater.PROP_FOLDER, svgFolder);
            if ~ischar(fldr)
                fldr=svgFolder;
            end
            File.mkDir(fldr);
            fldr=uigetdir(fldr, 'Pick folder for SVG output');
            if isempty(fldr) || isnumeric(fldr)
                return;
            end
            this.gt.tp.set(FcsTreeGater.PROP_FOLDER, fldr); 
            figure(this.fig);
            pu=PopUp;
            pu.setCancel(true);
            pu.initProgress(N);
            pu.setText(['Creating svg files for ' ...
                String.Pluralize2('1D PathFinder image', N)]);
            for i=1:N
                idx=this.getIdx(ids{i});
                [~,fig1D]=this.vectorGraphics(...
                    idx, true, false, false);
                file=this.getSvgFile(idx);
                file=fullfile(fldr, file);
                saveas(fig1D, file);
                svg_fix_viewbox(file, 100, file);
                pu.incrementProgress;
                if pu.cancelled
                    break;
                end
            end
            pu.close;
            if ismac
                system(['open ' String.ToSystem(fldr)]);
            else
                system(['explorer ' String.ToSystem(fldr)]);
            end
        end
        
        function motion(this)
            this.toolTip(this.findNearestNode( ...
                get(this.ax,'CurrentPoint'), true), 0, 0);
        end
        
        function [tip, bp]=toolTip(this, idx, x, y)
            tip='';
            bp=[];
            if idx==0
                return;
            end
            if ~isempty(this.mnuOpts) && this.mnuOpts.isVisible
                return;
            end
            if this.tipIsNew
                this.lastNear=0;
                return;
            end
            
            [~, dsc, freq, ttl]=this.unpackNode(idx, true);
            ids=this.nodes{idx};
            N=length(ids);
            if N==0
            elseif N==1
                key=ids{1};
            else
                key=BooleanGate.ToInfix(ids);
            end
            this.lastTipFnc=@()showTip;
            showTip;
            this.tipIsNew=true;
            MatBasics.DoLater(@(h,e)staleTip(), .71);
            function staleTip
                this.tipIsNew=false;
            end
            function showTip
                dtl=Gui.TipDetail(this.app);
                if dtl<2 && ~isnan(x)
                    tip='<html></html>';
                    bp=Gui.AddTipImgCheckBox(this.app,[],this.lastTipFnc);
                elseif dtl==2
                    tip='<html></html>';
                    if ispc
                        tip=['<html>' dsc ' ' freq '</html>'];
                        bp=Gui.AddTipImgCheckBox(this.app, [], ...
                            this.lastTipFnc);                        
                    else
                        bp=Gui.AddTipImgCheckBox(this.app, [], ...
                            this.lastTipFnc,  'North', 'South', ...
                            [dsc ' ' freq]);
                    end
                else
                    if this.app.highDef
                        scale=.3*this.app.toolBarFactor;
                    else
                        scale=.33;
                    end
                    [img, preExists]=this.getPngImg(key, '<hr>', scale,...
                        false);
                    if dtl==4
                        useImg=preExists;
                    else
                        useImg=false;
                    end
                    if N>1
                        key2=strrep(key, ' OR ', ' <b><i>or</i></b> ');
                        if N>4
                            if useImg
                                width='250';
                            else
                                width='200';
                            end
                            idTxt=['<table width="' width ...
                                'px"><tr><td><b>IDs</b>=' key2 ...
                                '</td></tr></table>'];
                        else
                            idTxt=['<b>IDs</b>=' key2];
                        end
                    else
                        idTxt=['<b>ID</b>=' key];
                    end
                    if useImg
                        tip=['<html><table cellspacing="11"><tr>'...
                            '<td align="center">' dsc ' ' freq ...
                            '<br>' idTxt '</td></tr><tr><td>' img...
                            '</td></tr></table></html>'];
                    else
                        tip=['<html><center>'  ttl '<br>' idTxt...
                            '<hr></center></html>'];
                    end
                    this.gt.app.closeToolTip;
                    if dtl==4
                        bp=this.getTipPanel(idx, ~preExists);
                    else
                        bp=this.getTipPanel(idx, false);
                    end
                end
                if ~isnan(x) && ~isnan(y)
                    if x==0 && y==0
                        if ispc
                            if this.app.highDef
                                xOff=-96;
                                yOff=1;
                            else
                                xOff=-56;
                                yOff=21;
                            end
                        else
                            xOff=-54;
                            yOff=40;
                        end
                        Gui.ShowToolTipHere(tip, this.jBtnGenie,...
                            this.ax, this.fig, this.gt.app, 5, bp, ...
                            dtl<3, xOff, yOff, this.xyIn(idx,:));
                    else
                        Gui.ShowToolTipHere(tip, this.jBtnGenie,...
                            this.ax, this.fig, this.gt.app, 5, bp, ...
                            dtl<3, x, y);
                    end
                end
            end
        end
        
        function [dsc, freq]=describe(this, idx, doHtml)
            if nargin<3
                doHtml=true;
            end
            ids=this.nodes{idx};
            N=length(ids);
            gid_=this.gid;
            gt_=this.gt;
            gtp=gt_.tp;
            if N==1
                dsc=this.nodeNames{idx};
                if doHtml
                    dsc=strrep(dsc, '^{', this.sup1);
                    dsc=strrep(dsc, '}', this.sup2);
                else
                    dsc=strrep(dsc,'^{','  (');
                    dsc=strrep(dsc, '}', ')');
                end
            elseif N>0
                dsc=sprintf('%d merged gates', N);
            else
                dsc='';
            end            
            cnt=this.qf.nodeSzs(idx);
            totalCnt=this.qf.treeSz;
            freq=[String.encodePercent(cnt, totalCnt, 1) ', '...
                gtp.encodeCount(cnt) ...
                ' events'];
            if doHtml
                if gtp.isEpp(gid_)
                    dsc=['<font color="blue">'  dsc this.eppImg '</font> '];
                else
                    dsc=['<font color="blue">'  dsc '</font> ' ];
                end
                freq=[this.sm1 ' <b>'  freq '</b>' this.sm2];                
            end
        end
        
        function optionsMenu(this, idx)
            this.app.closeToolTip;
            if isempty(this.mnuOpts)
                this.mnuOpts=PopUp.Menu;
            else
                this.mnuOpts.removeAll;
            end
            if nargin<2
                idx=0;
            end
            if idx>0 && ~this.selected.contains(idx)
                return;
            end
            cntAllLeaf=length(this.leafNames);
            [ids, cntLeaf]=this.getSelectedLeaves(true);
            cntBranch=this.getBranchCount;
            branchTxt=String.Pluralize2('branch', cntBranch, 'branches');
                
            if this.selected.size==0
                ttl=sprintf('%d subset(s) exist', cntAllLeaf);
                leafTxt=['all ' num2str(cntLeaf) ' subsets'];
            else
                leafTxt=String.Pluralize2('subset', cntLeaf);
                if cntBranch>0
                    ttl=sprintf('<b>Selections</b>: %s &amp; %s', ...
                        String.Pluralize2('subset', cntLeaf), ...
                        String.Pluralize2('branch', cntBranch, 'branches'));
                else
                    ttl=['<b>Selections</b>: ' ...
                        String.Pluralize2(  'subset', cntLeaf)];
                end
            end
            Gui.NewMenuItem(this.mnuOpts, Html.Wrap(['&nbsp;&nbsp;'...
                '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'...
                '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'...
                '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'...
                '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'...
                '&nbsp;&nbsp;&nbsp;<i>Close options</i>&nbsp;' ...
                Html.Img('rightArrow.png') '&nbsp;' ...
                Html.Img('close.gif')]), [], ...
                'blank.png');
            ml=Gui.NewMenuLabel(this.mnuOpts, ['<html>&nbsp;&nbsp;' ...
                ttl '</html>']);
            ml.setEnabled(true);
            this.mnuOpts.addSeparator;
            this.getPinComponent(idx, '', '', this.mnuOpts, ids, cntLeaf);
            mnu1D=Gui.NewMenu(this.mnuOpts, '1D PathFinders (1DPF)', ...
                [], 'pseudoBarHi.png');
            Gui.NewMenuItem(mnu1D, ['High res (' leafTxt ') '], ...
                @(h,e)browseSelectedHiRes1DPF(...
                this, [], true), 'pseudoBarHi.png');
            mi=Gui.NewMenuItem(mnu1D, ['High res (' ...
                branchTxt ' via boolean gates)'], ...
                @(h,e)browseSelectedHiRes1DPF(this, [], false), ...
                'pseudoBarHi.png');
            mi.setEnabled(cntBranch>0);
            mnu1D.addSeparator;
            Gui.NewMenuItem(mnu1D, ['Low res (' leafTxt ')'], ...
                @(h,e)browseSelectedLoRes1DPF(this, [], true), 'pseudoBarLo.png');
            mi=Gui.NewMenuItem(mnu1D, ['Low res (' ...
                branchTxt ' via boolean gates)'], ...
                @(h,e)browseSelectedLoRes1DPF(this, [], false), 'pseudoBarLo.png');
            mi.setEnabled(cntBranch>0);
            Gui.NewMenuItem(mnu1D, ['Heatmap + low res (' ...
                leafTxt ')'], @(h,e)browseAll(this, false), 'heatmap.png');
            gateMenuItems(ids);
            Gui.NewMenuItem(this.mnuOpts, [String.Pluralize2(...
                'vector graphics file', cntLeaf)],...
                    @(h,e)svg(this), 'svg.png');
                
            ids=this.getPhenogramId(ids, idx);
            if ~isempty(ids)
                Gui.NewMenuLabel(this.mnuOpts, ['<html>Boolean gate in GatingTree</html>']);
                gateMenuItems(ids);
            end
            if idx<1
                this.mnuOpts.addSeparator;
                Gui.TipMenu2(this.app, this.mnuOpts, 'blank.png');
                Gui.NewMenuItem(this.mnuOpts, ...
                    'See Heat map & 1DPF phenogram table', ...
                    @(h,e)browseAll(this, true), 'phenogram.png');
                Gui.NewMenuItem(this.mnuOpts, ...
                    'Create boolean gates for branches', ...
                    @(h,e)doBranches(this), 'tree.png',  [], true, ...
                    'Represent phenogram branches as boolean gates in GatingTree');
                Gui.NewMenuItem(this.mnuOpts, ...
                    'See median marker expression heat map', ...
                    @(h,e)showMarkerHeatMap(this, false), 'heatMap.png');
                Gui.FontMenu(this.mnuOpts, this.fig, QfTree.PROP_FONT, false);
                this.mnuOpts.addSeparator;            
                for on=1:-1:0
                    if on==1
                        mnu=Gui.NewMenu(this.mnuOpts,'Select', [], 'blank.png');
                    else
                        mnu=Gui.NewMenu(this.mnuOpts,'Deselect', [], 'blank.png');
                    end
                    doSelectionMenu(mnu, on, 0);
                    doSelectionMenu(mnu, on, 1);
                    doSelectionMenu(mnu, on, -1);
                end
                this.mnuOpts.addSeparator;
                if idx>0
                else
                    this.figures.doJavaMenu(this.mnuOpts, this.figs1D.values);
                end
                Gui.TipMenu2(this.app, this.mnuOpts, 'blank.png');
            end
            
            function doSelectionMenu(mnu, on, which)
                if which==0 %all
                    scope='All subsets';
                elseif which==1 %leaves
                    scope='Leaf subsets';
                elseif which==-1 %branches
                    scope='Phenogram branches';
                end
                Gui.NewMenuItem(mnu, ...
                    scope, @(h,e)doSelection(on, which), 'blank.png');
            end
            function doSelection(on, which)
                if which~=0
                    idxs=this.sortXyAndAx('Phenogram order', which);
                else
                    idxs=this.sortXyAndAx('north west', which);
                end
                nIdxs=length(idxs);
                this.selected.clear;
                if on
                    for i=1:nIdxs
                        this.selected.add( idxs(i) );
                    end
                else
                    for i=1:nIdxs
                        this.selected.remove( idxs(i) );
                    end
                end
                this.showSelected;
            end
            
            function gateMenuItems(ids)
                N=length(ids);
                Gui.NewMenuItem(this.mnuOpts, label(N, 'Gating sequence'), ...
                    @(h,e)seeComicStrip(this, ids, 0), ...
                    'comicStrip.png');                
                lblTree=Gui.NewMenuItem(this.mnuOpts, label(N, 'GatingTree position'), ...
                    @(h,e)treePosition(this, ids, 0), 'tree.png');
                if N==1
                    pid=this.gt.tp.getParentNode(ids{1});
                    lvlParent=this.gt.tp.getColorLevel(pid, false);
                    lvlChild=this.gt.tp.getColorLevel(ids{1}, true);
                    this.gt.app.treeIcon.setIconWith2ndColor(lblTree, ...
                        'tree.png', lvlChild, [0 0 1], lvlParent);                
                end
                Gui.NewMenuItem(this.mnuOpts, label(N, 'PlotEditor'), ...
                    @(h,e)plotEditor(this, ids, 0), 'eye.gif'); 
            end
            
            function item=label(N, item)
                if N>1
                    item=String.Pluralize2(item, N);
                end
            end
        end
        
        function [ids, cbBrowse, idxs, numCols, cancelled, totalIds]=chooseIds(...
                this, ids, title, suffix, supportBrowsing, chosen, allChosen)
            N=length(ids);
            totalIds=ids;
            if nargin<7
                allChosen=0:this.qf.numLeaves-1;
                if nargin<6
                    chosen=0:N-1;
                    if nargin<5
                        supportBrowsing=true;
                    end
                end
            end
            numCols=0;
            cb=Gui.CheckBox(['<html>Choose from <b><u>ALL</u></b> '...
                num2str(this.qf.numLeaves) ' subsets</html>'], ...
                false, [],[],@close);
            txtAll=Html.WrapHr(['<h3>Select one or more of the '...
                String.Pluralize2('subset', N) ' below <br>to see their '...
                suffix '...']);
            if N<this.qf.numLeaves
                txtMsg=Html.WrapHr(['<h3>Choose one or '...
                    ' more of the ' String.Pluralize2(['<u><font color'...
                    '="#222222">currently selected</font></u> subset'], N) ...
                    '<br>below to see their ' suffix '...']);
                msg_=Gui.BorderPanel;
                msg_.add(Gui.Label(txtMsg), 'North');
                msg_.add(cb, 'East');
            else    
                msg_=txtAll;
            end
            options=cell(1, N);
            idxs=zeros(1, N);
            for i=1:N
                [str, idxs(i)]=this.unpackId(ids{i});
                options{i}=['<html>' str '</html>'];
            end
            if supportBrowsing
                [south, cbBrowse, jt, prop]=...
                    GatingTree.BrowseSeqGui(this.app);
                [choices, cancelled]=mnuMultiDlg(msg_, ...
                    title, options, chosen, false, ...
                    true, south, 'south west buttons');
            else
                cbBrowse=[];
                [choices, cancelled]=mnuMultiDlg(msg_, ...
                    'Choose subsets...', options, chosen, false, true);
            end
            if cb.isSelected
                ids=this.getSelectedLeaves(false);
                totalIds=ids;
                [ids, cbBrowse, idxs, numCols, cancelled]=chooseIds(...
                    this, ids, title, suffix, supportBrowsing, allChosen);
                return;
            end
            if cancelled
                ids={};
            else
                ids=ids(choices);
                idxs=idxs(choices);
            end
            if supportBrowsing
                [ok, num]=Gui.SetNumberField(jt, this.app, prop, 0, 20,...
                    'Gates per table row', false, ...
                    '<br>(or 0 for unlimited gates per row)');
                if ok
                    numCols=num;
                else
                    ids={};
                end
            end
            
            function close(h,e)
                wnd=Gui.WindowAncestor(h);
                wnd.dispose;
            end
        end
        
        function [ttl, idx, sym, dsc, freq, key, N]=unpackId(this, id, doHtml)
            if nargin<3
                doHtml=true;
            end
            idx=this.getIdx(id);
            [sym, dsc, freq, ttl, key, N]=this.unpackNode(idx, doHtml);
        end
        
        function [sym, dsc, freq, ttl, key, N]=unpackNode(this, idx, doHtml)
            if nargin<3
                doHtml=false;
            end
            isLeaf=idx<=length(this.leafNames);
            if isLeaf
                sym=[this.getSubsetSymbol(idx, true) '&nbsp;&nbsp;'];
            else
                sym='';
            end
            [dsc, freq]=this.describe(idx, doHtml);
            ttl=[sym strtrim(dsc) ', ' freq ]; 
            if nargout>=5
                [key,~,N]=this.getIdOrInfix(idx);
            end
        end
        
        function contextMenu(this, idx, H)
            this.optionsMenu(idx);
            if nargin<3
                if idx==0
                    [X, Y]=Gui.GetCurrentXY(this.ax, this.fig,this.jBtnGenie);
                    this.mnuOpts.show(this.jBtnGenie, X, Y);
                else
                    xx=get(this.plotHs(idx), 'xdata');
                    yy=get(this.plotHs(idx), 'ydata');
                    [X, Y]=Gui.GetCurrentXY(this.ax, this.fig, ...
                        this.jBtnGenie, [xx yy]);
                    if idx<=this.nLeafs
                        this.mnuOpts.show(this.jBtnGenie, X-45, Y);
                    else
                        if this.app.highDef
                            xOff=335*this.app.toolBarFactor;
                        else
                            xOff=415;
                        end
                        this.mnuOpts.show(this.jBtnGenie, X-xOff, Y);
                        %this.mnuOpts.show(this.btnSvg, -235, -85);
                    end
                end
            else
                this.mnuOpts.show(H, 5, 15);
            end
        end
        
        function genieSearch(this)
            this.optionsMenu;
            mc=MenuAutoComplete(this.mnuOpts);
            toc;
            mc.popUpDlg([], this.fig, 'south west', ...
                'Search Phenogram functions...', ...
                'Find ANY operation', 9, 50, 40);
        end
        
        function wdcu(this, cp)
            idx=this.findNearestNode( cp, false);
            if idx>0
                this.highlightSubsetsNonModal;
            end
        end
        
        function wbu(this, cp)
            idx=this.findNearestNode( cp, false);
            isShift=edu.stanford.facs.swing.KeyState.IsShiftPressed;
            isMulti=this.doubleClicker.multiSelectWhenButtonDown;
            a=get(gcbf, 'SelectionType');
            modifiers = get(gcf,'currentModifier');
            disp('THE MODIFIERS & selection type ARE');
            disp(modifiers);
            disp(a);
            disp('ALL SHOWN');
            if ~isShift && ~isMulti && strcmp(a, 'alt') 
                this.contextMenu(0);
            elseif idx>0
                if ispc
                    multiKey='Ctrl';
                else
                    multiKey='command';
                end
                this.gt.app.closeToolTip;
                priorSelected=this.selected.size;
                contains=this.selected.contains(idx);
                selectTip=[' <font color="red">without</font> shift or '...
                    multiKey ' pressed'];
                if isShift && this.selected.size>0
                    selectTip=' with <font color="red">shift</font> pressed';
                    if contains
                        this.selected.remove(idx);
                    else
                        a=this.selected.toArray;
                        lIdx=a(end);
                        D1=pdist2(this.xyNormalized(lIdx,:), this.xyNormalized(idx,:));
                        for j=1:length(this.nodes)
                            if ~this.selected.contains(j)
                                D2=pdist2(this.xyNormalized(j,:), this.xyNormalized(idx,:));
                                D3=pdist2(this.xyNormalized(j,:), this.xyNormalized(lIdx,:));
                                if D1>D2 && D1>D3
                                    this.selected.add(j);
                                end
                            end
                        end
                        this.selected.add(idx);                        
                    end
                elseif isMulti
                    selectTip=[' with <font color="red">' ...
                        multiKey '</font></i> pressed'];
                    if contains
                        this.selected.remove(idx);
                    else
                        this.selected.add(idx);
                    end
                else
                    this.selected.clear;
                    if ~contains
                        this.selected.add(idx);
                    end
                end
                nowSelected=this.selected.size;
                if priorSelected>nowSelected
                    word='Deselection';
                else
                    word='Selection';
                end
                this.showSelected;
                this.contextMenu(idx);
                yOff=-12;
                this.app.closeToolTip;
                Gui.ShowToolTipHere(['<html>' this.app.smallStart ...
                    word selectTip  this.app.smallEnd '</html>'], ...
                    this.jBtnGenie, this.ax, this.fig, this.gt.app, ...
                    3, [], true, -110, yOff);
            end
        end
        
        function browseSelectedHiRes1DPF(this, idx, useLeafs)
            if nargin<3
                idxs=idx;
                N=1;
            else
                if useLeafs
                    [~, N,idxs]=this.getSelectedLeaves(true);
                else
                    [idxs, N]=this.getSelectedBranches;
                end
            end
            figure(this.fig);
            figs={};
            pu=PopUp;
            pu.setCancel(true);
            pu.initProgress(N);
            pu.setText(['Creating ' String.Pluralize2('1D PathFinder image', N)]);
            for i=1:N
                idx=idxs(i);
                [~,fig1D]=this.vectorGraphics(idx, true, false);
                figure(fig1D);
                pu.incrementProgress;
                figs{end+1}=fig1D;
                if pu.cancelled
                    break;
                end
            end
            if N>1                
                Gui.CascadeFigs(figs, false, true, 15);
            end
            pu.close;
        end
        
        function browseSelectedLoRes1DPF(this, idx, useLeafs)
            if nargin<3
                idxs=idx;
                N=1;
            else
                if useLeafs
                    [~, N,idxs]=this.getSelectedLeaves(true);
                else
                    [idxs, N]=this.getSelectedBranches;
                end
            end
            
            this.createImg(idxs, false);
            html='<html>';
            for i=1:N
                idx=idxs(i);
                key=this.getIdOrInfix(idx);
                img=this.getPngImg(key, [], .66, true);
                [~, ~, ~, ttl]=this.unpackNode(idx, true);
                html=[html ttl '<br>' img '<hr>'];
            end
            html=[html '</html>'];
            Html.Browse(html);
        end
        
        function showSelected(this)
            N=length(this.plotHs);
            if this.selected.size==0
                for i=1:N
                    H=this.plotHs(i);
                    set(H, 'markerEdgeColor', this.edgeColors{i});
                    set(H, 'lineWidth', this.lineWidths(i));
                end                
            else
                hs=java.util.HashSet;
                it=this.selected.iterator;
                while it.hasNext
                    idx=it.next;
                    if idx==0
                        continue;
                    end
                    ids=this.nodes{idx};
                    nIds=length(ids);
                    if nIds>1
                        for i=1:nIds
                            hs.add(ids{i});
                        end
                    end
                end
                for i=1:N
                    H=this.plotHs(i);
                    if i>this.nLeafs
                        selColor=[.25 .25 .25];
                        selSize=4;
                    else
                        selColor=[0 0 0];
                        selSize=4;
                    end
                    if this.selected.contains(i)
                        set(H, 'markerEdgeColor', selColor, 'lineWidth', selSize);
                    else
                        if i>this.nLeafs
                            selColor=[.35 .35 .85];
                            selSize=2;
                        end
                        ids=this.nodes{i};
                        nIds=length(ids);
                        highlighted=0;
                        for j=1:nIds
                            if hs.contains(ids{j})
                                highlighted=highlighted+1;
                            else
                                break;
                            end
                        end
                        if highlighted<nIds
                            set(H, 'markerEdgeColor', this.edgeColors{i});
                            set(H, 'lineWidth', this.lineWidths(i));
                        else
                            set(H, 'markerEdgeColor', selColor, ...
                                'lineWidth', selSize);
                        end
                    end
                end
            end
        end
        
        function str=getSubsetSymbol(this, idx, big)
            c=get(this.plotHs(idx), 'MarkerFaceColor');
            str=['<font  ' Gui.HtmlHexColor(c)...
                '>&#8903;</font>'];
            str=['<font size="35">' str '</font>'];
        end
        
        function showMarkerHeatMap(this, justBrowse, phenogramOrientation)
            try
                if ~isempty(this.jdMrker) && ~justBrowse
                    this.jdMrker.setVisible(true);
                    this.jdMrker.toFront;
                    return;
                end
                try
                    args.measurements=this.qf.measurements;
                    args.rawMeasurements=this.qf.rawMeasurements;
                catch ex
                    msg('Rebuild the phenogram first');
                    return;
                end
                args.app=this.app;
                args.names=this.originalNames;
                cns=this.qf.columnNames;
                args.measurementNames=cns;
                args.freqs=this.qf.nodeSzs/this.qf.treeSz;
                args.subsetSymbol=@symb;
                args.fnc1DPathFinder=@pathFinder;
                args.mouseEar=@mouseEar;
                if isempty(this.showMeasurements)
                    this.showMeasurements=false(1, length(args.measurementNames));
                end
                [R,C]=size(args.measurements);
                options=cell(1,C);
                for col=1:C
                    options{col}=['<html><font color="blue"><i>' ...
                        num2str(col) '.</i></font> ' ...
                        args.measurementNames{col}  '</html>'];
                end
                this.heatMapArgs=args;
                showNoMatch=true;
                heatMap=Gui.BorderPanel(10,1);
                [top, ~, cb]=MarkerHeatMap.DropDowns(...
                    @(h,e)setColorScheme, 'phenogramHeatMapOrder', ...
                    @(h,e)setMap, 1, ...
                    {'Phenogram order', 'Reverse phenogram order'}, ...
                    'South', 'North', ' low to high expression');
                heatMap.add(top, 'North');
                jd=[];
                tp1=[];
                maxN=0;
                I=[];
                setMap;
                if justBrowse
                    browse(phenogramOrientation);
                    return;
                end
                nItems=12;
                if maxN<7
                    nItems=9;
                elseif maxN<12
                elseif maxN<20
                    nItems=20;
                else
                    nItems=25;
                end
                figure(this.fig);
                browseBtn1=Gui.NewBtn('Browse', @(h,e)browse(false),...
                    'Click to see web page of 1D PathFinders', ...
                    'table.gif');
                browseBtn2=Gui.NewBtn('Browse', @(h,e)browse(true),...
                    'Click to see web page of 1D PathFinders', ...
                    'phenogram.png');
                
                sw=Gui.Panel;
                sw.add(browseBtn1);
                sw.add(browseBtn2);
                sw=Gui.AddTipImgCheckBox(this.app, sw, this.lastTipFnc, ...
                    'East', 'West');
                [~,~,jd]=mnuMultiDlg(struct(...
                    'allMsg', ['<html><b>All</b><br>' this.sm1 ...
                    '(Selections show in<br> phenogram labels)' this.sm2 ...
                    '<hr></html>'],...
                    'noCancel', true, ...
                    'where', 'west++',...
                    'icon', 'none', 'modal', false, ...
                    'msg', ['<html><h2>' num2str(C) ' median marker '...
                    'expressions for ' num2str(R) ...
                    ' subsets</h2></html>'], 'checkBoxFnc',...
                    @(h, e, idxs, checkBoxes)markerSelected(idxs)),...
                    'Median marker expression heat map', options, ...
                    find(this.showMeasurements)-1, false, true, sw, ...
                    'south west buttons', heatMap,'West',  nItems);
                this.jdMrker=jd;
            catch ex
                BasicMap.Global.reportProblem(ex);
            end
            
            function markerSelected(idxs)
                this.updateMrkrExpr(idxs)
                args.selectedMarkers=idxs;
                this.heatMapArgs=args;
                setMap;
            end
            
            function setColorScheme
                setMap;
            end
            
            function [html, nCols]=topHtml
                mns=args.measurementNames;
                nCols=length(mns);
                if QfTree.SHOW_MARKER_HTML
                    markerTable=['<table><tr><td><font color="blue">' ...
                        '<i>Markers</i></font><hr></td></tr>'];
                    for i=1:nCols
                        markerTable=[markerTable '<tr><td>' ...
                            num2str(i) '. ' mns{i} '</td></tr>'];
                    end
                    markerTable=[markerTable '</table>'];
                else
                    markerTable='';
                end
                nCols=length(I);
                this.saveMainPng;
                if ispc
                    scale=.7;
                else
                    scale=.46;
                end
                img=this.getMainPngImg(scale, true);
                header='<center><h1>';
                ax_=get(this.fig, 'currentaxes');
                txt=get(get(ax_, 'title'), 'String');
                if iscell(txt)
                    N=length(txt);
                    for ii=1:N
                        header=[header txt{ii}];
                        if ii<N
                            header=[header '<br>'];
                        end
                    end
                else
                    header=[header txt];
                end
                header=strrep(header, '^{', '<sup>');
                header=strrep(header, '}', '</sup>');
                header=[header '</h1></center><hr>'];
                htmlHeat=MarkerHeatMap.Html(args, this.ax, ...
                    showNoMatch, I, [], true, true);
                
                html=['<html>' header ' <table cellspacing="8"><tr>'...
                     '<td  valign="top">' htmlHeat '</td><td>&nbsp;'...
                     '&nbsp;&nbsp;</td><td valign="top">' markerTable ...
                    '</td><td valign="top">' img '</td></tr></table>'];
            end
            
            function browse(phenogramOrientation)
                figure(this.fig);
                pu=PopUp('Creating webpage html', 'west+');
                [html, N]=topHtml;
                if phenogramOrientation
                    this.getGridHtml(true, html);
                    pu.close;
                    return;
                end
                figure(this.fig);
                [ok, cancelled]=askYesOrNo(Html.WrapHr(['Show ' String.Pluralize2(...
                    '1D PathFinder', N) ' for <i>each</i> phenogram leaf?']), ...
                    'Browsing ...', 'west+', ...
                    false, 'phenogramBrowse');
                if cancelled
                    pu.close;
                    return;
                end
                if ok
                    selectedIds=StringArray(this.getSelectedLeaves(true));
                    if selectedIds.N>0
                        idxs=[];
                        for i=1:N
                            idx=I(i);
                            id=this.nodes{idx};
                            id=id{1};
                            debug=this.leafNames{idx}
                        
                            if selectedIds.contains(id)
                                idxs(end+1)=idx;
                            end
                        end
                    else
                        idxs=I;
                    end
                    this.createImg(idxs, false);
                    N=length(idxs);
                    for i=1:N
                        idx=idxs(i);
                        debug=this.leafNames{idx}
                        id=this.nodes{idx};
                        id=id{1};
                        img=this.getPngImg(id, [], .7, true);
                        html=[html img '<hr>'];
                    end
                end
                html=[html '</html>'];
                Html.Browse(html);
                pu.close;
            end
            
            function pathFinder(idx, where)
                this.vectorGraphics(idx);
            end
            
            function [tip, bp]=mouseEar(ax, ismotion, idx, x, y)
                if ismotion==-1
                    plotHs_=this.plotHs(idx);
                    lw=get(plotHs_, 'linewidth');
                    clr=get(plotHs_, 'markerEdgeColor');
                    set(plotHs_, 'markerEdgeColor', 'red', 'lineWidth', 7);
                    for id = 1:3 % Repeat 3 times
                        set(plotHs_, 'visible', 'off');
                        pause(0.15);
                        set(plotHs_, 'visible', 'on');
                        pause(0.15);
                    end
                    set(plotHs_, 'markerEdgeColor', clr, 'lineWidth', lw);
                else
                    [tip, bp]=this.toolTip(idx, nan, nan);
                end
            end
            
            function str=symb(idx)
                str=this.getSubsetSymbol(idx);    
            end
            
            function setMap
                subsetOrder=char(cb.getSelectedItem);
                if ~isempty(tp1)
                    heatMap.remove(tp1);
                end
                I=this.sortXyAndAx(subsetOrder);
                maxN=length(I);
                [~, tp1]=MarkerHeatMap.Html(args, this.ax, ...
                    showNoMatch, I, [], true, justBrowse);
                heatMap.add(tp1, 'Center');
                if ~isempty(jd)
                    jd.pack;
                end
            end
        end
        
        function missing=countMissingPngs(this, idxs, create, which)
            if nargin<4
                which=0;
                if nargin<3
                    create=false;
                    if nargin<2
                        idxs=[];
                    end
                end
            end
            if isempty(idxs)
                if which==0
                    idxs=[1:length(this.nodes)];
                elseif which==1
                    idxs=[1:length(this.leafNames)];
                else
                    idxs=[length(this.leafNames)+1:length(this.nodes)];
                end
            end
            N=length(idxs);
            if create
                figure(this.fig);
                pu=PopUp;
                pu.setCancel(true);
                pu.initProgress(N);
                pu.setText(['Creating ' String.Pluralize2('1D PathFinder image', N)]);
                drawnow;
            end
            missing=[];
            for i=1:N
                idx=idxs(i);
                if ~this.vectorGraphics(idx, create, false, false)
                    missing(end+1)=idx;
                end
                if create
                    pu.incrementProgress;
                    if pu.cancelled
                        break;
                    end
                end
            end
            if create
                pu.close;
            end
        end
        
        function grid=getGrid(this)
            factor=10;
            done=false;
            while ~done
                done=true;
                factor=factor+1;
                N=length(this.nodes);
                uRows=unique(ceil(this.xyIn(:,2)*factor)+1);
                uCols=unique(ceil(this.xyIn(:,1)*factor)+1);
                nRows=max(uRows);
                nCols=max(uCols);
                grid=zeros(nRows, nCols);
                nL=length(this.leafNames);
                for i=1:N-1
                    if i<=nL
                        debug=this.leafNames{i};
                    else
                        debug=this.qf.branchNames{i-nL};
                    end
                    xy=ceil(this.xyIn(i,:)*factor)+1;
                    fprintf('%d %d %s\n', xy(2), xy(1), debug);
                    if grid(xy(2), xy(1))~=0
                        done=false;
                        disp('COLLISION');
                        break;
                    end
                    grid(xy(2), xy(1))=i;
                end
            end
            root=zeros(nRows, 1);
            xy=ceil(this.xyIn(N,:)*factor)+1;
            root(xy(2))=N;
            grid=[root grid];
        end
        
        function browseAll(this, phenogramOrientation)
            this.showMarkerHeatMap(true, phenogramOrientation);
        end
        
        function html=getGridHtml(this, browseToo, html)
            if nargin<3
                html='<html>';
                if nargin<2
                    browseToo=true;
                end
            end
            this.createImg([], false);
            grid=this.getGrid;
            ttl=['<h2>' String.Pluralize2('1D PathFinder', length(this.nodes))...
                ' arranged by phenogram''s associations</h2>(<b><i>Scroll vertically ' ...
                '& horizontally to see all)<hr>'];
            html=[html ttl '<table border="0"  cellpadding="0" cellspacing="0">'];
            startHtml='<td><table border="1" cellpadding="0" cellspacing="0"><tr><td>';
            endHtml='</td></tr></table></td>';
            [R,C]=size(grid);
            
            for r=1:R
                if any(grid(r,:))
                    html=[html '<tr>'];
                    for c=1:C
                        if any(grid(:,c))
                            idx=grid(r,c);
                            if idx==0
                                html=[html '<td></td>'];
                            else
                                [~, ~, ~, ttl, key, N]=...
                                    this.unpackNode(idx, true);
                                ttl=strrep(ttl, this.app.smallStart, '<small>');
                                ttl=strrep(ttl, this.app.smallEnd, '</small>');
                                img=this.getPngImg(key,[], .22, true);
                                if N==1
                                    html=[html startHtml ttl '<hr>'...
                                        img endHtml];
                                else
                                    html=[html startHtml ttl '<hr>' img endHtml];
                                end
                            end
                        end
                    end
                    html=[html '</tr>'];
                end
            end
            html=[html '</table></html>'];
            if browseToo
                Html.Browse(html)
            end
        end
        
        function I=sortXyAndAx(this, distFromWhere, which, rowsIdxs)
            if nargin<4
                rowsIdxs=[]; 
                if nargin <3
                    which=1;
                    if nargin<2
                        distFromWhere=[];
                        which=0;
                    end
                end
            end
            N1=length(this.leafNames);
            N0=length(this.nodes);
            if isempty(rowsIdxs)
                if which==1
                    rowIdxs=1:N1;
                elseif which==-1
                    rowIdxs=N1+1:N0;
                else
                    rowIdxs=1:N0;
                end
            end
            if isempty(distFromWhere)
                distFromWhere='north west';
            end
            xy=this.xyIn(rowIdxs,:);                
            if String.StartsWith(distFromWhere, 'Ascending')
                I=sort_(this.leafNames(rowIdxs));
            elseif String.StartsWith(distFromWhere, 'Descending')
                I=sort_(this.leafNames(rowIdxs));
                I=flip(I);
            elseif isequal('Phenogram order', distFromWhere)
                [~,I]=sort(xy(:,2));
            elseif isequal('Reverse phenogram order', distFromWhere)
                [~,I]=sort(xy(:,2));
                I=flip(I);
            elseif String.EndsWith(distFromWhere, 'expression')
                I=MarkerHeatMap.SortByMeasurement(distFromWhere, this.heatMapArgs);                
            else
                I=MarkerHeatMap.SortDists(this.ax, xy, distFromWhere);
            end
            %rowIdxs=find(rowIdxs);
            I=rowIdxs(I);
            try
                this.leafNames(I)
            catch ex
            end
            %rowIdxs=rowIdxs(I);
            function I=sort_(strs)
                out=StringArray.ToLower(strs);
                N=length(out);
                for i=1:N
                    out{i}=strrep(out{i}, '\bf', '');
                end
                [out, I]=sort(out);
            end
        end
        function findBranch(this, pid, brIdx)
            if nargin<2
                pid=this.gid;
                brIdx=length(this.nodes);
                this.gtMap=cell(1, brIdx);                
            end
            ids=this.nodes{brIdx};
            infix=BooleanGate.ToInfix(ids) ;
            fprintf('%s pid=%s\n\t\t boolean gate=%s\n', ...
                String.Pad('', brIdx*2, '>'), pid, infix);
            [ok, newPid]=this.gt.hasSpecificSampleSelections(pid, infix);
            if ok
                this.gtMap{brIdx}=newPid;
                this.gt.tp.setPhenogram(newPid)
                nextBrIdx=brIdx-this.qf.numLeaves;
                if nextBrIdx>0
                    subIdxs=this.qf.phyTree(nextBrIdx,:);
                    N=length(subIdxs);
                    for i=1:N
                        this.findBranch(newPid, subIdxs(i));
                    end
                end
            end
        end

        function doBranches(this, pid, brIdx, pu)
            if nargin<2
                figure(this.fig);
                if ~askYesOrNo(Html.WrapHr(['Represent phenogram branches '...
                        '<br>as boolean gates in GatingTree?']))
                    return;
                end
                figure(this.fig);
                pu=PopUp('Creating/confirming tree branches',...
                    'One moment...');
                brIdx=length(this.nodes);
                pu.initProgress(brIdx);
                pid=this.gid;
                Gui.setEnabled(this.fig, false);
                this.tb.setEnabled(false);
            end
            ids=this.nodes{brIdx};
            nIds=length(ids);
            if nIds==0
                msg('Rebuild phenogram, gates have changed', 4);
                pu.close;
                return;
            elseif nIds==1
                word='leaf';
            else
                word='branch';
            end
            txt=sprintf('Phenogram %s  #%d, %s', word, brIdx, ...
                String.Pluralize2('subset', length(ids)));
            pu.setText2(txt);
            pu.incrementProgress;
            infix=BooleanGate.ToInfix(ids) ;
            fprintf('%s pid=%s\n\t\t boolean gate=%s\n', ...
                String.Pad('', brIdx*2, '>'), pid, infix);
            [ok, newPid]=this.gt.hasSpecificSampleSelections(pid, infix);
            if ~ok
                if nIds==1
                    gn=this.describe(brIdx, false);
                else
                    gn=txt;
                end
                newPid=BooleanGate.Compute(this.gt, false, infix, true, ...
                    pid, gn);
                this.gt.ensureVisible(newPid, true);
            end
            nextBrIdx=brIdx-this.qf.numLeaves;
            if nextBrIdx>0
                subIdxs=this.qf.phyTree(nextBrIdx,:);
                N=length(subIdxs);
                for i=1:N
                    if i==2
                        if nargin<2
                            disp(newPid);
                        end
                    end
                    this.doBranches(newPid, subIdxs(i), pu);
                end
            end
            if nargin<2
                if ok
                    this.gt.ensureVisible(newPid, true);
                end
                this.findBranch;
                Gui.setEnabled(this.fig, true, true);
                this.tb.setEnabled(true);
                pu.close;
            end
        end
        
        function updateMrkrExpr(this, idxs)
            J=MarkerHeatMap.ColorScheme;
            nJ=size(J,1);
            if nargin>1
                this.showMeasurements(:)=false;
                this.showMeasurements(idxs)=true;
            end
            [R,C]=size(this.qf.measurements);            
            for row=1:R
                H=this.textHs(row);
                txts=get(H, 'String');
                if length(txts)>=4
                    while length(txts)>=4
                        txts(end)=[];
                    end
                end
                if any(this.showMeasurements)
                    strs={''};
                    done=0;
                    for col=1:C
                        if ~this.showMeasurements(col)
                            continue;
                        end
                        done=done+1;
                        pk=this.qf.measurements(row,col);
                        cIdx=floor(pk*nJ);
                        if cIdx<1
                            cIdx=1;
                        end
                        c=J(cIdx,:);
                        if mod(done, 7)==0
                            strs{end+1}='';
                        end
                        strs{end}=[strs{end} '\color[rgb]{' num2str(c(1)) ...
                            ',' num2str(c(2)), ',' num2str(c(3)) '}' ...
                            num2str(col) ' '];
                    end
                    txts=[txts' strs];
                end
                set(H, 'String', txts);
            end
        end
        
        function browseLowRes1DPF(this, idx, sym, dsc, freq)
            this.createImg(idx, false);
            key=this.getIdOrInfix(idx);
            img=this.getPngImg(key, [], .66, false);
            msg(['<html>' sym dsc freq '<hr>' img '</img>'], 0, 'west++', ...
                edu.stanford.facs.swing.Basics.RemoveXml(dsc), 'none');
        end
        
        function [key, ids, N]=getIdOrInfix(this, idx)
            ids=this.nodes{idx};
            N=length(ids);
            if N==1
                key=ids{1};
            else
                key=BooleanGate.ToInfix(ids);
            end
        end
        
        function [pngExists, fig1D, key]=...
                vectorGraphics(this, idx, create, showProgress, keepVisible)
            pngExists=false;
            if nargin<5
                keepVisible=true;
                if nargin<4
                    showProgress=true;
                    if nargin<3
                        create=true;
                    end
                end
            end
            pu=[];
            ids=this.nodes{idx};
            N=length(ids);
            tempGate=false;
            if N==1
                id=ids{1};
                if showProgress
                    figure(this.fig);
                    pu=PopUp(['<html><center>Building '...
                        ' 1D  PathFinder for <br>'...
                        this.gt.describe(id) ...
                        '<hr></center></html>']);
                end
                key=id;
            else
                infix=BooleanGate.ToInfix(ids);
                if ~this.figs1D.contains(infix)
                    if showProgress
                        figure(this.fig);
                        pu=PopUp(['<html><center>Building '...
                            ' 1D  PathFinder for '...
                            String.Pluralize2('gate', N) ...
                            '<br><b>' this.app.smallStart ...
                            infix this.app.smallEnd ...
                            '</b><hr></center></html>']);
                    end
                    id=this.gtMap{idx};
                    if (isempty(id) || ~this.gt.tp.nodeExists(id)) && create
                        tempGate=true;
                        id=BooleanGate.Compute(this.gt, false, ...
                            infix, true);
                    end
                end
                key=infix;
            end
            if ~create
                [fp, file]=this.getPngPath(key);
                pngExists=exist(fp, 'file');
            else
                fig1D=this.figs1D.get(key);
                if isempty(fig1D)
                    [saved, fig1D]=this.create1DPF(id, key, idx, keepVisible);
                    if ~saved                       
                        pngExists=true;
                    end
                elseif showProgress
                    figure(fig1D);
                end
            end
            if ~isempty(pu)
                pu.close;
            end
            if tempGate %more than one gate ... BOOLEAN gate
                this.gt.removeGate(id, true, false);
            end
        end
        
        function [cancelled, fig1D]=create1DPF(this, id, key, idx, keepVisible)  
            if nargin<5
                keepVisible=true;
            end
            cancelled=false;
            oneGate=strcmp(id,key);
            if isempty(this.fg)
                this.fg=this.gt.goToNode(id);
            else
                this.fg.goToNode(id);
            end
            if isstruct(this.qf)
                if isfield(this.qf, 'columnNames')
                    colNames=this.qf.columnNames;
                else
                    colNames=[];
                end
            else
                colNames=this.qf.columnNames;
            end
            prior=this.gt.multiProps.getAll(StatPlot.PROP12);
            this.gt.multiProps.set(StatPlot.PROP12, '3');
            if isempty(colNames)
                cancelled=false;
                fig1D=StatPlot.Go(this.fg);
            else
                info=QfTreeInfo(colNames, true, false);
                fig1D=StatPlot.Go(this.fg, info);
            end
            Gui.ImageLabel(Html.Wrap(['&nbsp;&nbsp;' this.getSubsetSymbol(idx)]), [],[],[], fig1D)
            this.gt.multiProps.setAll(StatPlot.PROP12, prior);
            if keepVisible
                this.figs1D.set(key, fig1D);
            end
            sz=this.qf.nodeSzs(idx);
            treeSz=this.qf.treeSz;
            ttl=this.gt.describeForTex(id, sz, treeSz);
            if String.Contains(ttl{1}, ' | ')
                ttl{1}=['Gate IDs: ' strrep(ttl{1}, ' | ', ' or ')];
            end
            if idx<=length(this.leafNames)
                who=this.gt.who(id);
                String.IndexOf(ttl{1}, who)
                ttl{1}=strrep(ttl{1}, who, ['\bf' this.leafNames{idx}]);
            end
            ax_=get(fig1D, 'currentaxes');
            T=title(ax_, ttl);
            if ~isempty(colNames)
                nD=length(colNames);
                if nD<7
                    factor=.9;
                else
                    factor=.94;
                end
                P=get(T, 'position');
                set(T, 'Position', [P(1) P(2)*factor P(3)])
            end
            Gui.UpdateFontName(fig1D, QfTree.PROP_FONT)
            Gui.UpdateAndAdjustFontSize(fig1D, QfTree.PROP_FONT, ...
                BasicMap.Global,-3)
            this.gt.otherFigs{end+1}=fig1D;
            pid=this.gt.tp.getParentNode(id);
            lvlParent=this.gt.tp.getColorLevel(pid, false);
            lvlChild=this.gt.tp.getColorLevel(id, true);
            if keepVisible
                Gui.SetToRight(fig1D, this.fig, false, 10);
                Gui.FitFigToScreen(fig1D);
                Gui.RepositionOnSameScreenIfRequired(fig1D);
                Gui.SetFigVisible(fig1D);
                set(fig1D, 'CloseRequestFcn', @(h, e)hush2(key, h), 'visible', 'on');
                drawnow;
                tb_=ToolBar(findall(fig1D,'tag','FigureToolBar'));
                Gui.Locate(fig1D, this.fig, 'west++')
                
                if oneGate
                    tbTree=ToolBarMethods.addButton(tb_, 'tree.png', ...
                        'See position in gating tree', ...
                        @(h,e)treePosition(this, id, 0));
                    this.gt.app.treeIcon.setIconWith2ndColor(tbTree, ...
                        'tree.png', lvlChild, [0 0 1], lvlParent);
                    
                    ToolBarMethods.addButton(tb_, 'eye.gif', ...
                        'Open PlotEditor', ...
                        @(h,e)plotEditor(this, id, 0));
                    ToolBarMethods.addButton(tb_, 'comicStrip.png', ...
                        'See full gating sequence', ...
                        @(h,e)seeComicStrip(this, id, 0));
                    P=get(fig1D, 'OuterPosition');
                    if ~isdeployed
                        %ensure all buttons are visible
                        set(fig1D, 'OuterPosition', [P(1) P(2) P(3)*1.05 P(4)*1.1]);
                    else
                        set(fig1D, 'OuterPosition', [P(1) P(2) P(3) P(4)*1.1]);
                    end
                end
                svgFileName=this.getSvgFile(idx);
                ToolBarMethods.addButton(tb_, 'svg.png', ...
                    'Create SVG File', @(h,e)svg1D(this, fig1D, ...
                    svgFileName));
                Gui.CascadeFigs(this.figs1D.values, true);
                
            else
                set(fig1D, 'visible', 'off', 'CloseRequestFcn',...
                    @(h, e)hush2(key, h));
            end
            fp=this.getPngPath(key);
            if ~exist(fp, 'file')
                disp(get(fig1D, 'outerposition'));
                saveas(fig1D, fp);
                %Gui.SavePng(fig1D, fp, 400);
            end
            
            function hush2(key, h)
                this.figs1D.remove(key);
                delete(h);
                if ishandle(this.fig)
                    figure(this.fig);
                end
            end
        end
        
        function [ids, cmp, N]=getPinComponent(this, idx, start, stop, jMenu, ids, N)
            if nargin<6
                [ids,N]=this.getSelectedLeaves(true);
            end
            isMenu=nargin>=5;
            cmp=[];
            doubleClickTxt=[this.gt.app.smallStart ...
                ' (<b><i>double click</i></b>)' this.gt.app.smallEnd];
            if N>1&&isMenu
                Gui.NewMenuItem(jMenu, ['<html>SubsetHighlighter: ' num2str(N) ...
                    ' subsets' doubleClickTxt '</html>'], ...
                    @(h,e)highlightSubsetsNonModal(this, ids), ...
                    'pinFlashlightTransparent.png');
                jMenu.addSeparator;
            else
                if idx==0
                    idx=this.getIdx(ids{1});
                    gid_=ids{1};
                else
                    gid_=this.nodes{idx};
                    if length(gid_)>1
                        cmp=[];
                        return;
                    end
                    gid_=gid_{1};
                end
                newPinColor=get(this.plotHs(idx), 'markerFaceColor');
                if isMenu
                    bckg2=[1 1 1];
                else
                    bckg=javax.swing.UIManager.get('ToolTip.background');
                    bckg2=[bckg.getRed/255 bckg.getGreen/255 bckg.getBlue/255];
                end
                on=CellHighlighter.Is(this.gt, gid_);
                if on
                    if ~isMenu
                        lblPin='Stop !';
                    else
                        lblPin=['<html>Stop highlighting' doubleClickTxt '</html>'];
                    end
                    fncHighlight=@(h,e)highlight(this, gid_, [], isMenu);
                else
                    if ~isMenu
                        lblPin='High-<br>light';
                    else
                        lblPin=['<html>Start highlighting' doubleClickTxt '</html>'];
                    end
                    fncHighlight=@(h,e)highlight(this, gid_, newPinColor, isMenu);
                end
                lbl=[start lblPin stop];
                if ~isMenu
                    [~, cmp]=Gui.ImageLabel(lbl, 'pinFlashlight.png',...
                        [], fncHighlight);
                    cmp.setForeground(java.awt.Color.blue);
                else
                    cmp=Gui.NewMenuItem(jMenu, lbl, fncHighlight);
                    jMenu.addSeparator;
                end
                if sum(newPinColor)>2.8
                    Gui.SetIcon2(cmp, 'pinFlashlight.png', ...
                        [1 1 0], [1 1 0], [1 1 1], bckg2);
                else
                    Gui.SetIcon2(cmp, 'pinFlashlight.png', ...
                        [1 1 0], newPinColor, [1 1 1], bckg2);
                end  
            end
        end
        
        function bp=getTipPanel(this, idx, showCreateBtn)
            btnPnl=Gui.Panel;
            start=['<html><center>&nbsp;<u>' this.app.smallStart];
            stop=[this.app.smallEnd '</u></center></html>'];
            [ids, btnPin]=this.getPinComponent(idx, start, stop);
            if ~isempty(btnPin)
                btnPnl.add(btnPin);
            end
            lbl=[start '<u>High<br>res 1D' stop];
            [~, btnHi]=Gui.ImageLabel(lbl, 'pseudoBarHi.png', ...
                [],@(h,e)browseSelectedHiRes1DPF(this, idx));
            btnHi.setForeground(java.awt.Color.blue);
            btnPnl.add(btnHi);
            lbl=[start '<u>Low<br>res 1D' stop];
            [~, btnLo]=Gui.ImageLabel(lbl, 'pseudoBarHi.png', ...
                [],@(h,e)browseSelectedLoRes1DPF(this, idx));
            btnLo.setForeground(java.awt.Color.blue);
            btnPnl.add(btnLo);            
            lbl=[start '<u>Gate<br>seq.' stop];
            [~, btnGs]=Gui.ImageLabel(lbl, 'comicStrip.png', ...
                [],@(h,e)seeComicStrip(this, ids, idx));
            btnGs.setForeground(java.awt.Color.blue);
            btnPnl.add(btnGs);
            extra=false;
            if extra
                lbl=[start '<u>Plot<br>editor' stop];
                [~, btnGs]=Gui.ImageLabel(lbl, 'eye.gif', ...
                    [], @(h,e)plotEditor(this, ids, idx));
                btnGs.setForeground(java.awt.Color.blue);
                btnPnl.add(btnGs);
                lbl=[start '<u>Tree<br>position' stop];
                [~,btnTree]=Gui.ImageLabel(lbl, 'tree.png', ...
                    [], @(h,e)treePosition(this, ids, idx));
                btnPnl.add(btnTree);
            end
            if showCreateBtn
                btnCreate=Gui.NewBtn(['<html>' this.app.smallStart ...
                    'Create 1D <br>PathFinder(s)' ...
                    this.app.smallEnd '</html>'], ...
                    @(h,e)createImg(this),'', 'pseudoBarLo.png');
                bp=Gui.AddTipImgCheckBox(this.app, btnPnl, this.lastTipFnc, ...
                    'North', 'South', btnCreate);                
            else
                bp=Gui.AddTipImgCheckBox(this.app, btnPnl, this.lastTipFnc);
            end
        end
        
        function createImg(this, idxs, ask, which)
            if nargin<4
                which=0;
                if nargin<3
                    ask=true;
                    if nargin<2
                        idxs=[];
                    end
                end
            end
            missing=this.countMissingPngs(idxs, false, which);
            nMissing=length(missing);
            if nMissing>1 && ask
                str=String.Pluralize2('missing low res 1D <br>PathFinder image',...
                    nMissing);
                if ~askYesOrNo(Html.WrapC(['Create the ' str '?']), ...
                        str, 'south', true, 'missingLowRes1DPF')
                    return;
                end
            end
            if nMissing>0
                Gui.setEnabled(this.fig, false);
                this.tb.setEnabled(false);
                try
                    this.countMissingPngs(missing, true, which);
                catch ex
                    disp(ex);
                end
                Gui.setEnabled(this.fig, true, true);
                this.tb.setEnabled(true);
            end
        end

        function ids=getPhenogramId(this, ids, idx)            
            if ischar(ids)
                ids{1}=ids;
            end
            N=length(ids);
            if N>1 && idx>0
                if length(this.gtMap)>=idx
                    id=this.gtMap{idx};
                    if ~isempty(id)
                        ids={id};
                        return;
                    end
                end
            end
            ids=[];
        end
        
        function [ids, N]=resolveIds(this, ids, idx)            
            if ischar(ids)
                ids{1}=ids;
            end
            N=length(ids);
            if N>1 && idx>0
                if length(this.gtMap) >= idx && ~isempty(this.gtMap{idx})
                    ttl=String.Pluralize2('merged gate', N);
                    [yes, cancelled]=askYesOrNo(['Use boolean gate of ' ...
                        ttl '?'],ttl, 'North', false, 'usePhenoBoolean');
                    if yes
                        ids={this.gtMap{idx}};
                        N=1;
                    end
                    if cancelled
                        N=0;
                    end
                end
            end
        end
        
        function treePosition(this, ids, idx)
            [ids, N]=this.resolveIds(ids,idx);
            figure(this.gt.fig);
            ok=this.gt.gtt.areGatingTreeSelectionsInSync;
            if ok
                this.gt.gtt.keepGatingTreeSelectionsInSync(false);
            end
            
            this.gt.ensureVisible(ids{1}, true);
            for i=2:N
                this.gt.ensureVisible(ids{i}, true, true, true);
            end
            this.gt.syncGatingTreeTableSelections(true);
            if ok
                this.gt.gtt.keepGatingTreeSelectionsInSync(ok);
            end
            this.gt.selectHistory.remember;
        end
        
        
        function highlightSubsetsNonModal(this, ids)
            if nargin<2
                ids=this.getSelectedLeaves(true);
            end
            sid=this.gt.getSampleNode(ids{1});
            fld=FlashLightDlg({this.gt}, {sid}, {ids});
            N=length(ids);
            dscs=cell(1,N);
            labels=cell(1,N);
            colors=cell(1,N);
            for i=1:N
                idx=this.getIdx(ids{i});
                [~, dscs{i}, ~, labels{i}]=this.unpackNode(idx, true);
                colors{i}=get(this.plotHs(idx), 'markerFaceColor');
            end
            fld.setLabels({labels});
            fld.setPinLabels({dscs});
            fld.setColors({colors});
            fld.show('south');
            fld.askAboutUnsetColors('phenogram');
        end
        
        function highlightSubsetsModal(this)
            [selectedIds, N]=this.getSelectedLeaves(true);
            if N==0
                return;
            end
            chosen=[];
            for i=1:N
                on=CellHighlighter.Is(this.gt, selectedIds{i});
                if on
                    chosen(end+1)=i-1;
                end
            end
            allChosen=[];
            [allIds, N2]=this.getSelectedLeaves(false);
            for i=1:N2
                on=CellHighlighter.Is(this.gt, allIds{i});
                if on
                    allChosen(end+1)=i-1;
                end
            end
            [~, ~, idxs, ~, cancelled, finalIds]=this.chooseIds(selectedIds, ...
                'Subset highlighting', ...
                'highlighting (AKA back-gating)', false, chosen, allChosen);
            if cancelled
                return;
            end
            figure(this.fig);
            pu=PopUp('Adjusting subset highlighting');
            N=length(finalIds);
            pu.initProgress(N);
            this.gt.manyPins=true;
            for i=1:N
                gid_=finalIds{i};
                dsc=this.gt.tp.getDescription(gid_);
                idx=this.getIdx(gid_);
                turnOn=~isempty(find(idx==idxs, 1));
                if ~turnOn
                    newPinColor=[];
                else
                    newPinColor=get(this.plotHs(idx), 'markerFaceColor');
                end
                this.highlight(gid_, newPinColor, true);
                pu.incrementProgress;
            end
            this.gt.manyPins=false;
            this.gt.notifyHighlightChanges;                
            pu.close;
        end
        
        function [clrs, edgeClrs, lineWidths]=getLeafPlotElements(this)
            N=this.nLeafs;
            clrs=zeros(N,3);
            edgeClrs=zeros(N,3);
            lineWidths=zeros(1, N);
            for i=1:N
                clrs(i,:)=get(this.plotHs(i), 'MarkerFaceColor');
                edgeClrs(i,:)=get(this.plotHs(i), 'MarkerEdgeColor');
                lineWidths(i)=get(this.plotHs(i), 'LineWidth');
            end
        end
        
        function [ids, N, idxs]=getSelectedLeaves(this, useSelected)
            if useSelected
                nSelected=this.selected.size;
            else
                nSelected=0;
            end
            if nSelected==0
                nSelected=length(this.leafNames);
                selectedIdxs=this.sortXyAndAx('Phenogram order', 1);
            else
                selectedIdxs=this.selected.toArray;
            end
            if selectedIdxs(1)==0
                nSelected=length(this.leafNames);
                selectedIdxs=this.sortXyAndAx('Phenogram order', 1);
            end                
            ids=this.nodes{selectedIdxs(1)};
            for i=1:nSelected
                ids=[ids this.nodes{selectedIdxs(i)}];
            end
            ids=StringArray.RemoveDuplicates(ids, java.util.LinkedHashSet);
            N=length(ids);
            if nargout>2
                idxs=zeros(1, N);
                for i=1:N
                    idxs(i)=this.getIdx(ids{i});
                end
            end            
        end
        
        function [idxs, N]=getSelectedBranches(this)
            nSelected=this.selected.size;
            selectedIdxs=this.selected.toArray;
            idxs=[];
            for i=1:nSelected
                ids=this.nodes{selectedIdxs(i)};
                if length(ids)>1
                    idxs(end+1)=selectedIdxs(i);
                end
            end 
            N=length(idxs);
        end
 
        function branchCnt=getBranchCount(this)
            nSelected=this.selected.size;
            selectedIdxs=this.selected.toArray;
            branchCnt=0;
            for i=1:nSelected
                idx=selectedIdxs(i);
                if idx<1
                    continue;
                end
                ids=this.nodes{idx};
                cnt=length(ids);
                if cnt>1
                    branchCnt=branchCnt+1;
                end
            end 
        end
 
        function seeComicStrip(this, ids, idx)
            if nargin<2
                ids=this.getSelectedLeaves(true);
            else
                ids=this.resolveIds(ids,idx);
            end
            [ids, cbBrowse, idxs, numCols]=this.chooseIds(...
                ids, 'Gate sequences', 'gating sequences', true);
            N=length(ids);
            inBrowser=cbBrowse.isSelected;
            if N==0
                return;
            end
            if ~inBrowser
                BasicMap.Global.set('showPopUp', 'false');
                try
                    for i=1:N
                        jd=this.gt.seeComicStrip(ids{i}, 0, 'north++');
                        Gui.SetJavaVisible(jd);
                        Gui.CascadeFromNorth(jd, i);
                    end
                catch ex
                    disp(ex);
                end
                BasicMap.Global.set('showPopUp', 'true')
            else
                figure(this.fig);
                pu=PopUp('Preparing gating sequences');
                pu.initProgress(N);

                html=['<html><h2>' String.Pluralize2('gating sequence', ...
                    N) '</h2><h3>Picked in phenogram</h3><table>'];
                if numCols==0
                    td='<td><hr>';
                else
                    td=['<td colspan="' num2str(numCols) ...
                        '" bgcolor="#FFFFDD"><hr>'];
                end
                for i=1:N
                    [~, ~, ~, ttl]=this.unpackNode(idxs(i), true);
                    html=[html '<tr>' td ttl  '</td></tr><tr><td>' ...
                        this.gt.describe(ids{i}, false, true, false, ...
                        'blue', true, true, true, true, true, numCols)...
                        '</td></tr>'];
                    pu.incrementProgress;
                end
                html=[html '</table></html>'];
                Html.Browse(html);
                pu.close;
            end
        end
            
        function plotEditor(this, ids, idx)
            [ids, N]=this.resolveIds(ids,idx);
            if N>0
                figure(this.fig);
                pu=PopUp(['Opening PlotEditor ' String.Pluralize2(...
                    'window', N)], 'center', 'Opening', true, true);
                pu.initProgress(N);
                for i=1:N
                    GatingTree.NavigateTo(this.gt, ids{i}, 'new');
                    drawnow;
                    if pu.cancelled
                        break;
                    end
                    pu.incrementProgress;
                end
                pu.close;
            end
        end
        
        function svg1D(this, fig1D, fileName)
            Exporter.SaveVg(fig1D, this.gt, ...
                sprintf('QFM_%s_allD', fileName));
        end

        function idx=findNearestNode(this, cp, isMotion)
            x=cp(1, 1);
            y=cp(1, 2);
            xy=NdSee.Normalize(this.ax, [x y]);
            x=xy(1);
            y=xy(2);
            near2=pdist2(this.xyNormalized, xy);
            [minDist,idx_]=min(near2);
            idx=0;
            if isMotion
                limit=.155;
            else
                limit=.03;
            end
            if idx_>0
                if  minDist>limit
                    this.app.closeToolTip;
                    this.lastNear=0;
                elseif ~isMotion || this.lastNear ~= idx_
                    idx=idx_;
                    this.lastNear=idx_;
                end
            end
            %fprintf('%d not found\n', idx);
        end
        
        function highlight(this, gid_, newPinColor, isMenu)
            busy=Gui.ShowBusy(this.fig,'<h3>Updating highlights<hr></h3>');
            this.tb.setEnabled(false);
            gt_=this.gt;
            if ~isempty(newPinColor)
                if sum(newPinColor)>2.8
                    CellHighlighter.Start(gt_, gid_);
                else
                    CellHighlighter.Start(gt_, ...
                        gid_, newPinColor);
                end
            else
                CellHighlighter.Stop(gt_, gid_, false);
            end
            this.tb.setEnabled(true);
            Gui.HideBusy(this.fig, busy);
            if ~isMenu
                feval(this.lastTipFnc);
            end
        end
        
    end
    
    methods(Static)
        
        function [fig, ax, nodes,xyNormalized, figs1D, nodeNames, ...
                xyIn, textHs, plotHs, prop]=New(qf, ttl, props, propPfx, ...
                colors, edgeColors, lineWidths, tNames, qfObj)
            app=BasicMap.Global;
            fig=[];
            ax=[];
            nodes={};
            xyNormalized=[];
            xyIn=[];
            textHs=[];
            plotHs=[];
            figs1D=[];
            showIDs=app.is(QfTree.PROP_IDS, false);
            showDs=app.is(QfTree.PROP_D, false);
            try
            figs1D=MatLabMap;
            catch ex
            end
            priorFig=get(0, 'currentFig');
            nColors=size(colors,1);
            nEdgeColors=size(edgeColors, 1);
            MIN=7; % 10 pixels ( sqrt(100)==10 )
            MAX=60; % 40 pixels
            nSubsets=length(qf.tIds);
            nLeaves=qf.numLeaves;
            nBranches=length(qf.branchNames);
            nLeafs=length(tNames);
            
            cntr=edu.stanford.facs.swing.Counter;
            phyNames=cell(1, nLeaves+nBranches);
            for iNode=1:nLeaves
                key=tNames{iNode};
                cntr.count(key);
                n=cntr.getCount(key);
                if n==1
                    phyNames{iNode}=key;
                else
                    phyNames{iNode}=String.PadEnding(key, n, ' ');
                end
            end
            if ~showIDs
                for iNode=iNode+1:nLeaves+nBranches
                    phyNames{iNode}=num2str(iNode);
                end
            else
                for iNode=iNode+1:nLeaves+nBranches
                    phyNames{iNode}=['                         ' num2str(iNode)];
                end
            end
            nodeNames=[tNames qf.branchNames];
            nNodes=length(nodeNames);
            prop=[propPfx '.' num2str(nNodes)];
            try
                pt=phytree(qf.phyTree, qf.nodeQfs', phyNames);
            catch ex
                msg(Html.WrapHr(['No phenogram can be done.<br><br>'...
                    'Ensure that you have installed the <br>'...
                    'Bioinformatics toolbox from MathWorks....']));
                ex.getReport
                fig=[];
                return;
            end
            typeValue='square';
            plot(pt, 'type', typeValue, 'BranchLabels', true, ...
                'LeafLabels', true, 'TerminalLabels', false);
            fig=gcf;
            set(fig, 'visible', 'off', 'UserData', 'phenogram');
            op=get(fig, 'OuterPosition');
            
            ax=get(fig, 'currentaxes');
            hold(ax, 'on')
            axis(ax, 'off')
            ttl{1}=[num2str(nSubsets) ' subsets in ' ttl{1}];
            title(ax, ttl);
            set(ax, 'Xticklabel', []);
            try
                prefix=CytoGate.Abbreviate('Phenogram');
            catch
                prefix='QF tree';
            end
            set(fig, 'Name', [prefix qf.distanceType ...
                ' distance, ' num2str(length(qf.columnNames)) 'D'],...
                'NumberTitle', 'off', 'CloseRequestFcn', @(h, e)hush);
            txtObjs=findall(fig, 'type', 'text');
            nTxtObjs=length(txtObjs);
            strsArray=StringArray(phyNames);
            op=str2num(props.get(prop,'0'));
            if length(op)==4
                op=Gui.RepositionOnSameScreenIfRequired(op);
                set(fig, 'OuterPosition', op);
            end
            set(fig, 'color', 'white');
            treeSz_=qf.treeSz;
            xl=xlim(ax);
            nudgeX=(xl(2)-xl(1))/50;
            nudgeXb=(xl(2)-xl(1))/55;
            yl=ylim(ax);
            nudgeY=(yl(2)-yl(1))/50;
            fs1='';
            fs2='';
            fs3=QfTree.FONT_SIZE_QF_SCORE;
            nodes=cell(1,strsArray.N);
            xyIn=zeros(strsArray.N, 2);
            textHs=zeros(1, strsArray.N);
            plotHs=zeros(1, strsArray.N);
            for i=1:nTxtObjs
                txtObj=txtObjs(i);
                if isequal('off', get(txtObj, 'visible'))
                    continue;
                end
                txt=get(txtObj, 'String');
                idx=strsArray.indexOf(txt);
                if idx>0
                    isRoot=idx==nNodes;
                    P=get(txtObj, 'Position');
                    sz=qf.nodeSzs(idx);
                    isLeaf=idx<=nLeaves;
                    if nNodes>25
                        if nLeafs>40
                            if isLeaf
                                fs1='\fontsize{8}';
                                fs2='\fontsize{6}';
                            else
                                fs1='\fontsize{8}';
                                fs2='\fontsize{6}';
                            end
                        elseif nLeafs>25
                            if isLeaf
                                fs1='\fontsize{9}';
                                fs2='\fontsize{8}';
                            else
                                fs1='\fontsize{9}';
                                fs2='\fontsize{9}';
                            end
                        end
                    else
                        if isLeaf
                            fs1='\fontsize{11}';
                            fs2='\fontsize{10}';
                        else
                            fs1='\fontsize{11}';
                            fs2='\fontsize{11}';
                        end
                    end
                    if isLeaf
                        ids={num2str(qf.tIds(idx))};
                        qfScore=0;
                    else
                        if showIDs
                            txt=nodeNames{idx};
                            ids=strsplit(txt, ' ');
                        else
                            txt=nodeNames{idx};
                            ids=strsplit(txt, ' ');
                            txt='';
                            %ids={};
                        end
                        qfScore=qf.branchQfs(idx-nLeaves);
                    end
                    txt=strrep(txt, '\bf', ' \bf \color{black}');
                    freq=sz/treeSz_;
                    if isLeaf
                        perc=[fs2 '  \bf\color[rgb]{0 0 .8}'...
                            String.encodePercent(sz, treeSz_)];
                    else
                        perc=[fs2 '\bf\color[rgb]{0 0 .8}'...
                            String.encodePercent(sz, treeSz_) '  '];
                    end
                    %perc='';
                    if showDs
                        sqf=[fs2 '\color{magenta}D=' ...
                            String.encodeRounded(qfScore, 2)];
                    else
                        sqf='';
                    end
                    if showDs || showIDs || isLeaf
                        if ~isLeaf
                            if showIDs
                                toks=ids;
                                nToks=length(toks);
                            else
                                nToks=0;
                            end
                            MAX_TOK=3;
                            if nToks>MAX_TOK
                                toks{end+1}=[' ' perc ' '];
                                toks{end+1}=sqf;
                                nToks=nToks+2;
                                str='';
                                strs={};
                                for j=1:nToks
                                    tok=toks{j};
                                    at4=mod(j, MAX_TOK)==0;
                                    if at4
                                        strs{end+1}=[str tok];
                                        str='';
                                    else
                                        str=[str tok ' '];
                                    end
                                end
                                if ~at4
                                    strs{end+1}=str;
                                end
                            else
                                strs={fs1 txt, ['\color{black}' perc ' ' sqf]};
                            end
                            fontColor=[.35 .35 .4];
                        else
                            if showIDs
                                strs={fs1 txt, [num2str(qf.tIds(idx))...
                                    ' ' perc]};
                            else
                                if nSubsets>20
                                    strs={ [fs1 txt perc]};
                                else
                                    strs={fs1 txt, perc};
                                end
                            end
                            fontColor=[.4 .4 .45];
                        end
                        set(txtObj, 'FontName', QfTree.FONT_NAME, ...
                            'Color', fontColor, ...
                            'FontWeight', 'bold', ...
                            'String', strs, 'Interpreter', 'tex');
                    else
                        set(txtObj, 'FontName', QfTree.FONT_NAME, 'String', perc, 'Interpreter', 'tex');
                    end
                    ms=floor((freq+(.2*(1-freq)))*MAX);
                    if ms<MIN
                        ms=MIN;
                    end
                    ratio=ms/MAX;
                    if isLeaf
                        symbol='o';
                        nudge=nudgeX;
                        %set(txtObj, 'Position', ...
                        %    [P(1)+(nudgeX*ratio)+(nudgeX*.35) P(2)-nudgeY P(3)]);

                    else
                        ratio=ratio+(.8*(1-ratio));
                        nudge=0-(nudgeXb*ratio);
                        
                        ratio=ms/MAX;
                        if ratio>.97
                            ratio=2;
                        elseif ratio>.6
                            ratio=1.4;
                        end
                        if isRoot
                            symbol='o';
                            %ms=.7*MAX;
                            if ~showIDs %&& ~showQFs
                                txtPos=[P(1)+(nudgeX*2.2)+(nudgeX*ratio) ...
                                    P(2)-(.2*nudgeY) P(3)];
                            else
                                txtPos=[P(1)+(nudgeX*7)+(nudgeX*ratio) ...
                                    P(2) P(3)];
                            end
                            ms=ms*.9;
                        else
                            symbol='o';
                            if ~showIDs  %&& ~showQFs
                                txtPos=[P(1)+(nudgeX*1.8)+(nudgeX*ratio) ...
                                    P(2)+nudgeY P(3)];
                            else
                                txtPos=[P(1)-(nudgeX*ratio) P(2)-nudgeY P(3)];
                            end
                            ms=ms*.9;
                        end
                        if ~showIDs
                            %set(txtObj, 'Position', txtPos, 'HorizontalAlignment', 'left');
                        else
                            set(txtObj, 'Position', txtPos);
                        end
                    end
                    try
                        draggable(txtObj, 'endFcn', @(h)setMouse(qfObj));
                    catch ex
                    end
                    textHs(idx)=txtObj;
                    if isLeaf
                        if idx<=nColors && lineWidths(idx)~=0
                            faceColor=colors(idx,:);
                            if sum(faceColor)>2.95
                                symbol='s';
                            end
                        else
                            faceColor=Gui.HslColor(idx, nNodes);
                        end
                        if idx<=nEdgeColors && lineWidths(idx)~=0
                            edgeColor=edgeColors(idx, :);
                            lineWidth=lineWidths(idx);
                        else
                            edgeColor='blue';
                            lineWidth=1;
                        end
                    else
                        if isRoot
                            faceColor='none';
                            lineWidth=2;
                        else
                            ratio2=(idx-nLeaves)/nBranches;
                            rgb=.45+(ratio2*.4);
                            faceColor=[rgb rgb rgb];
                            lineWidth=1;
                            edgeColor='black';
                        end
                    end
                    plotHs(idx)=plot(ax, P(1)-nudge, P(2), symbol, ...
                        'markerSize', ms, ...
                        'markerFaceColor', faceColor,...
                        'markerEdgeColor', edgeColor, ...
                        'lineWidth', lineWidth);
                    xyIn(idx, :)=[P(1)-nudge P(2)];
                    nodes{idx}=ids;
                    uistack(txtObj, 'top');
                end
            end
            xyNormalized=Gui.Normalize(ax, xyIn);
            function hush
                op=get(fig, 'OuterPosition');
                props.set(prop, num2str(op));
                if ~isempty(figs1D)
                    names=figs1D.names;
                    nNames=length(names);
                    for k=1:nNames
                        f=figs1D.get(names{k});
                        delete(f);
                    end
                end
                delete(fig);
                if ~isempty(priorFig) && ishandle(priorFig)
                    if isequal('on', get(priorFig, 'visible'))
                        figure(priorFig);
                    end
                end
            end        
            set(ax, 'units', 'normalized', 'Position', [.03 .01 .95 .925]);
        end
        
        function this=Load(filePath, ttl, gt, gid, mp)
            load(filePath, 'QF');
            if nargin<2
                ttl=QF.ttl;
            end
            if nargin>2
                names=QfHiDM.GetNames(QF, gt.tp);
                [clrs, edgeClrs, lineWidths]=QfHiDM.GetColors(QF, gt.tp);
                this=QfTree(QF, ttl, mp, '', clrs, edgeClrs, lineWidths,...
                    names);
            else
                this=QfTree(QF, ttl, BasicMap.Global, '', QF.colors, ...
                    QF.edgeColors, QF.lineWidths,...
                    QF.names);
            end
            if nargin>2
                if ~isempty(this.fig)
                    this.setGt(gt, gid);
                end
            end
        end
        
        function sz=Size(this, ids)
            sz=0;
            N=length(ids);
            for i=1:N
                id=ids(i);
                idx=find(this.tIds==id,1);
                if ~isempty(id)
                    sz=sz+this.tSizes(idx);
                end
            end
        end
        
        function [qf, qfTree]=RunWithGt(gt, gid, fcs, fcsIdxs, visible, pu)
            qf=[];
            qfTree=[];
            if isempty(fcsIdxs)
                return;
            end
            newPu=nargin<6 || isempty(pu);
            if newPu
                pu=PopUp('Generating phenogram', 'north', 'Note...', false, ...
                    true, 'genieSearch.png');
                pu.setTimeSpentTic(tic);
            end
            [finalData, teachCols,teachColNames, teachIds, teachData, teachCompData]...
                =QfHiDM.GetRequiredData2(fcs, fcsIdxs, gt, gid, pu, visible);
            if isempty(finalData)
                if newPu
                    pu.close;
                end
                return;
            end
            try           
                maxDeviantParameters=gt.multiProps.getNumeric(...
                    QfHiDM.PROP_MAX_DEVIANT_PARAMETERS, 0);
                if QfHiDM.F_MEASURE_OPTIMIZE
                    matchStrategy=gt.multiProps.getNumeric(...
                        QfHiDM.PROP_F_MEASURE_OPTIMIZE, 0);
                else
                    matchStrategy=0;
                end
                mergeStrategy=gt.multiProps.getNumeric(QfHiDM.PROP_MERGE_STRATEGY, 1);
                bins=0;
                qfTreeCode=1;
                [done, qf]=NdTreeGater.Qf(teachColNames, gt, gid, gt, [], ...
                    teachData, teachCompData, teachColNames, teachCols, ...
                    teachIds, pu, [], [], qfTreeCode, bins,...
                    QfHiDM.BIN_STRATEGY, maxDeviantParameters, matchStrategy, ...
                    mergeStrategy);
                if done
                    [~,fig_, qfTree]=QfTree.ViewWithGt(qf, gt, gid, false);
                    if ~isempty(fig_)
                        gt.otherFigs{end+1}=fig_;
                        if visible
                            set(fig_, 'visible', 'on');
                        end
                    end
                end
            catch ex
                ex.getReport
            end
            if newPu
                pu.close;
            end
        end
        
        function [rc,fig, qfTree]=ViewWithGt(qf, gt, gid, usePrior)
            fig=[];
            qfTree=[];
            gtp=gt.tp;
            ttl=GatingTree.DescribeGateSample(gtp, gid, true);
            mp=MultiProps(BasicMap.Global, gtp);
            str=gtp.getNode(gid, MatchInfo.PROP_LAST_TITLE);
            if ~isempty(str)
                ttl=String.DecodeStrs(str);
            end
            if usePrior
                [filePath, exists]=GatingTreeTable.GetQfTreeFile(gtp, gid);
                if exists
                    [~, yes,cancelled]=questDlg(struct('msg', ...
                        'Use prior phenogram ?', ...
                        'where', 'north east+', 'remember', 'phenogramReuse'),...
                        'Please confirm...', 'Yes', 'No, rebuild',  'Cancel', 'Yes');
                    if cancelled
                        rc=-2;
                        return;
                    end
                    if ~yes
                        rc=0;
                        [filePath, exists]=GatingTreeTable.GetQfTreeFile(gtp, gid);
                        if exists
                            delete(filePath);
                        end
                        return;
                    end
                    
                    qfTree=QfTree.Load(filePath,ttl, gt, gid, mp);
                    fig=qfTree.fig;
                    rc=-1;
                end
            end
            if ~isempty(qf)
                rc=1;
                [clrs, edgeClrs, lineWidths]=QfHiDM.GetColors(qf, gtp);
                qf.prepareMedianMarkerExpressions(gt, gid);
                qfTree=qf.viewTree(ttl, mp, clrs,edgeClrs, lineWidths, ...
                    QfHiDM.GetNames(qf, gtp));
                fig=qfTree.fig;
                if ~isempty(qfTree.fig)
                    qfTree.setGt(gt, gid);
                    if usePrior
                        qfTree.save(filePath, gt, gid)
                    end
                end
            else
                rc=0;
            end
        end
    end
end
