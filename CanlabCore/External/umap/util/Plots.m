%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

classdef Plots < handle
    properties(Constant)
        DOING_ALL='doingAll';
        DEBUGGING=false;
    end
    
    properties(SetAccess=private)
        is3D;
        legendHs;
        mdns;
        stds;
        legendH;
        otherPlots={};
    end
    
    properties
        Hs;
        otherHs;
        cnts;
        names;
        N;
        distFactor;
        mx;
        mxI;
        mn;
        mnI;
        ax;
        javaLegend;
        otherPlotMap;
    end
    
    methods
        function this=Plots()
            this.Hs=[];
        end
        
        function setHs(this, Hs)
            this.Hs=Hs;
            this.N=length(this.Hs);
            this.cnts=zeros(1, this.N);
            if isempty(this.names) 
                this.names=cell(1, this.N);
                for i=1:this.N
                    this.cnts(i)=length(get(this.Hs(i), 'XData'));
                    this.names{i}=get(this.Hs(i), 'DisplayName');
                end
            else
                for i=1:this.N
                    this.cnts(i)=length(get(this.Hs(i), 'XData'));
                end
            end
            this.setMinMax;
            this.ax=get(Hs(1), 'Parent');
        end
        
        function setMinMax(this)
            [this.mx, this.mxI]=max(this.cnts);
            [this.mn, this.mnI]=min(this.cnts);
        end
        
        function setOtherHs(this, Hs)
            this.otherHs=Hs;
        end
        
        function setLegendHs(this, Hs)
            this.legendHs=Hs;
        end
        
        function initStats(this)
            if isempty(this.mdns)
                if isempty(this.Hs)
                    return;
                end
                [this.mdns, this.stds, this.is3D]=Plots.Stats(this.Hs);
                if this.is3D
                    this.distFactor=6;
                else
                    this.distFactor=3;
                end
            end
        end        
        
        function setNames(this, names)
            this.names=names;
        end
        
        function names=getNames(this)
            names=this.names;
        end
        
        function toggleVisibility(this, idx)
            if ~isempty(this.otherHs)
                R=size(this.otherHs,1);
                for r=1:R
                    plotH=this.otherHs(r, idx);
                    if strcmp('off', get(plotH, 'visible'))
                        set(plotH(r), 'visible', 'on');
                    else
                        set(plotH(r), 'visible', 'off');
                    end
                end
            end
        end
        
        function clear(this)
            try
                N_=length(this.Hs);
                for i=1:N_
                    delete(this.Hs(i));
                end
                N_=length(this.otherHs);
                for i=1:N_
                    delete(this.otherHs(i));
                end
            catch ex
                %ex.getReport
            end
        end
        
        function addOtherPlots(this, other)
            this.otherPlots{end+1}=other;
        end
        
        function flashOtherPlots(this, idx, on)
           if ~isempty(this.otherPlots)
               name=edu.stanford.facs.swing.Basics.RemoveXml(String.RemoveTex(this.names{idx}));
               lIdx=String.LastIndexOf(name, ' training ');
               if lIdx>0
                   name=name.substring(0,lIdx-1);
               end
               nOthers=length(this.otherPlots);
               for i=1:nOthers
                   other=this.otherPlots{i};
                   nNames=length(other.names);
                   for j=1:nNames
                       if strcmp(name, edu.stanford.facs.swing.Basics.RemoveXml(String.RemoveTex(other.names{j})))
                           perc=other.cnts(j)/sum(other.cnts);
                           otherH=other.Hs(j);
                           if nargin<3
                               if isequal(get(otherH, 'visible'), 'on')
                                   set(otherH, 'visible', 'off');
                                   times=3;
                               else
                                   set(otherH, 'visible', 'on');
                                   times=4;
                               end
                           else
                               if ~on
                                   set(otherH, 'visible', 'off');
                                   times=3;
                               else
                                   set(otherH, 'visible', 'on');
                                   times=4;
                               end
                           end
                           Plots.Flash(otherH, perc, times);
                           break;
                       end
                   end
               end
           end
           if nargin>2
               this.setOtherPlotMapVisible(idx, on)
           end
        end
        
        function yes=isInOtherPlotMap(this, idx)
            if ~isempty(this.otherPlotMap)
                H=this.Hs(idx);
                yes=this.otherPlotMap.isKey(H);
            else
                yes=false;
            end
        end
        
        function setOtherPlotMapVisible(this, idx, on)
            if ~isempty(this.otherPlotMap)
                H=this.Hs(idx);
                if on
                    vis='on';
                else
                    vis='off';
                end
                if this.otherPlotMap.isKey(H)
                    v=this.otherPlotMap(H);
                    if iscell(v)
                        v=v{1};
                    end
                    if any(v<=0)
                        v(v<0)=v;
                    end
                    set(v, 'visible', vis);
                end
            end
        end
        
        function setOtherPlotVisible(this, idx, on)
           if ~isempty(this.otherPlots)
               name=edu.stanford.facs.swing.Basics.RemoveXml(String.RemoveTex(this.names{idx}));
               lIdx=String.LastIndexOf(name, ' training ');
               if lIdx>0
                   name=name.substring(0,lIdx-1);
               end
               nOthers=length(this.otherPlots);
               for i=1:nOthers
                   other=this.otherPlots{i};
                   nNames=length(other.names);
                   for j=1:nNames
                       if strcmp(name, edu.stanford.facs.swing.Basics.RemoveXml(String.RemoveTex(other.names{j})))
                           perc=other.cnts(j)/sum(other.cnts);
                           otherH=other.Hs(j);
                           if ~on
                               set(otherH, 'visible', 'off'); 
                           else
                               set(otherH, 'visible', 'on');
                           end
                       end
                   end
               end
           end
           this.setOtherPlotMapVisible(idx, on)
        end
    end
    
    methods(Static)
        
        function [mdns, stds, is3D]=Stats(plotHs)
            N=length(plotHs);
            mdns=zeros(N,3);
            stds=zeros(N,3);
            is3D=true;
            for ii=1:N
                plotH=plotHs(ii);
                if ~ishandle(plotH)
                    return;
                end
                x_=get(plotH, 'XData');
                y_=get(plotH, 'YData') ;
                z_=get(plotH, 'ZData');
                if isempty(z_)
                    [RR,CC]=size(get(plotHs(ii), 'XData'));
                    z_=ones(RR, CC);
                    is3D=false;
                end
                d=[x_; y_; z_]';
                mdns(ii,:)=median(d);
                stds(ii,:)=std(d);
            end
        end
        
        function Flash(plotH, perc, times)
            vi=get(plotH, 'Visible');
            ms=get(plotH, 'MarkerSize');
            if perc<.02
                set(plotH, 'MarkerSize', ms+7)
            elseif perc<.05
                set(plotH, 'MarkerSize', ms+5)
            elseif perc<.1
                set(plotH, 'MarkerSize', ms+3);
            elseif perc<.2
                set(plotH, 'MarkerSize', ms+2);
            end
            for i=1:times
                %MatBasics.RunLater(@(h,e)go, .1);
                go
            end
            %MatBasics.RunLater(@(h,e)done, .51);
            done
            function go
                set(plotH, 'visible', 'off');
                pause(0.1);
                set(plotH, 'visible', 'on'); 
                pause(0.05);
            end
            function done
                set(plotH, 'visible', vi);
                set(plotH, 'MarkerSize', ms);
            end
        end
        
        
        function [HL, fcnMotion, btns, sortI, plots]=Legend(plotHs, names, ...
                idxsWithName, xMargin, yMargin, doMotion, freq, priorFcn,...
                doubleClicker, doJavaLegend, oldJavaBtns, southComponent,...
                selectedCallback)
            
            if nargin<13
                selectedCallback=[];
                if nargin<12
                    southComponent=[];
                    if nargin<11
                        oldJavaBtns=[];
                        if nargin<10
                            if nargin<9
                                doubleClicker=[];
                                if nargin<8
                                    priorFcn=[];
                                    if nargin<7
                                        freq=[];
                                        if nargin<6
                                            doMotion=true;
                                            if nargin<5
                                                yMargin=[];
                                                if nargin<4
                                                    xMargin=[];
                                                    if nargin<3
                                                        idxsWithName=[];
                                                        if nargin<2
                                                            names={};
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            doJavaLegend=false;
                        end
                    end
                end
            end
            app=BasicMap.Global;
            if doJavaLegend
                sup1=app.supStart;
                sup2=app.supEnd;
            end
            btns=[];
            if isa(plotHs, 'Plots')
                that=plotHs;
                plots=plotHs;
                N=plots.N;
            else
                that=[];
                plots=Plots;
                plots.setNames(names);
                if iscell(plotHs)
                    N=length(plotHs);
                    plotHs2=zeros(1,N);
                    for i=1:N
                        plotHs2(i)=plotHs{i};
                    end
                    plotHs=plotHs2;
                else
                    [R,N]=size(plotHs);
                    if R>1
                        otherPlotHs=plotHs(2:end,:);
                        plotHs=plotHs(1,:);
                        plots.setOtherHs(otherPlotHs);
                    end
                end
                plots.setHs(plotHs);
                names=plots.names;
            end
            if isempty(doJavaLegend)
                HL=[];
                fcnMotion=[];
                sortI=[];
                return;
            end
            if N<1
                return;
            end              
            if isempty(idxsWithName)
                idxsWithName=1:N;
            end
            assert(N>=length(idxsWithName), ...
                '# of plotHs must be >= # of idxsWithName');
            ax=get(plots.Hs(1), 'Parent');
            fig=get(ax, 'Parent');
            curOver=0;
            overH=[];
            try
                set(fig,'WindowButtonMotionFcn', @legMotion);
            catch ex
            end
            fcnMotion=@legMotion;
            X=xlim(ax);
            Y=ylim(ax);
            assert( length(names) == length(idxsWithName), ...
                'Must be as many namedPlotIdxs as names');

            mx=plots.mx;
            cnts=plots.cnts;
            if isempty(freq)
                Cnt=sum(cnts);
            else
                if length(freq)==length(cnts)
                    if any(isnan(freq))
                        freq(isnan(freq))=-1;
                    end
                    Cnt=sum(freq(freq>=0));
                    cnts=freq;
                    mx=max(cnts);
                else
                    Cnt=freq;
                end
            end
            nNames=length(names);
            cnts2=zeros(1, nNames);
            cntPref=app.getNumeric('cntPref', 2);
            if ~doJavaLegend
                legendHs=zeros(1, nNames);
                for i=1:nNames
                    idx=idxsWithName(i);
                    cnt=cnts(idx);
                    cnts2(i)=cnt;
                    if cnt>=0
                        strCnt=['    ^{\color{blue}' ...
                            String.encodePercent(cnt, Cnt, 1) '  '...
                            '\color[rgb]{.2 .5 .2} ' ...
                            String.encodeCount(cnt, cntPref) '} '];
                        mkSz=6 + (cnt/mx*15);
                        clr=get(plots.Hs(idx), 'MarkerEdgeColor');
                    
                        legendHs(i)=plot(ax, X(1)-X(2), Y(1)-Y(2), 's', ...
                            'MarkerSize', mkSz, 'Color',  clr, ...
                            'MarkerFaceColor', clr, 'LineStyle', 'none');
                    else
                        if isa(plotHs, 'Plots')
                            clr=get(plots.Hs(idx), 'MarkerEdgeColor');
                            legendHs(i)=plot(ax, X(1)-X(2), Y(1)-Y(2), 's', ...
                                'MarkerSize', 6, 'Color',  clr, ...
                                'MarkerFaceColor', clr, 'LineStyle', 'none');
                        else
                            legendHs(i)=plotHs(i);
                        end
                        names{i}=names{i};
                    end
                    
                end
                xlim(ax,X);
                ylim(ax,Y);
                [~,sortI]=sort(cnts2, 'descend');
                HL=legend(ax, legendHs(sortI), names(sortI), ...
                    'Location', 'northeast','AutoUpdate', 'off');
                if ~isempty(xMargin) && ~isempty(yMargin)
                    p=get(HL, 'position');
                    p2=Gui.GetNormalized(ax);
                    xNudge=1-(p2(1)+p2(3));
                    xNudge=xNudge-xMargin;
                    yNudge=1-(p2(2)+p2(4));
                    yNudge=yNudge-yMargin;
                    set(HL, 'position', [p(1)+xNudge, p(2)+yNudge p(3) p(4)]);
                end
                plots.setLegendHs(legendHs);
                HL.ItemHitFcn=@(h,event)Plots.LegendClick(plots, event,...
                    idxsWithName, cnts2);
            else
                for i=1:nNames
                    idx=idxsWithName(i);
                    cnt=cnts(idx);
                    cnts2(i)=cnt;
                    clr=get(plots.Hs(idx), 'MarkerEdgeColor');
                    if cnt>=0
                        strSymbol=['<font  ' Gui.HtmlHexColor(clr)...
                            '>&bull;</font>'];
                        f=11+(cnt/mx*15);
                        f= ceil((f-4)/4)+3;
                        if f>=10
                            f=10;
                        end
                        f=String.encodeInteger(f);
                        str=['<font size="' f '">' strSymbol '</font>'];
                        strCnt=[sup1 ' <font color="blue">' ...
                            String.encodePercent(cnt,Cnt,1) '</font>'...
                            '  <font ' Html.Color([.2 .5 .2]) ' ><i>' ...
                            String.encodeCount(cnt, cntPref) '</i></font>' ...
                            sup2 ];
                    else
                        ls=get(plots.Hs(idx), 'LineStyle');
                        if strcmpi(ls, 'none')
                            str=['<font size="7" ' Gui.HtmlHexColor(clr)...
                                '>&bull;</font>&nbsp;'];
                        else
                            str=['<font size="8" ' Gui.HtmlHexColor(clr)...
                                '>&minus;</font>&nbsp;'];
                        end
                        strCnt='';
                    end
                    names{i}=['<html>' str names{i} '&nbsp;&nbsp;'...
                        strCnt Html.EncodeSort('name', strtrim(char(...
                        edu.stanford.facs.swing.Basics.RemoveXml(...
                        lower(names{i}))))) ...
                        Html.EncodeSort('frequency', cnt) ...
                        '</html>'];
                end
                [~, sortI]=sort(cnts2, 'descend');
                [outerPnl, ~, btns,sortGui]=...
                    Radio.Panel3(names(sortI), 1:nNames, 13, ...
                    @(h,e, i)Plots.JavaLegendClick(plots, i, sortI,...
                    idxsWithName, plots.Hs, cnts2, h, selectedCallback), ...
                    true,Html.Wrap(['<i>Deselecting '...
                    '<b>hides</b> items in plot & selecting '...
                    '<b>unhides</b></i>']));
                sortGui.setProperties(app, 'javaLegend.');
                if ~isempty(oldJavaBtns)
                    it=oldJavaBtns.iterator;
                    while it.hasNext
                        btn=it.next;
                        if ~btn.isSelected
                            k=StringArray.IndexOf(names,char(btn.getText));
                            m=find(sortI==k, 1);
                            if ~isempty(m)
                                b=btns.get(m-1);
                                b.setSelected(false);
                                set(plots.Hs(k), 'visible', 'off');
                                disp('was OFF');
                            end
                        end
                    end
                    sortGui.setAllChbText;
                elseif isequal(idxsWithName, 1:nNames)
                    for i=1:nNames
                        if isequal('off', get(plots.Hs(i), 'visible'))
                            b=btns.get(i-1);
                            b.setSelected(false);
                        end
                    end
                end
                if ~isempty(southComponent)
                    bp=Gui.BorderPanel(0,0);
                    bp.add(outerPnl, 'Center');
                    bp.add(southComponent, 'South');
                    outerPnl=bp;
                end
                if ~isempty(oldJavaBtns) && oldJavaBtns.size>0
                    HL=Gui.WindowAncestor(oldJavaBtns.get(0));
                    HL.setContentPane(outerPnl);
                    HL.pack;
                else
                    pu=showMsg(outerPnl, 'Legend', 'north east+', false, ...
                        false, 0, false);
                    HL=pu.dlg;
                end
                sortGui.dlg=HL;
                btns.get(0).getParent.scrollRectToVisible(java.awt.Rectangle(0,0,1,1))
                if ismac
                    if 1<size(get(0, 'MonitorPositions'), 1)
                        if isequal('on', get(fig, 'visible'))
                            MatBasics.RunLater(@(h,e)relocate(getjframe(fig)), .52)
                        else
                            disp('huh');
                        end
                    end
                end
            end
            if ~isempty(that)
                that.legendH=HL;
            end
            function relocate(ref)
                Gui.LocateJava(HL, ref, 'north east+');
                HL.setVisible(true);
            end
            function legMotion(hObj, event)
                if ~isvalid(ax)
                    return;
                end
                cp=get(ax, 'CurrentPoint');
                if isempty(cp)
                    return;
                end
                C=size(cp,2);
                x=cp(2,1);
                y=cp(2,2);
                z=cp(2,3);
                try
                    e=Gui.GetPixels(HL);
                catch ex
                    e=[];
                end
                %[normX normY normFigX normFigY e]
                if ~isempty(e)
                    cp2=get(get(HL, 'Parent'), 'currentpoint');
                    if cp2(2)<=e(2)+e(4) && cp2(2)>=e(2)
                        if cp2(1)>=e(1) && cp2(1)<= e(1)+e(3)
                            if ~isempty(doubleClicker)
                                doubleClicker.stopListening;
                            end
                            return;
                        end
                    end
                end
                if ~isempty(doubleClicker)
                    doubleClicker.startListening;
                end
                if ~doMotion
                    if ~isempty(priorFcn)
                        feval(priorFcn, hObj, event);
                    end
                    return;
                end
                plots.initStats;
                if plots.is3D
                    [D,II]=pdist2(plots.mdns, [x y z], 'euclidean', 'Smallest', 1);
                    limit=pdist2(plots.mdns(II,:), ...
                        plots.mdns(II,:)+(plots.stds(II,:)*plots.distFactor));
                else
                    [D,II]=pdist2(plots.mdns(:,1:2), [x y], 'euclidean', 'Smallest', 1);
                    limit=pdist2(plots.mdns(II,1:2), plots.mdns(II,1:2)+...
                        (plots.stds(II,1:2)*plots.distFactor));
                end
                if Plots.DEBUGGING
                    %is3D
                    [x y z ; II D limit]
                    plots.names{II}
                    %cp
                    %[plots.mdns(II,:);plots.stds(II,:)]
                    if II==plots.mxI
                        disp('hey');
                    end
                end
                if D<=limit
                    over=II;
                else
                    over=0;
                end
                if ~ishandle(overH)
                    return;
                end
                if curOver ~= over
                    curOver=over;
                    if over>0
                        %disp(['Over ' names{II}]);                        
                        nameIdx=find(idxsWithName==II,1);
                        if nameIdx>0
                            nm=names{nameIdx};
                            if doJavaLegend
                                nm=strrep(nm, '<sup>', '_{{');
                                nm=strrep(nm, '</sup>', '}');
                                nm=char(...
                                    edu.stanford.facs.swing.Basics.RemoveXml(nm));
                                nm=strrep(nm, '_{{', '^{');
                                nm=strrep(nm, '&bull;','');                                
                            end
                            if isempty(overH)
                                if plots.is3D
                                    overH=text(ax, x, y, z, nm, ...
                                        'fontSize', 9, ...
                                        'color', [0 0 .5],...
                                        'EdgeColor', 'red', ...
                                        'FontName', 'Arial', ...
                                        'backgroundColor', [255/256 252/256 170/256]);
                                else
                                    overH=text(ax, x, y, nm, ...
                                        'fontSize', 9, 'color', [0 0 .5], 'EdgeColor',...
                                        'red', 'FontName', 'Arial', ...
                                        'backgroundColor', [255/256 252/256 170/256]);
                                end
                            else
                                if plots.is3D
                                    set(overH, 'visible', 'on', 'Position', ...
                                        [x y z], 'String', nm);
                                    %uistack(overH, 'top');
                                else
                                    set(overH, 'visible', 'on', 'Position', ...
                                        [x y 0], 'String', nm);
                                end
                            end
                        else
                            disp('background');
                        end
                    else
                        %disp('Over NOTHING');
                        if ~isempty(overH)
                            set(overH, 'visible', 'off');
                        end
                    end
                end
                if ~isempty(priorFcn)
                    feval(priorFcn, hObj, event);
                end
            end
        end
       
        function LegendClick(this, event, idxsWithName, nums)
            H=get(event.Peer, 'UserData');
            idx_=find(this.legendHs==event.Peer, 1);
            idx=idxsWithName(idx_);
            if ~isempty(H)
                set(event.Peer, 'UserData', []);
                delete(H);
                this.toggleVisibility(idx);
                return;
            end
            this.initStats;
            plotH=this.Hs(idx);
            perc= nums(idx_)/max(nums);
            Plots.Flash(plotH, perc, 4);
            X=this.mdns(idx, 1);
            Y=this.mdns(idx, 2);
            ax=get(this.Hs(idx), 'Parent');
            str=get(event.Peer, 'DisplayName');
            fsz=11;
            if this.is3D
                Z=this.mdns(idx, 3);
                Plots.ShowLegendTip(ax, event.Peer, X, Y, str, fsz, 0, Z);
                [X Y Z]
            else
                Plots.ShowLegendTip(ax, event.Peer, X, Y, str, fsz, 0);
            end
            this.toggleVisibility(idx);   
            this.flashOtherPlots(idx);            
        end
        
        function ok=JavaLegendClick(this, idx_, sortI, idxsWithName, plotHs, ...
                nums, h, selectedCallback)
            idx=idxsWithName(sortI(idx_));
            plotH=plotHs(idx);
            if ~ishandle(plotH)
                msg('Supervisors not showing', 5, 'east+');
                ok=false;
                return;
            end
            ok=true;
            hasOtherMap=this.isInOtherPlotMap(idx);
            doingAll=nargin==8 && isequal(Plots.DOING_ALL, ...
                    char(h.getActionCommand));
            if ~h.isSelected
                set(plotH, 'visible', 'off');
                times=3;
            else
                if doingAll && hasOtherMap
                    %nope
                    if this.otherPlotMap.size<length(this.Hs)
                        h.setSelected(false)
                        return;
                    end
                end
                set(plotH, 'visible', 'on');
                times=4;
            end
            if doingAll
                if ~isempty(selectedCallback)
                    feval(selectedCallback, h, idx, true);
                end
                this.setOtherPlotVisible(idx, h.isSelected);
                return;
            end
            perc=nums(idx)/max(nums);
            Plots.Flash(plotH, perc, times);
            if ~isempty(selectedCallback)
                feval(selectedCallback, h, idx, false);
            end
            this.flashOtherPlots(idx, h.isSelected);            
        end
        
        function ShowLegendTip(ax, hLegend, posX, posY, str, fsz, secs, posZ)
            if nargin<8
                [normX, normY]=Gui.DataToAxNorm(ax, posX, posY);
                tipH=text(normX+.03, normY+.03, str, ...
                    'fontSize', fsz, 'color', [0 0 .5], 'EdgeColor',...
                    'red', 'FontName', 'Arial', 'parent', ax, ...
                    'FontAngle', 'italic', 'units', 'normalized',...
                    'backgroundColor', [255/256 252/256 170/256]);
            else
                tipH=text(ax, posX, posY, posZ, str, ...
                    'fontSize', fsz, 'color', [0 0 .5], 'EdgeColor',...
                    'red', 'FontName', 'Arial', 'FontAngle', 'italic',...
                    'backgroundColor', [255/256 252/256 170/256]);
                
            end
            set(tipH, 'ButtonDownFcn', @(h,e)freeze(h, hLegend));
            if secs>0
                MatBasics.RunLater(@(h,e)closeTip, secs);
            else
                freeze(tipH, hLegend);
            end
            
            function freeze(hText, hLegend)
                if isempty(get(tipH, 'UserData'))
                    set(hText, 'UserData', true);
                    set(hText, 'FontSize', fsz-2);
                    set(hText, 'FontAngle', 'normal');
                    set(hLegend, 'UserData', hText);
                else
                    set(hLegend, 'UserData', []);
                    delete(hText);
                end
            end
            
            function closeTip
                if isempty(get(tipH, 'UserData'))
                    delete(tipH);
                end
            end
        end
        
    end
end