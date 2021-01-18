%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef SortGui < handle
    properties(Constant)
        PROP_ORDER='sortGui.sortOrder';
        PROP_ASCENDING='sortGui.sortAscending';
        PROP_SEARCH='sortGui.search';
        PROP_RESTRICT='sortGui.restrict';
    end
    
    properties
        fncRefresh;
        dlg;
    end
    
    properties(SetAccess=private)
        givenTip=false;
        path;
        strings;
        allChb;
        allMsg;
        sortIdxs;
        visibleIdxs;
        sortTypes;
        sortKeys;
        btnUp;
        btnDown;
        btnSortDir;
        btnSearch;
        options;
        sortAscending=true;
        nOptions;
        checkBoxes
        innerPnl;
        cb;
        map;
        propPrefix;
        allChbPnl;
        searchPnl;
        searchJtf;
        searching=false;
        found;
        foundLabel;
        foundIdx;
        btnRestrict;
        originalSize;
        originalHeight;
        firstCb;
        sortLabel;
    end
    
    methods
        function this=SortGui(dlg, allChb, allMsg, allChbPnl, options, ...
                checkBoxes, innerPnl)
            pnl=Gui.Panel;
            this.dlg=dlg;
            this.allChb=allChb;
            this.allMsg=allMsg;
            this.options=options;
            this.nOptions=length(options);
            this.sortIdxs=1:this.nOptions;
            this.checkBoxes=checkBoxes;
            this.innerPnl=innerPnl;
            [this.sortKeys, this.sortTypes]=Html.DecodeSort(options{1});
            this.sortLabel=Gui.Label('Sort:');
            if ~isempty(this.sortKeys)
                pnl.add(this.sortLabel);
                fl=javaObjectEDT('java.awt.FlowLayout', 0,0,0);
                pnl.setLayout(fl);
                labels=['<html><i>sort by</i></html>' this.sortKeys];
                cb=Gui.Combo(labels, 1,[],[],@(h,e)sortItems(this, h), ...
                    'Change order of list');
                pnl.add(cb);
                this.cb=cb;
                allChbPnl.add(Gui.Label('  '), 'Center');
            end
            this.path=BasicMap.Path;
            this.btnSortDir=Gui.ImageButton(...
                fullfile(this.path, 'sortDown.png'), ...
                'List order is ascending, click for DESCENDING?', ...
                @(h,e)toggleAscending(this));
            pnl.add(this.btnSortDir);
            this.btnSearch=Gui.ImageButton(...
                fullfile(this.path, 'find16.gif'), ...
                'Search items in this list', ...
                @(h,e)toggleSearch(this));
            pnl.add(this.btnSearch);
            this.initSearch;
            allChbPnl.add(pnl, 'East');
            this.allChbPnl=allChbPnl;
            this.seeSortCues;
        end
        
        function setSortDir(this)
            if this.sortAscending
                this.btnSortDir.setIcon(Gui.Icon('sortDown.png'));
                this.btnSortDir.setToolTipText(...
                    'Order is ascending, click for DESCENDING.');
            else
                this.btnSortDir.setIcon(Gui.Icon('sortUp.png'));
                this.btnSortDir.setToolTipText(...
                    'Order is descending, click for ASCENDING.');
            end
        end
        
        function toggleAscending(this)
            this.sortIdxs=flip(this.sortIdxs);
            this.refreshSearchAndSort(true);
            this.seeSortCues; 
        end
        
        function seeSortCues(this)
            sorting=~isempty(this.cb) && this.cb.getSelectedIndex>0;
            this.btnSortDir.setVisible(sorting);
            this.sortLabel.setVisible(sorting);
        end
        
        function refreshSearchAndSort(this, sorting)
            if ~isempty(this.fncRefresh)
                if this.btnRestrict.isSelected
                    pnl=feval(this.fncRefresh);
                else
                    was=this.visibleIdxs;
                    this.visibleIdxs=[];
                    pnl=feval(this.fncRefresh);
                    this.visibleIdxs=was;
                end
                if ~isempty(this.dlg)
                    if isempty(this.originalSize)
                        sz=this.dlg.getSize;
                        if sz.height>0 && sz.width>0
                            this.originalSize=sz;
                            this.originalHeight=this.originalSize.height*.8;
                        end
                    end
                    this.dlg.pack;
                    if ~isempty(this.originalSize)
                        sz=this.dlg.getSize;
                        if ~sz.equals(this.originalSize)
                            this.dlg.setSize(this.originalSize);
                        end
                    end
                end
            else
                this.innerPnl.removeAll;
                if isempty(this.visibleIdxs) || ~this.btnRestrict.isSelected
                    for jj=1:this.nOptions
                        idx=this.sortIdxs(jj);
                        this.innerPnl.add(this.checkBoxes.get(idx-1));
                    end
                else
                    for jj=1:this.nOptions
                        idx=this.sortIdxs(jj);
                        if isempty(find(this.visibleIdxs==idx,1))
                            continue;
                        end
                        this.innerPnl.add(this.checkBoxes.get(idx-1));
                    end
                end
                if nargin>1 && sorting && ~isempty(this.dlg)
                    this.dlg.pack;
                else
                    this.innerPnl.repaint;
                end
            end
            this.sortAscending=~this.sortAscending;
            this.setAscending;
            this.setSortDir;
            this.handleFound;
            this.seeFirstCb;
            drawnow;
        end
        
        function setAscending(this)
            if ~isempty(this.map)
                if this.sortAscending
                    this.map.set(this.propAscending, 'true');
                else
                    this.map.set(this.propAscending, 'false');
                end
            end
        end
        
        function prop=propOrder(this)
            prop=[this.propPrefix SortGui.PROP_ORDER];
        end
        
        function prop=propAscending(this)
            prop=[this.propPrefix SortGui.PROP_ASCENDING];
        end
        
        function prop=propSearch(this)
            prop=[this.propPrefix SortGui.PROP_SEARCH];
        end
        
        function prop=propRestrict(this)
            prop=[this.propPrefix SortGui.PROP_RESTRICT];
        end

        function setProperties(this, map, propPrefix, searchDflt)
            if nargin<4
                searchDflt=false;
            end
            this.map=map;
            this.propPrefix=propPrefix;
            if map.has(this.propOrder)
                idx=map.getNumeric(this.propOrder, 0);
                if isnumeric(idx) && idx>0
                    this.cb.setSelectedIndex(idx-1);
                    %this.sort(idx);
                end
            end
            if ~map.is(this.propAscending, true)
                this.toggleAscending;
                this.setSortDir;
            end
            if map.is(this.propSearch, searchDflt)
                MatBasics.RunLater(@(h,e)toggleSearch(this), .12);
            end
            if map.is(this.propRestrict, true)
                this.btnRestrict.setSelected(true);
                this.restrict;
            end
        end
        
        function sortItems(this, h)
            idx=h.getSelectedIndex+1;
            if ~isempty(this.map)
                this.map.set(this.propOrder, idx);
            end
            this.sort(idx);
        end
        
        function sort(this, idx)
            if idx<=1
                this.sortIdxs=1:this.nOptions;
            else
                isNum=isequal('N', this.sortTypes{idx-1});
                if isNum
                    rm=realmin;
                    values=zeros(1, this.nOptions);
                else
                    values=cell(1, this.nOptions);
                end
                key=this.sortKeys{idx-1};
                for jj=1:this.nOptions
                    v=Html.DecodeSortValues(this.options{jj}, key);
                    if isNum
                        if isempty(v)
                            values(jj)=rm;
                        else
                            values(jj)=str2double(v{1});
                        end
                    else
                        if isempty(v)
                            values{jj}='';
                        else
                            values{jj}=v{1};
                        end
                    end
                end
                [~, this.sortIdxs]=sort(values);
            end
            if this.sortAscending
                this.sortIdxs=flip(this.sortIdxs);
            end
            this.sortAscending=~this.sortAscending;
            this.toggleAscending;
            this.setAscending;
        end
        
        function handleFound(this)
            if ~isempty(this.foundLabel)
                if isempty(this.found)
                    this.btnDown.setVisible(false);
                    this.btnUp.setVisible(false);
                    if this.searchJtf.getText.length==0
                        this.foundLabel.setText('');
                    else
                        this.foundLabel.setText(...
                            ['<html><b><font color="red">0</font></b>/' ...
                            num2str(this.nOptions) '</html>']);
                    end
                else
                    N=length(this.found);
                    if ~this.btnRestrict.isSelected
                        s=[num2str(this.foundIdx) '/' num2str(N)];
                        this.btnDown.setVisible(true);
                        this.btnUp.setVisible(true);
                    else
                        s=[num2str(N) '/' num2str(this.nOptions)];
                        this.btnDown.setVisible(false);
                        this.btnUp.setVisible(false);
                    end
                    this.foundLabel.setText(s);
                    if ispc
                        this.dlg.pack;
                    end
                end
            end
        end
        
        function initSearch(this)
            jtf=javaObjectEDT('javax.swing.JTextField');
            jtf.setColumns(15)
            %jtf.setHorizontalAlignment(jtf.RIGHT);
            jj=handle(jtf, 'CallbackProperties');
            set(jj, 'KeyTypedCallback', @(h,e)search(this, h));
            this.searchPnl=Gui.BorderPanel;
            pnl=Gui.Panel;
            fl=javaObjectEDT('java.awt.FlowLayout', 0,0,0);
            pnl.setLayout(fl);
            emptyBorder=javax.swing.BorderFactory.createEmptyBorder(0,0,0,0);
            pnl.setBorder(emptyBorder);
            pnl.add(javaObjectEDT('javax.swing.JLabel', ' '));
            this.btnRestrict=Gui.CheckBox(Html.Wrap(Html.Img(...
                'hide.png')), false, [],...
                [], @(h,e)restrict(this), ...
                Html.Wrap(['Select ' Html.Img('magnify.png') ...
                ' to restrict list to found items']));
            pnl.add(this.btnRestrict);

            this.btnUp=Gui.ImageButton(...
                fullfile(this.path, 'upArrow.png'), ...
                'Find next', ...
                @(h,e)nextSearch(this, -1));
            pnl.add(this.btnUp);
            this.btnDown=Gui.ImageButton(...
                fullfile(this.path, 'downArrow.png'), ...
                'Find previous', ...
                @(h,e)nextSearch(this, 1));
            pnl.add(this.btnDown);
            this.foundLabel=javaObjectEDT('javax.swing.JLabel');
            pnl.add(this.foundLabel);
            pnl.add(jtf);
            btnClear=Gui.ImageButton(...
                fullfile(BasicMap.Path, 'cancel.gif'), ...
                'Clear/close search', @(h,e)clearSearch(this));
            pnl.add(btnClear);
            pnl.add(javaObjectEDT('javax.swing.JLabel', '     '));
            this.searchPnl.add(pnl, 'East');
            this.searchJtf=jtf;
            this.btnDown.setPreferredSize(java.awt.Dimension(12, 16))
            this.btnUp.setPreferredSize(java.awt.Dimension(12, 16))
        end
        
        function restrict(this)
            if this.btnRestrict.isSelected
                this.btnUp.setVisible(false);
                this.btnDown.setVisible(false); 
                if ~isempty(this.map)
                    this.map.set(this.propRestrict, 'true');
                end
            else
                this.btnUp.setVisible(true);
                this.btnDown.setVisible(true);
                if ~isempty(this.map)
                    this.map.set(this.propRestrict, 'false');
                end
            end
            this.sortAscending=~this.sortAscending;
            this.refreshSearchAndSort;
        end
        
        function toggleSearch(this)
            if this.searching
                this.allChbPnl.remove(this.searchPnl);
                if ~isempty(this.map)
                    this.map.set(this.propSearch, 'false');
                end
            else
                sz=this.innerPnl.getSize;
                if sz.width==0
                    sz=this.checkBoxes.get(0).getParent.getSize;
                end
                if sz.width>0 && sz.width<500
                    where='North';
                    if sz.width<200
                        this.searchJtf.setColumns(9);
                    end
                else
                    where='Center';
                end
                
                this.allChbPnl.add(this.searchPnl, where);
                this.searchJtf.requestFocus;
                if ~isempty(this.map)
                    this.map.set(this.propSearch, 'true');
                end
                
            end
            this.allChbPnl.repaint;
            this.searching=~this.searching;
            this.dlg.pack;
        end
        
        function nextSearch(this, i)
            N=length(this.found);
            this.foundIdx=this.foundIdx+i;
            if this.foundIdx>N
                this.foundIdx=1;
            elseif this.foundIdx<1
                this.foundIdx=N;
            end
            cb_=this.found{this.foundIdx};
            b=cb_.getBounds;
            b.height=b.height*2.5;
            cb_.getParent.scrollRectToVisible(b);
            this.handleFound;
        end
        
        function clearSearch(this)
            if this.searchJtf.getText.length==0
                this.toggleSearch;
            else
                this.searchJtf.setText('');
                this.found={};
                this.search(this.searchJtf);
                this.handleFound;
            end
        end
        
        function search(this, h)
            this.found={};
            this.visibleIdxs=[];
            N=this.nOptions-1;
            s=lower(char(h.getText));
            N=this.nOptions-1;
            if isempty(this.strings)
                this.strings={};
                for i=0:N
                    this.strings{end+1}=string(lower(...
                        char(edu.stanford.facs.swing.Basics.RemoveXml(...
                        this.checkBoxes.get(i).getText))));
                end
            end
            this.firstCb=this.checkBoxes.get(0);
            if isempty(s)
                for i=0:N
                    cb_=this.checkBoxes.get(i);
                    cb_.setEnabled(true);
                end
            else
                first=false;
                N=this.nOptions;
                for i=1:N
                    idx=this.sortIdxs(i);
                    cb_=this.checkBoxes.get(idx-1);
                    if this.strings{idx}.contains(s)
                        cb_.setEnabled(true);
                        if ~first
                            this.firstCb=cb_;
                            first=true;
                        end
                        this.found{end+1}=cb_;
                        this.visibleIdxs(end+1)=idx;
                    else
                        cb_.setEnabled(false);
                    end
                end
            end
            this.foundIdx=1;
            if this.btnRestrict.isSelected
                this.sortAscending=~this.sortAscending;
                this.refreshSearchAndSort;
            else
                this.handleFound;
                this.seeFirstCb;
                if ~this.givenTip && h.getText.length==1
                    this.givenTip=true;
                    MatBasics.RunLater(@(h,e)tip, .51);
                end
            end
            %this.dlg.pack;
            function tip
                try
                    BasicMap.Global.showToolTip(this.btnRestrict,...
                        [], 2, -50);
                catch ex
                    disp(ex);
                end
            end
        end
        
        function seeFirstCb(this)
            cb_=[];
            if ~isempty(this.firstCb)
                cb_=this.firstCb;
            else
                for jj=0:this.nOptions-1
                    if this.checkBoxes.get(jj).isSelected
                        cb_=this.checkBoxes.get(jj);
                    end
                end
            end
            if ~isempty(cb_)
                b=cb_.getBounds;
                b.height=b.height*2.5;
                cb_.getParent.scrollRectToVisible(b);
            end
        end
        
        function [idxs_, N]=setAllChbText(this)
            [idxs_, N]=Gui.SetAllChb(this.allChb, this.allMsg, this.checkBoxes);
        end
        
        function idxs=reorder(this, idxs)
            idxs=MatBasics.Sort1D(idxs, this.sortIdxs);
        end
    end
    
    methods(Static)
        function items=EncodeLabelsAndList(labels, list, gtp)
            if iscell(list)
                N=length(list);
            else
                N=list.size;
            end
            app=BasicMap.Global;
            items=cell(1, N);
            for i=1:N
                if iscell(list)
                    item=list{i};
                else
                    item=list.get(i-1);
                end
                frequency=sum(labels==i);
                cue1=Html.EncodeSort('frequency', frequency);
                cue2=Html.EncodeSort('name', item);
                if nargin>2
                    frequency=gtp.encodeCount(frequency) ;
                else
                    frequency=String.encodeCount(frequency);
                end
                    
                items{i}=['<html>' item app.smallStart...
                        '<font ' Html.HexColor([.33 .33 .52]) '>'...
                        '<b> &lt;' frequency ...
                        ' events&gt;</b></font>'  app.smallEnd ...
                        cue1 cue2 '</html>'];
            end
        end
    end
end


