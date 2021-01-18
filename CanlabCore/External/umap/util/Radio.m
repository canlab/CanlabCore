%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef Radio
    methods(Static)
        function [jp, bg]=Panel(varargin)
            [~, jsa, ~, defaultIdx]=getMsgTypeAndOptions(...
                javax.swing.JOptionPane.QUESTION_MESSAGE, varargin);
            [jp,bg]=Radio.Panel2(jsa, defaultIdx, 11);
        end
        
        function [jp, bg]=PanelAndCallback(callback, noDefault, varargin)
            [~, jsa, ~, defaultIdx]=getMsgTypeAndOptions(...
                javax.swing.JOptionPane.QUESTION_MESSAGE, varargin);
            if noDefault
                defaultIdx=0;
            end
            [jp,bg]=Radio.Panel2(jsa, defaultIdx, 11, false, callback);
        end

        function [outerPnl, bg, checkBoxes, innerPnl]=Panel2(javaArray, ...
                dfltIdxs, addScrollAfter, checkBox, callback, ...
                priorCmps, sortIdxs, visibleIdxs)
            if nargin<4 || ~checkBox
                checkBox=false;
                jObj='javax.swing.JRadioButton';
            else
                jObj='javax.swing.JCheckBox';
            end
            N=javaArray.length;
            isJPanels=N>0 && isa(javaArray(1),'javax.swing.JPanel');
            innerPnl=javaObjectEDT('javax.swing.JPanel');
            innerPnl.setLayout(javax.swing.BoxLayout(innerPnl, 1));
            if checkBox
                bg=javaObjectEDT('javax.swing.JPanel');
                if ~isJPanels
                    columns=1;
                    horizontal=true;
                    for i=1:N
                        if javaArray(i).length>2
                            horizontal=false;
                            break;
                        else
                        end
                    end
                else
                    columns=1;
                    if N>1
                        columns=javaArray(1).getComponentCount;
                    end
                    horizontal=false;
                end
                if ~horizontal
                    if ~isJPanels
                        innerPnl.setLayout(javaObjectEDT(...
                            'java.awt.GridLayout', N, 1));
                    else
                        innerPnl.setLayout(javaObjectEDT(...
                            'java.awt.GridBagLayout'));
                        gbc=javaObjectEDT('java.awt.GridBagConstraints');
                        gbc.fill=2;
                    end
                else
                    innerPnl.setLayout(javaObjectEDT(...
                        'java.awt.GridLayout', 1, N));
                end
            else
                bg=javax.swing.ButtonGroup;
            end
            checkBoxes=java.util.ArrayList;
            idx=1;
            scrollToHeight=0;
            first=true;
            if ~isJPanels
                %javaArray is char or strings
                for i=1:N
                    rb=javaObjectEDT(jObj, javaArray(i));
                    rObj=rb;
                    pb=handle(rb, 'CallbackProperties');
                    if nargin>=5
                        set(pb,'ActionPerformedCallback', callback);
                    end
                    checkBoxes.add(rb);
                    bg.add(rb);
                    innerPnl.add(rObj);
                    if ~isempty(find(dfltIdxs==idx, 1))
                        rb.setSelected(true);
                        if first
                            d=rb.getPreferredSize;
                            scrollToHeight=d.height*(i-1);
                            first=false;
                        end
                    end
                    idx=idx+1;
                end                
            else
                try
                    if nargin>5
                        priorCnt=length(priorCmps);
                        columns=priorCnt/N;
                    end
                    if nargin>6
                        if nargin<8 || isempty(visibleIdxs)
                            for i=1:N
                                idx2=(2*sortIdxs(i))-1;
                                rb=priorCmps{idx2};
                                checkBoxes.add(rb);
                                gbc.gridy=i-1;
                                gbc.gridx=0;
                                gbc.weightx=.8;
                                innerPnl.add(rb, gbc);
                                for j=2:columns
                                    gbc.gridx=j-1;
                                    gbc.weightx=.2;
                                    jp=priorCmps{(idx2+j)-1};
                                    innerPnl.add(jp, gbc);
                                end
                                idx=idx+1;
                            end
                        else
                            for i=1:N
                                si=sortIdxs(i);
                                if isempty(find(visibleIdxs==si, 1))
                                    continue;
                                end
                                idx2=(2*si)-1;
                                rb=priorCmps{idx2};
                                checkBoxes.add(rb);
                                gbc.gridy=i-1;
                                gbc.gridx=0;
                                gbc.weightx=.8;
                                innerPnl.add(rb, gbc);
                                for j=2:columns
                                    gbc.gridx=j-1;
                                    gbc.weightx=.2;
                                    jp=priorCmps{(idx2+j)-1};
                                    innerPnl.add(jp, gbc);
                                end
                                idx=idx+1;
                            end
                        end
                    else
                    for i=1:N
                            rb=javaObjectEDT(jObj, javaArray(i).getComponent(0).getText);
                            pb=handle(rb, 'CallbackProperties');
                            if nargin>=5
                                set(pb,'ActionPerformedCallback', callback);
                            end
                            checkBoxes.add(rb);
                            gbc.gridy=i-1;
                            gbc.gridx=0;
                            gbc.weightx=.8;
                            innerPnl.add(rb, gbc);
                            for j=2:columns
                                gbc.gridx=j-1;
                                gbc.weightx=.2;
                                jp=javaObjectEDT('javax.swing.JPanel');
                                jp.add(javaArray(i).getComponent(j-1));
                                innerPnl.add(jp, gbc);
                            end
                            if ~isempty(find(dfltIdxs==idx, 1))
                                rb.setSelected(true);
                                if first
                                    d=rb.getPreferredSize;
                                    scrollToHeight=d.height*(i-1);
                                    first=false;
                                end
                            end
                            idx=idx+1;
                        end
                    end
                catch ex
                    disp(['If javaArray contains JPanels then '...
                        'Radio.Panel2' newline 'expects component 0'...
                        ' to be JLabel and any 2nd component!' newline]);
                    ex.getReport
                end
            end
            outerPnl=javaObjectEDT('javax.swing.JPanel');
            outerPnl.setBorder(...
                javax.swing.BorderFactory.createEmptyBorder(1,8,1,4));
            if addScrollAfter<N
                d=innerPnl.getPreferredSize;
                d.width=floor(d.width*1.1);
                d.height=floor( (addScrollAfter/N)*d.height);
                innerPnl.setBorder([]);
                scroll=javaObjectEDT('javax.swing.JScrollPane', innerPnl);
                if scrollToHeight>0
                    scroll.getViewport.setViewPosition(...
                        java.awt.Point(0, scrollToHeight));
                end
                scroll.setPreferredSize(d);
                outerPnl=scroll;
            else                
                cb=javax.swing.BorderFactory.createCompoundBorder(...
                    javax.swing.BorderFactory.createBevelBorder(...
                    javax.swing.border.BevelBorder.RAISED), ...
                    javax.swing.BorderFactory.createEmptyBorder(8, 8, 8, 8));
                innerPnl.setBorder(cb);
                outerPnl.add(innerPnl);
            end
        end
        
        function [outerPnl, bg, btns, sortGui]=Panel3(strs, dfltIdxs, ...
                addScrollAfter, callback, checkBox, advice)
            isRadio=nargin<5||~checkBox;
            bg=[];
            sortGui=[];
            if isRadio
                jObj='javax.swing.JRadioButton';
            else
                jObj='javax.swing.JCheckBox';
            end
            N=length(strs);
            innerPnl=javaObjectEDT('javax.swing.JPanel');
            innerPnl.setLayout(javax.swing.BoxLayout(innerPnl, 1));
            if isRadio
                bg=javax.swing.ButtonGroup;
            end
            btns=java.util.ArrayList;
            scrollToHeight=0;
            for i=1:N
                rb=javaObjectEDT(jObj, strs{i});
                rObj=rb;
                pb=handle(rb, 'CallbackProperties');
                set(pb,'ActionPerformedCallback', @(h,e)innerCallback(h,e,i));
                btns.add(rb);
                if isRadio
                    bg.add(rb);
                end
                innerPnl.add(rObj);                
            end
            first=true;
            if ~isempty(dfltIdxs)
                for i=1:N
                    if ~isempty(find(dfltIdxs==i, 1))
                        rb=btns.get(i-1);
                        rb.setSelected(true);
                        if first
                            d=rb.getPreferredSize;
                            scrollToHeight=d.height*(i-1);
                            first=false;
                        end
                    end
                end
            end
            outerPnl=javaObjectEDT('javax.swing.JPanel');
            outerPnl.setBorder(...
                javax.swing.BorderFactory.createEmptyBorder(1,8,1,4));
            if addScrollAfter<N
                d=innerPnl.getPreferredSize;
                d.width=floor(d.width*1.1);
                d.height=floor( (addScrollAfter/N)*d.height);
                innerPnl.setBorder([]);
                scroll=javaObjectEDT('javax.swing.JScrollPane', innerPnl);
                if scrollToHeight>0
                    scroll.getViewport.setViewPosition(...
                        java.awt.Point(0, scrollToHeight));
                end
                scroll.setPreferredSize(d);
                outerPnl=scroll;
            else                
                cb=javax.swing.BorderFactory.createCompoundBorder(...
                    javax.swing.BorderFactory.createBevelBorder(...
                    javax.swing.border.BevelBorder.RAISED), ...
                    javax.swing.BorderFactory.createEmptyBorder(8, 8, 8, 8));
                innerPnl.setBorder(cb);
                outerPnl.add(innerPnl);
            end
            if ~isRadio
                bp=Gui.BorderPanel;                
                north=Gui.BorderPanel;
                bp.add(north, 'North');
                bp.add(outerPnl, 'Center');
                if nargin==6
                    south=Gui.Panel;
                    bp.add(south, 'South');
                    south.add(Gui.Label(advice));
                end
                allChb=Gui.CheckBox('All', true, [],[], @(h,e)doAll(h,e));
                north.add(allChb, 'West');
                outerPnl=bp;
                sortGui=SortGui([], allChb, 'All', north, strs, btns, innerPnl);
                sortGui.setAllChbText;
            end
            
            
            function innerCallback(h,e,i)
                if ~isempty(sortGui)
                    sortGui.setAllChbText;
                end
                feval(callback, h, e, i);
            end
            
            
            function doAll(h, e)
                state=h.isSelected;
                if state
                    Gui.ShowBusy(sortGui.dlg, 'Showing all', ...
                        'demoicon.gif', 7);
                else
                    Gui.ShowBusy(sortGui.dlg, 'Hiding all', 'hide.png', 5);
                end
                it=btns.iterator;
                idx=1;
                while it.hasNext
                    b=it.next;
                    if ~isempty(Gui.WindowAncestor(b)) 
                        if state ~= b.isSelected
                            b.setSelected(state);
                            try
                                priorActionCommand=b.getActionCommand;
                                b.setActionCommand(Plots.DOING_ALL);
                                ok=feval(callback, b, e, idx);
                                b.setActionCommand(priorActionCommand);
                                if ~ok
                                    break;
                                end
                            catch
                            end
                        end
                    end
                    idx=idx+1;
                end
                drawnow;
                Gui.HideBusy(sortGui.dlg)
            end
            
            
        end
        
        function [pnl, checkBoxes]=VerticalCheckBoxes(title, ...
                labels, dfltIdxs)
            jObj='javax.swing.JCheckBox';
            N=length(labels);
            pnl=Gui.SetTitledBorder(title);
            pnl.setLayout(javaObjectEDT('java.awt.GridLayout', N, 1));
            checkBoxes={};
            for i=1:N
                cbx=javaObjectEDT(jObj, labels{i});
                if ~isempty(find(dfltIdxs==i, 1))
                    cbx.setSelected(true);
                end
                checkBoxes{end+1}=cbx;
                pnl.add(cbx);
            end            
        end

        function [pnl, bg, radioButtons]=Vertical(title, ...
                labels, dfltIdx)
            jObj='javax.swing.JRadioButton';
            N=length(labels);
            pnl=Gui.SetTitledBorder(title);
            pnl.setLayout(javaObjectEDT('java.awt.GridLayout', N, 1));
            radioButtons={};
            bg=javax.swing.ButtonGroup;
            for i=1:N
                rb=javaObjectEDT(jObj, labels{i});
                if dfltIdx==i
                    rb.setSelected(true);
                end
                radioButtons{end+1}=rb;
                pnl.add(rb);
                bg.add(rb);
            end            
        end

        function [pnl, bg, radioButtons]=Horizontal(title, labels,...
                dfltIdx)
            jObj='javax.swing.JRadioButton';
            N=length(labels);
            if ~isempty(title)
                pnl=Gui.SetTitledBorder(title);
            else
                pnl=Gui.Panel;
            end
            
            bg=javax.swing.ButtonGroup;
            radioButtons={};
            for i=1:N
                rb=javaObjectEDT(jObj, labels{i});
                radioButtons{end+1}=rb;
                bg.add(rb);
                pnl.add(rb);
                if dfltIdx==i
                    rb.setSelected(true);
                end
            end
        end
        
        function Clear(bg)
            e=bg.getElements;
            while e.hasMoreElements
                rb=e.nextElement;
                if rb.isSelected
                    rb.setSelected(false);
                end
            end
        end
        
        function idx=Choice(bg)
            idx=1;
            e=bg.getElements;
            while e.hasMoreElements
                rb=e.nextElement;
                if rb.isSelected
                    return;
                end
                idx=idx+1;
            end
            idx=0;
        end
        
        
        function Disable(bg, idx)
            idx_=1;
            e=bg.getElements;
            while e.hasMoreElements
                rb=e.nextElement;
                if idx==idx_
                    rb.setEnabled(false);
                    return;
                end
                idx_=idx_+1;
            end
            idx_=0;
        end
        
        function [outerPnl, bg, checkBoxes, innerPnl]=Panel4(...
                javaArray, dfltIdxs, addScrollAfter, checkBox, gridBag)
            if ~checkBox
                checkBox=false;
                jObj='javax.swing.JRadioButton';
            else
                jObj='javax.swing.JCheckBox';
            end
            usingBag=false;
            N=javaArray.length;
            isJPanels=N>0 && isa(javaArray(1),'javax.swing.JPanel');
            innerPnl=javaObjectEDT('javax.swing.JPanel');
            innerPnl.setLayout(javax.swing.BoxLayout(innerPnl, 1));
            columns=1;
            if checkBox
                bg=javaObjectEDT('javax.swing.JPanel');
                if ~isJPanels
                    horizontal=true;
                    for i=1:N
                        if javaArray(i).length>2
                            horizontal=false;
                            break;
                        else
                        end
                    end
                else
                    columns=javaArray(1).getComponentCount;
                    horizontal=false;
                end
                if ~horizontal
                    if ~gridBag && ~isJPanels
                        innerPnl.setLayout(javaObjectEDT(...
                            'java.awt.GridLayout', N, 1));
                    else
                        innerPnl.setLayout(javaObjectEDT(...
                            'java.awt.GridBagLayout'));
                        gbc=javaObjectEDT('java.awt.GridBagConstraints');
                        %gbc.weighty=1;
                        gbc.anchor=gbc.WEST;
                        gbc.fill=0;
                        usingBag=true;
                    end
                else
                    innerPnl.setLayout(javaObjectEDT(...
                        'java.awt.GridLayout', 1, N));
                end
            else
                bg=javax.swing.ButtonGroup;
            end
            checkBoxes=java.util.ArrayList;
            idx=1;
            scrollToHeight=0;
            first=true;
            if ~checkBox || ~gridBag
                %javaArray is char or strings
                for i=1:N
                    if isJPanels
                        rb=javaObjectEDT(jObj, javaArray(i).getComponent(0).getText);
                    else
                        rb=javaObjectEDT(jObj, javaArray(i));
                    end
                    rObj=rb;
                    checkBoxes.add(rb);
                    bg.add(rb);
                    if ~usingBag
                        innerPnl.add(rObj);
                    else
                        gbc.gridy=i-1;
                        gbc.gridx=0;
                        innerPnl.add(rObj, gbc);
                    end
                    for j=2:columns
                        gbc.gridx=j-1;
                        gbc.weightx=.2;
                        jp=javaObjectEDT('javax.swing.JPanel');
                        jp.add(javaArray(i).getComponent(j-1));
                        innerPnl.add(jp, gbc);
                    end
                    if ~isempty(find(dfltIdxs==idx, 1))
                        rb.setSelected(true);
                        if first
                            d=rb.getPreferredSize;
                            scrollToHeight=d.height*(i-1);
                            first=false;
                        end
                    end
                    idx=idx+1;
                end                
            else
                try
                    for i=1:N
                        if isJPanels
                            rb=javaObjectEDT(jObj, javaArray(i).getComponent(0).getText);
                        else
                            rb=javaObjectEDT(jObj, javaArray(i));
                        end
                        checkBoxes.add(rb);
                        gbc.gridy=i-1;
                        gbc.gridx=0;
                        gbc.weightx=.8;
                        fprintf('y=%d x=%d wx=%d wy=%d fill=%d\n', ...
                            gbc.gridy, gbc.gridx, gbc.weighty, ...
                            gbc.weightx, gbc.fill);
                        innerPnl.add(rb, gbc);
                        for j=2:columns
                            gbc.gridx=j-1;
                            gbc.weightx=.2;
                            jp=javaObjectEDT('javax.swing.JPanel');
                            jp.add(javaArray(i).getComponent(j-1));
                            innerPnl.add(jp, gbc);
                        end
                        if ~isempty(find(dfltIdxs==idx, 1))
                            rb.setSelected(true);
                            if first
                                d=rb.getPreferredSize;
                                scrollToHeight=d.height*(i-1);
                                first=false;
                            end
                        end
                        idx=idx+1;
                    end
                catch ex
                    disp(['If javaArray contains JPanels then '...
                        'Radio.Panel2' newline 'expects component 0'...
                        ' to be JLabel and any 2nd component!' newline]);
                    ex.getReport
                end
            end
            outerPnl=javaObjectEDT('javax.swing.JPanel');
            outerPnl.setBorder(...
                javax.swing.BorderFactory.createEmptyBorder(1,8,1,4));
            if addScrollAfter<N
                d=innerPnl.getPreferredSize;
                d.width=floor(d.width*1.1);
                d.height=floor( (addScrollAfter/N)*d.height);
                innerPnl.setBorder([]);
                scroll=javaObjectEDT('javax.swing.JScrollPane', innerPnl);
                if scrollToHeight>0
                    scroll.getViewport.setViewPosition(...
                        java.awt.Point(0, scrollToHeight));
                end
                scroll.setPreferredSize(d);
                outerPnl=scroll;
            else                
                cb=javax.swing.BorderFactory.createCompoundBorder(...
                    javax.swing.BorderFactory.createBevelBorder(...
                    javax.swing.border.BevelBorder.RAISED), ...
                    javax.swing.BorderFactory.createEmptyBorder(8, 8, 8, 8));
                innerPnl.setBorder(cb);
                outerPnl.add(innerPnl);
            end
        end
 
    end
end