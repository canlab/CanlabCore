%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef ToolBarMethods
    methods(Static)
        
        function comboBox=addMnu(this, hObject, add2End, resetTo0, useJp)
            items=get(hObject, 'String');
            cb=@(c,e)ToolBarMethods.mnuCallback(hObject, c,e, resetTo0);
            comboBox=ToolBarMethods.addComboBox(this, items, cb, add2End, useJp);
        end

        function addComponent(this, comp)
            this.jToolbar.add(comp);
        end
        
        function pb=addButton(this, icon, tip, callback, text)
            if ~isempty(icon)
                p=fileparts(icon);
                if isempty(p)
                    icon=fullfile(BasicMap.Global.contentFolder, icon);
                end
                pb=javaObjectEDT('edu.stanford.facs.swing.ImageButton', ...
                    icon, tip);
            else
                pb=javaObjectEDT('javax.swing.JButton');
                pb.setMargin(java.awt.Insets(0,0,0,0));
                pb.setIconTextGap(0);
                pb.setBorderPainted(false);
                pb.setBorder([]);
                pb.setText([]);
                pb.setOpaque(false);
                %setBackground(UIManager.getColor("Panel.background"));
                pb.setFocusPainted(false);
                pb.setToolTipText(tip);
            end
            if nargin>4
                pb.setText(text);
                if ~isempty(icon)
                    pb.setIconTextGap(2);
                end
            end
            if String.EndsWith(class(this), 'MJToolBar')
                this.add(pb);
            else
                this.jToolbar.add(pb);
            end
            if ispc
                try
                    if String.EndsWith(class(this), 'MJToolBar')
                        this.add(javax.swing.JLabel(' '));
                    else
                        this.jToolbar.add(javax.swing.JLabel(' '));
                    end
                catch ex
                end
            end
            pb=handle(pb,'CallbackProperties');
            set(pb, 'ActionPerformedCallback',callback)
        end
        
        function addSeparator(this)
            this.jToolbar.addSeparator(java.awt.Dimension(5,5));
        end
        
        function addHelpButton(tb, helpID)
            box=javaObjectEDT('javax.swing.Box', 1);
            zz=box.createHorizontalGlue();
            ToolBarMethods.addComponent(tb,zz);
            ToolBarMethods.addButton(tb, ...
                fullfile(CytoGate.Get.contentFolder, 'help2.png'),...
                'Help me',...
                @(h,e)ToolBarMethods.showHelp(helpID));
        end     
        
        function showHelp(helpID)
             CytoGate.setHelp(helpID);
        end
        
        function comboBox=addComboBox(this, items, callback,...
                add2End, useJp, protoTypeItem, tip, startingIdx)
            if nargin<8
                startingIdx=0;
                if nargin<7
                    tip=[];
                    if nargin<6
                        protoTypeItem=0;
                    end
                end
            end
            drawnow;            
            comboBox=javaObjectEDT('javax.swing.JComboBox',items);
            width=comboBox.getPreferredSize().width;
            height=comboBox.getPreferredSize().height;
            comboBox.setMaximumSize(java.awt.Dimension(width,height));
            comboBox.setSelectedIndex(startingIdx);
            comboBox=handle(comboBox, 'CallbackProperties');
            set(comboBox, 'ActionPerformedCallback', callback);
            if nargin>4 && useJp
                jp=javaObjectEDT('javax.swing.JPanel');
                bl=java.awt.BorderLayout;
                jp.setLayout(bl);
                jp.add(comboBox, 'West');
                if ~isempty(this)
                if add2End
                    this.jToolbar(1).add(jp);
                else
                    this.jToolbar(1).add(jp,1);
                end
                end
            else
                if ~isempty(this)
                    if add2End
                        this.jToolbar(1).add(comboBox);
                    else
                        this.jToolbar(1).add(comboBox,1);
                    end
                end
            end
            %cb.setPrototypeDisplayValue(items(1));
            if ~isempty(this)
                this.firstItems{end+1}=items(1);
                if isempty(this.comboBoxes)
                    this.comboBoxes=comboBox;
                else
                    this.comboBoxes(end+1)=comboBox(1);
                end
            end
            if protoTypeItem>0
                comboBox.setPrototypeDisplayValue(items{protoTypeItem})
            end
            if ~isempty(tip)
                comboBox.setToolTipText(tip);
            end
        end

        function comboBox=addAutoComboBox(this, items, callback,...
                add2End, useJp, columns)
            drawnow;            
            comboBox=javaObjectEDT(...
                'edu.stanford.facs.swing.AutoComboBox', columns);
            try
                comboBox.setItems(items(2:end));
            catch ex
                items2={};
                N=length(items);
                for i=1:N
                    if ~isempty(items{i})
                        items2{end+1}=items{i};
                    end
                end
                comboBox.setItems(items2);
            end
            comboBox.setReadOnlyPrompt(items{1});
            comboBox=handle(comboBox, 'CallbackProperties');
            comboBox.setFindStartsWith(false);
            jb=javaObjectEDT('javax.swing.JButton');
            jbH = handle(jb, 'CallbackProperties');
            set(jbH, 'ActionPerformedCallback', callback);
            comboBox.setExitPressedButton(jb);
            if nargin>4 && useJp
                jp=javaObjectEDT('javax.swing.JPanel');
                bl=java.awt.BorderLayout;
                jp.setLayout(bl);
                jp.add(comboBox, 'West');
                if add2End
                    this.jToolbar(1).add(jp);
                else
                    this.jToolbar(1).add(jp,1);
                end
            else
                if add2End
                    this.jToolbar(1).add(comboBox);
                else
                    this.jToolbar(1).add(comboBox,1);
                end
            end
            this.firstItems{end+1}=items(1);
            if isempty(this.comboBoxes)
                this.comboBoxes=comboBox;
            else
                this.comboBoxes(end+1)=comboBox(1);
            end
        end


        % Drop-down (combo-box) callback function
        function mnuCallback(hObject, hCombo,~, resetTo0)
            itemIndex = get(hCombo,'SelectedIndex');  % 0=topmost item
            set(hObject, 'value', itemIndex+1);
            f=get(hObject,'Callback');
            feval(f, hObject, []);
            if resetTo0
                set(hCombo, 'SelectedIndex', 0);
            end
        end
        
        function resize(this)
            N=length(this.comboBoxes);
            for idx=1:N
                jCombo=this.comboBoxes(idx);
                jCombo.setPrototypeDisplayValue(this.firstItems(idx));
            end
            ToolBarMethods.refresh(this);
        end
        
        function refresh(this)
            this.jToolbar(1).repaint;
            this.jToolbar(1).revalidate;
        end
        
        function setSelectedIndex(cb, mnu, idx, resetTo0)
            if ~isempty(cb)
                set(cb, 'ActionPerformedCallback', []);
                set(cb, 'SelectedIndex', idx-1);
                set(cb, 'ActionPerformedCallback', ...
                    @(c,e)ToolBarMethods.mnuCallback(...
                    mnu, c,e, resetTo0));
            end
        end
        
        function setItemsFromMnu(cb, mnu)
            if ~isempty(cb)
                strs=get(mnu, 'String');
                ToolBarMethods.setItems(cb,strs);
            end
        end
        
        function setItems(comboBox,strs)
            apc=get(comboBox, 'ActionPerformedCallback');
            set(comboBox, 'ActionPerformedCallback', []);
            drawnow;
            %            javaMethodEDT('removeAllItems', comboBox);
            comboBox.removeAllItems();
            for i=1:length(strs)
                %javaMethodEDT('addItem', comboBox, strs{i});
                comboBox.addItem(strs{i});
            end
            drawnow;
            set(comboBox, 'ActionPerformedCallback', apc);
            comboBox.revalidate();
            comboBox.repaint();
        end
        
        function j=getJ(ht)
            j=get(get(ht,'JavaContainer'),'ComponentPeer');
        end
    end    
end