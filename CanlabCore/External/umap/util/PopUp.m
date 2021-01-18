%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

classdef PopUp < handle
    methods(Static)
        function l=CLOSE_LABEL
            l=['<html>' BasicMap.Global.smallStart ...
                'Close' BasicMap.Global.smallEnd '</html>'];
        end
    end
    properties(Constant)
        PANE_TIMEOUT=1231;
        PANE_LEAVE_BTN=java.awt.Color(.88, .88, .95);
    end
    properties(GetAccess=private)
        busy=[];
        cancelFunction=[];
        msg1=[];
        msg2=[];
    end
    
    properties    
        cancelled=false;
        showTimeRemaining=true;
    end
    
    properties(SetAccess=private)
        priorPane;
        main=[];
        taskName;
        titledBorder;
        progressRate;
        answer;
        tock;
        overallTock;
        progressState=0;
        dlg=[];
        biggestD;
        cancelBtn;
        pb;
        label=[];
        priorFig;
        lastSpentSecs;
        lastReportedRemainingSecs=0;
        lastReportedSpentSecs=0;
    end
    methods
        function incrementProgress(this, by, progressDescription)
            if ~isempty(this.pb)
                if nargin<2
                    by=1;
                end
                if isempty(this.progressRate)
                    if this.pb.getMaximum>1000
                        this.progressRate=floor(this.pb.getMaximum/200);
                    else
                        this.progressRate=1;
                    end
                end
                if this.progressRate>2
                    this.progressState=this.progressState+by;
                    if this.progressState<this.progressRate
                        return;
                    end
                    by=this.progressState;
                    this.progressState=0;
                end
                this.pb.setValue(by+this.pb.getValue);
                if this.showTimeRemaining
                    secs=toc(this.tock);
                    [remaining, strRemaining]=PopUp.Remaining(secs, ...
                        this.pb.getValue/this.pb.getMaximum, ...
                        this.lastReportedRemainingSecs);
                    if ~isempty(strRemaining)
                        if remaining>0
                            if isempty(this.taskName)
                                this.pb.setString([ strRemaining ' remaining']);
                            else
                                this.pb.setString([ strRemaining ' left in '...
                                    this.taskName]);
                            end
                        end
                        this.lastReportedRemainingSecs=remaining;
                    end
                end
                if ~isempty(this.overallTock)
                    this.showTimeSpent;
                end
            end
        end
        
        function setTimeSpentTic(this, tock)
            if nargin<2
                tock=tic;
            end
            this.overallTock=tock;
            fs=11;
            try
                app=BasicMap.Global;
                if app.highDef
                    fs=18;
                end
            catch ex
            end
            
            [~, tp]=Gui.SetTitledBorder('Time spent', this.main,...
                java.awt.Font('Courier', java.awt.Font.PLAIN, fs),...
                java.awt.Color.BLACK);
            tp.setTitlePosition(tp.BOTTOM);
            tp.setTitleJustification(tp.RIGHT);
            this.titledBorder=tp;
            this.dlg.pack;
        end
        
        function showTimeSpent(this)
            secs=toc(this.overallTock);
            if this.pb.getValue/this.pb.getMaximum>.97
                strSpent=String.TimeReport(secs, 0);
            else
                strSpent=String.TimeReport(secs, this.lastReportedSpentSecs);
                if isempty(strSpent)
                    return;
                end
            end
            this.titledBorder.setTitle([strSpent ' spent']);
            this.dlg.getComponent(0).repaint;
            this.lastReportedSpentSecs=secs;
        end
        
        function initProgress(this, cnt, taskName)
            if nargin>2
                idx=String.IndexOf(taskName, '<hr>');
                if idx>0
                    taskName=taskName(1:idx-1);
                end
                taskName=char(...
                    edu.stanford.facs.swing.Basics.RemoveXml(taskName));
                if ~isequal(taskName, this.taskName)
                    this.taskName=taskName;
                else
                    this.taskName=[];
                end
            else
                this.taskName=[];
            end
            if cnt<=0
                cnt=1;
            end
            this.progressRate=[];
            if ~isempty(this.pb)
                this.pb.setMaximum(cnt);
                this.pb.setValue(0)
            else
                %color of string painted inside bar
                javax.swing.UIManager.put(...
                    'ProgressBar.selectionForeground',...
                    java.awt.Color(.04, .11, .75));
                %color of BAR
                javax.swing.UIManager.put(...
                    'ProgressBar.foreground',...
                    java.awt.Color(.72, 1, .72));
                %color of string painted OUTSIDE bar
                javax.swing.UIManager.put(...
                    'ProgressBar.selectionBackground',...
                    java.awt.Color(.41, .41, .91));
                %color of ?????
                javax.swing.UIManager.put(...
                    'ProgressBar.background',...
                    java.awt.Color(.8, .1, .2));
                
                this.pb=javaObjectEDT('javax.swing.JProgressBar', 0, cnt);
                if ismac
                    this.pb.setUI( javax.swing.plaf.basic.BasicProgressBarUI );
                end
                D=this.pb.getPreferredSize;
                if ispc
                    if BasicMap.Global.highDef
                    else
                        rightHeight=20;
                        if D.height<rightHeight
                            D.height=rightHeight;
                            this.pb.setPreferredSize(D);
                        end
                    end
                end
                jp=javaObjectEDT('javax.swing.JPanel', ...
                    javaObjectEDT('java.awt.BorderLayout', 2, 8));
                jp.add(this.label, 'Center');
                jp.add(this.pb, 'South');
                this.main.add(jp, 'Center');
                
                fs=12;
                try
                    app=BasicMap.Global;
                    if app.highDef
                        fs=18;
                    end
                catch ex
                end
                this.pb.setFont(java.awt.Font('Courier', 1, fs));
                
                this.stop;
                this.dlg.pack;
            end
            this.tock=tic;
            this.lastReportedRemainingSecs=0;
            this.pb.setStringPainted(true);
        end
        
        function cancel(this)
            this.cancelled=true;
            this.setText('Cancelling...');
            if ~isempty(this.pb)
                mx=this.pb.getMaximum;
                this.pb.setMaximum(mx-2);
            end
        end

        function this=PopUp(msg, where, title, showBusy, cancelBoolOrFnc, ...
                icon, modal, priorPu)
            if nargin<8
                priorPu=[];
                if nargin<7
                    modal=false;
                    if nargin<6
                        icon=[];
                        if nargin<5
                            cancelBoolOrFnc=[];
                            if nargin<4
                                showBusy=true;
                                if nargin<3
                                    title='Note ....';
                                    if nargin<2
                                        where='center';
                                        if nargin<1
                                            msg='One moment please ...';
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            this.priorFig=get(0, 'currentFigure');
            this.main=Gui.BorderPanel(5, 5);
            if isempty(priorPu)
                this.main.setBorder(javax.swing.BorderFactory.createEmptyBorder (12,12,12,12));
                jd=javaObjectEDT('javax.swing.JDialog', Gui.ParentFrame);
                if ~isempty(title)
                    jd.setTitle(title);
                end
                jd.setContentPane(this.main);
                this.dlg=jd;
            else
                jd=priorPu.dlg;
                this.dlg=jd;
                this.pushPopUp(title);
            end
            if ischar(msg)
                jl=javaObjectEDT('javax.swing.JLabel');
                jl.setText(msg);
                this.msg1=String.RemoveXml(msg);
                this.msg2=msg;
                this.main.add(jl, 'Center');
                this.label=jl;
            else
                this.main.add(msg, 'Center');
            end
            if showBusy
                busy=Gui.JBusy('..');
                this.main.add(busy.getComponent, 'West');
            else
                busy=[];
            end
            this.setCancel(cancelBoolOrFnc);
            if ~showBusy && isempty(icon)
                try
                    icon=Gui.Icon('smallGenie.png');
                catch ex
                    icon=[];
                end
            end
            if ~isempty(icon)
                this.setIcon(icon);
            else
                jd.pack;
            end
            try
                setAlwaysOnTopTimer(jd);
            catch
            end
            if isempty(priorPu)
                try
                    if ~isempty(this.priorFig) && Gui.IsVisible(this.priorFig)
                        Gui.LocateJava(jd, Gui.JFrame(this.priorFig), where);
                    else
                        Gui.LocateJava(jd, [], where);
                    end
                catch
                end
            end
            if ~modal
                Gui.SetJavaVisible(jd);
            else
                jd.setModal(true);
            end
            if ~isempty(busy)
                busy.start;
            end
            this.busy=busy;
            drawnow;
        end
        
        function pushPopUp(this, title)
            this.priorPane=this.dlg.getContentPane;
            whole=Gui.BorderPanel;
            whole.add(this.priorPane, 'North');
            center=Gui.Panel;
            if ~isempty(title)
                Gui.SetTitledBorder(title, center);
            else
                Gui.SetTitledBorder('Sub tasks...', center);
            end
            center.add(this.main);
            whole.add(center, 'Center');
            this.dlg.setContentPane(whole);
        end
        
        function popPopUp(this)
            this.dlg.setContentPane(this.priorPane);
        end
        
        function stop(this)
            if ~isempty(this.busy)
                this.busy.stop;
                this.main.remove(this.busy.getComponent);
            end
        end
        
        function close(this, force, refreshOld)
            if nargin<2
                force=false;
            end
            try
                globalPu=BasicMap.Global.pu;
            catch ex
                globalPu=[];
            end
            if force || isempty(globalPu) || this.dlg ~= globalPu.dlg
                this.stop;
                if isempty(this.priorPane)
                    this.dlg.dispose;
                else
                    this.popPopUp;
                    this.dlg.pack;
                end
                if (nargin<3 || refreshOld) && ~isempty(this.priorFig)
                    if ~ishandle(this.priorFig)
                        this.priorFig=get(0, 'currentFigure');
                    end
                    if ~isempty(this.priorFig)
                        if strcmpi('off', get(this.priorFig, 'visible'))
                            %disp('prior figure is invibisble');
                        else
                            figure(this.priorFig);
                        end
                    end
                end
            end
        end
        
        function packIfNeedBe(this)
            this.label.setPreferredSize([])
            drawnow;
            d=this.label.getPreferredSize;
            if isempty(this.biggestD)
                this.biggestD=d;
                this.dlg.pack;
            elseif d.width>this.biggestD.width || ...
                    d.height>this.biggestD.height
                if d.width>this.biggestD.width 
                    this.biggestD.width=d.width;
                end
                if d.height>this.biggestD.height
                    this.biggestD.height=d.height;
                end
                this.label.setPreferredSize(this.biggestD);
                this.dlg.pack;
            else
                difs=[abs(d.width-this.biggestD.width) ...
                    abs(d.height-this.biggestD.height)];
                mx=max(difs);
                if mx<10
                    this.dlg.pack;
                else
                    this.main.invalidate;
                end
            end
        end
        
        function setText(this, msg)
            this.msg1=String.RemoveXml(msg);
            this.msg2=msg;
            this.label.setText(msg);
            this.packIfNeedBe;
        end

        function setText2(this, msg2)
            this.label.setText(['<html>' this.msg1 '<hr><br>' msg2 '</html>']);
            this.packIfNeedBe;
            drawnow;
        end
        
        function setText3(this, msg3)
            N=length(this.msg2);
            if N>15 && strcmpi('</center></html>', this.msg2(end-15:end))
                this.label.setText(['<html>' this.msg2(7:end-16) ...
                    ' ' msg3 '</html>']);
            elseif N>6 && strcmpi('</html>', this.msg2(end-6:end))
                this.label.setText(['<html>' this.msg2(7:end-7)  ...
                    ' ' msg3 '</html>']);
            else
                this.label.setText(['<html>' this.msg2 '<hr>' ...
                    msg3 '</html>']);
            end
            this.packIfNeedBe;
            drawnow;
        end

        function setIcon(this, icon)
            if ischar(icon)
                icon=Gui.Icon(icon);
            end
            if ~isempty(this.label)
                this.label.setIcon(icon);
                this.label.setIconTextGap(9);
            end
            this.dlg.pack;
        end

        function setVisible(this, on)
            this.dlg.setVisible(on);
        end
    
        function addCloseBtn(this, txt)
            if nargin<2
                txt='Ok';
            end
            jp=Gui.Panel;
            btn=Gui.NewBtn(txt, @(h,e)shut);
            jp.add(btn);
            this.main.add(jp, 'South');
            this.dlg.pack;
            function shut
                this.dlg.dispose;
            end
        end
        
        function addYesNo(this, cancelToo, defaultAnswer)
            if nargin<2
                cancelToo=true;
            end
            jp=Gui.Panel;
            btnYes=Gui.NewBtn('Yes', @(h,e)close(1));
            jp.add(btnYes);
            btnNo=Gui.NewBtn('No', @(h,e)close(0));
            jp.add(btnNo);
            btnCancel=Gui.NewBtn('Cancel', @(h,e)close(-1));
            root=this.dlg.getRootPane;
            if cancelToo
                jp.add(btnCancel);
                Gui.RegisterEscape(root, btnCancel);
            end
            this.main.add(jp, 'South');
            if nargin>2
                if defaultAnswer==1
                   root.setDefaultButton(btnYes);
                elseif defaultAnswer==-1
                    root.setDefaultButton(btnCancel);
                else
                    root.setDefaultButton(btnNo);
                end
            end
            this.dlg.pack;
            function close(answ)
                this.answer=answ;
                this.dlg.dispose;                
            end
        end
        
        function setCancel(this, cancelBoolOrFnc)
            if islogical(cancelBoolOrFnc) 
                if cancelBoolOrFnc
                    cancelBoolOrFnc=@(h,e)cancel(this);
                else
                    return;
                end
            end
            this.cancelFunction=cancelBoolOrFnc;
            if ~isempty(this.cancelFunction)
                c=javaObjectEDT('javax.swing.JButton', 'Cancel');
                ch=handle(c,'CallbackProperties');
                set(ch, 'ActionPerformedCallback', ...
                    this.cancelFunction);
                jp=javaObjectEDT('javax.swing.JPanel');
                jp.add(c);
                this.main.add(jp, 'South');
                this.cancelBtn=c;
            end
        end
        
        function pack(this)
            this.dlg.pack;
        end
        
        function setAlwaysOnTop(this, ok)
            javaMethodEDT( 'setAlwaysOnTop', this.dlg, ok);
            tmr=timer;
            tmr.StartDelay=1.5;
            tmr.TimerFcn=@(h,e)act;
            start(tmr);
            
            function act
                javaMethodEDT( 'setAlwaysOnTop', this.dlg, ok);
            end
        end
    end
    
    methods(Static)
        function [remaining, strRemaining]=Remaining(secs, percentDone,...
                lastReportedRemainingSecs)
            totalSecs=secs/percentDone;
            remaining=totalSecs-secs;
            guessFactor=abs(percentDone-1)/2;
            remaining=remaining*(1+guessFactor);
            strRemaining=String.TimeReport(remaining, ...
                lastReportedRemainingSecs, true);
        end
        
        function KeepPaneBtnText(pane)
            pane.setBackground(PopUp.PANE_LEAVE_BTN)
        end
        
        function TimedClose(jd, pauseSecs, pane, stripActions)
            if nargin<4
                stripActions=true;
                if nargin<3
                    pane=[];
                end
            end
            closeNow=false;            
            title='';
            firstBtnAl=[];
            clickCount=0;
            
            btn=javaObjectEDT('javax.swing.JButton');
            btnClass=btn.getClass;
            if pauseSecs>2
                tmr=timer;
                tmr.StartDelay=pauseSecs;
                tmr.TimerFcn=@(h,e)closeit;
                start(tmr);
                
                tmr=timer;
                if isdeployed || ~ismac
                    tmr.StartDelay=.5;
                else
                    tmr.StartDelay=.8;
                end
                pauseSecs=pauseSecs-2;
                tmr.Period=1;
                tmr.TasksToExecute=pauseSecs;
                tmr.ExecutionMode='fixedRate';
                tmr.TimerFcn=@(h,e)countDown;
                title=char(jd.getTitle);
                start(tmr);
            end
            
            firstBtn=[];
            txt1='';
            txt2='';
            try
                if BasicMap.Global.is('showPopUp', true)
                    Gui.SetJavaVisible(jd);
                end
            catch ex
            end
            drawnow;
            function closeit
                if ~closeNow
                    if ~isempty(firstBtnAl)
                        try
                            feval(firstBtnAl);
                        catch ex
                            feval(firstBtnAl, firstBtn, []);                        
                        end
                    end
                    if ~isempty(pane)
                        pane.setValue(PopUp.PANE_TIMEOUT);
                    end
                     jd.setVisible(false);
                end
            end
            function click
                clickCount=clickCount+1;
                if clickCount==1
                    stop(tmr);
                    delete(tmr);
                    if ismac
                        firstBtn.setText('<html><font color="#ADFF2F"><b>Close</b></font></html>');
                    else
                        firstBtn.setText('Close');
                    end
                    closeNow=true;
                    jd.setTitle(title);
                else
                    if ~isempty(firstBtnAl)
                        feval(firstBtnAl);
                    end
                    jd.setVisible(false);
                end
            end
            function countDown
                pauseSecs=pauseSecs-1;
                secs=[num2str(pauseSecs) ' secs'];
                jd.setTitle([ title ' (closes in ' secs ')']);
                if ~isempty(pane)
                    if isequal(PopUp.PANE_LEAVE_BTN, pane.getBackground)
                        return;
                    end
                end
                firstCount=false;
                if isempty(firstBtn)
                    firstCount=true;
                    firstBtn=Gui.FindFirst(jd, btnClass, 'Ok');
                    if ~isempty(firstBtn)
                        try
                            app=BasicMap.Global;
                        catch ex
                            app.smallStart='<small>';
                            app.smallEnd='</small>';
                        end
                        txt1=firstBtn.getText;
                        try
                            fbc=handle(firstBtn, 'CallbackProperties');
                            set(fbc, 'MouseEnteredCallback', ...
                                @(h,e)stopCountingDown());
                            jcmp=jd.getComponent(0).getComponent(1);
                            fbc=handle(jcmp, ...
                                'CallbackProperties');
                            set(fbc, 'MouseEnteredCallback', ...
                                @(h,e)stopCountingDown());
                        catch ex
                        end
                        if strcmpi('ok', char(txt1))
                            txt1='Closing in';
                        end
                        if isempty(txt1)
                            txt1='';
                        else
                            Html.remove(txt1);
                        end
                        if ismac
                            txt1=['<html><font color="#8DDF2F">' ...
                                txt1 ' <b>' app.smallStart];
                            txt2=[app.smallEnd '</font></b></html>'];
                        else
                            txt1=['<html>' txt1 ' <b>' app.smallStart];
                            txt2=[app.smallEnd '</b></html>'];
                            
                        end
                    end
                    firstBtnAl=get(handle(firstBtn, 'CallbackProperties'), ...
                        'ActionPerformedCallback');
                    if stripActions
                        firstBtnAls=firstBtn.getActionListeners;
                        nFirstBtnAls=length(firstBtnAls);
                        for i=1:nFirstBtnAls
                            firstBtn.removeActionListener(firstBtnAls(i));
                        end
                        drawnow;
                        set(handle(firstBtn, 'CallbackProperties'), ...
                            'ActionPerformedCallback', @(h,e)click);
                    end
                end
                if ~isempty(firstBtn)
                    firstBtn.setText([txt1 secs txt2]);
                    if isdeployed || ~firstCount || ~ismac || ...
                            1==size(get(0, 'MonitorPositions'), 1)
                        jd.pack;
                    else
                        old=jd.getLocation;
                        MatBasics.RunLater(@(h,e)handleMac(old), .2);
                        jd.pack;
                        drawnow;
                    end
                end
            end
            
            function handleMac(old)
                if ~old.equals(jd.getLocation) %wierd issue with java 7 and MAC 2 phycical screens
                    jd.setLocation(old);
                end
            end
            function stopCountingDown()
                if clickCount==0
                    firstBtn.setToolTipText(Html.WrapHr([...
                        'Automatic close halted '...
                        '<br>... click "Close" to close']));
                    try
                        CytoGate.Get.showToolTip(firstBtn);
                    catch ex
                    end
                    firstBtn.requestFocus;
                    click;
                end
            end
        end
        
        function this=New(msg, where, title, showBusy, cancelFnc, icon)
            this=BasicMap.Global.pu;
            if isempty(this) || ~this.dlg.isVisible
                if nargin<6
                    icon=[];
                    if nargin<5
                        cancelFnc=[];
                        if nargin<4
                            showBusy=true;
                            if nargin<3
                                title='Note ....';
                                if nargin<2
                                    where='center';
                                    if nargin<1
                                        msg='One moment please ...';
                                    end
                                end
                            end
                        end
                    end
                end
                if ~showBusy && isempty(icon)
                    icon=Gui.Icon('smallGenie.png');
                end
                this=PopUp(msg,where,title,showBusy,cancelFnc,icon);
            else
                if nargin>=5
                    this.setCancel(cancelFnc);
                end
                if nargin>=3
                    this.dlg.setTitle(title);
                end
                if nargin>=6
                    this.label.setIcon(icon);
                end
                if nargin>=1
                    this.label.setText(msg);
                    this.dlg.pack;
                end
            end
        end
        
        function jMenu=Menu
            jMenu = javaObjectEDT(javax.swing.JPopupMenu);
            ff=jMenu.getFont;
            if ispc
                f2=java.awt.Font('Arial', ff.getStyle, ff.getSize+1);
            else
                f2=java.awt.Font('Arial', ff.getStyle, ff.getSize);
            end
            jMenu.setFont(f2);
        end
        
        function jd=Pane(pane, title, where, javaWin, modal, pauseSecs, ...
                suppressParent, stripActions, widthLimit)
            if nargin<8
                stripActions=true;
            end
            isFreeFloating=false;
            if nargin<4 || isempty(javaWin)
                if nargin<7 || ~suppressParent
                    javaWin=Gui.ParentFrame;
                else
                    isFreeFloating=true;
                end
            end
            jd=javaMethodEDT('createDialog', pane, javaWin, title);
            if nargin>8
                sz=jd.getSize;
                if sz.width>widthLimit
                    jd.setSize(java.awt.Dimension(widthLimit, sz.height))
                end
            end
            jd.setResizable(true);
            app=BasicMap.Global();
            if nargin>4
                jd.setModal(modal);
            end
            %disp(['app.parentCmpForPopup=' app.parentCmpForPopup]);
            if ~isempty(app.parentCmpForPopup)
                javaMethodEDT( 'setLocationRelativeTo', jd, ...
                    app.parentCmpForPopup);
                if nargin>2
                    disp('YES');
                    Gui.LocateJava(jd, app.parentCmpForPopup, where);
                end
            elseif ~isempty(javaWin)
                javaMethodEDT( 'setLocationRelativeTo', jd, javaWin);
                if nargin>2
                    Gui.LocateJava(jd, javaWin, where);
                end
            elseif isFreeFloating
                if nargin>2
                    Gui.LocateJava(jd, Gui.ParentFrame, where);
                end
            end
            if ~ispc %Florian and others note this is a problem 
                setAlwaysOnTopTimer(jd);
            end
            if nargin>=6 && ~isempty(pauseSecs)
                PopUp.TimedClose(jd, pauseSecs, pane, stripActions);
            else
                try
                    if BasicMap.Global.is('showPopUp', true)
                        Gui.SetJavaVisible(jd);
                    end
                catch ex
                end
                
            end
            BasicMap.Global.closeToolTip;
        end
        
        
        function SetText(pu, txt)
            if ~isempty(pu)
                pu.setText(txt);
            end
        end
        
        function InitProgress(pu, N)
            if ~isempty(pu)
                pu.initProgress(N);
            end
            
        end
        
        function Increment(pu)
            if ~isempty(pu)
                pu.incrementProgress;
            end
        end
    end
end
