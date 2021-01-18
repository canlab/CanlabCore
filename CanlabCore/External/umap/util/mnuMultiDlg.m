%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function [idxs, cancelled, mainDlg, raIfRemembered, checkBoxes, allChb, sortGui]=...
    mnuMultiDlg(args, title, options, defaults, singleOnly, ...
    radioOrCheckBox, xtraCmp1, xtraWhere1, xtraCmp2, xtraWhere2, ...
    itemsPerScrollWindow,  where, subTitle)
cancelled=true;
raIfRemembered=[];
allChb=[];
if nargin<13
    subTitle=[];
    if nargin<12
        where='center';
        if nargin<11
            itemsPerScrollWindow=11;
            if nargin <10
                xtraWhere2='South';
                if nargin<9
                    xtraCmp2=[];
                    if nargin<8
                        xtraWhere1='South';
                        if nargin<7
                            xtraCmp1=[];
                            if nargin<6
                                radioOrCheckBox=false;
                                if nargin<5
                                    singleOnly=false;
                                    if nargin<4
                                        defaults=[];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
nOptions=length(options);
checkBoxes=[];
noCancel=isstruct(args) && isfield(args, 'noCancel') && args.noCancel;
disableChoices=isstruct(args) && isfield(args, 'disableChoices') ...
    && args.disableChoices;
if isstruct(args) && isfield(args, 'disableIdxs') 
    disableIdxs=args.disableIdxs;
else
    disableIdxs=[];
end
if isstruct(args) && isfield(args, 'allMsg')
    allMsg=args.allMsg;
elseif nOptions>1
    allMsg='All';
end
[guide, where, property, properties, defaults, icon, javaWindow, ...
    choiceTitle, checkFnc, modal, ~, ~, rememberId, ~, checkBoxFnc]=...
    decodeMsg(args, defaults, where);
if ~isempty(rememberId)
    ra=RememberedAnswers;
    idxs=ra.getIdxs(rememberId);
    if ra.isRemembered(rememberId)
        mainDlg=[];
        cancelled=false;
        raIfRemembered=ra;
        return;
    end
end
cancelled=true;
if ischar(options{1})
    jsa=javaArray('java.lang.String', nOptions);
    for i=1:nOptions
        jsa(i)=java.lang.String(options{i});
    end
else
    jsa=javaArray('java.lang.Object', nOptions);
    for i=1:nOptions
        jsa(i)=options{i};
    end
end
hasIcon=isempty(icon) || ~strcmp('none', icon);
if hasIcon
    ly=javaObjectEDT('java.awt.BorderLayout', 1, 1);
    pnl=javaObjectEDT('javax.swing.JPanel', ly);
    if radioOrCheckBox
        pnl.setBorder(javax.swing.BorderFactory.createEmptyBorder(4, 7, 4, 5));
    else
        pnl.setBorder(javax.swing.BorderFactory.createEmptyBorder(4, 7, 4, 5));
    end
    if isempty(icon)
        pnl.add(javaObjectEDT('javax.swing.JLabel', Gui.Icon('facs.gif')), 'West');
    else
        pnl.add(javaObjectEDT('javax.swing.JLabel', Gui.Icon(icon)), 'West');
    end
    pnl.add(javaObjectEDT('javax.swing.JLabel', Gui.Icon('blank.png')), 'East');
else
    ly=javaObjectEDT('java.awt.BorderLayout', 1, 1);
    pnl=javaObjectEDT('javax.swing.JPanel', ly);
    pnl.setBorder(javax.swing.BorderFactory.createEmptyBorder(2, 4, 2, 4));
end
if jsa.length<itemsPerScrollWindow
    itemsPerScrollWindow=jsa.length;
end
if radioOrCheckBox
    dfltIdxs=0;
    if ~isempty(defaults)
        if singleOnly
            dfltIdxs=defaults(1);
        else
            dfltIdxs=defaults+1;
        end
    end
    if isfield(args, 'sortProps') 
        [jsc,bg,checkBoxes, innerPnl]=Radio.Panel4(jsa, dfltIdxs, ...
            itemsPerScrollWindow, ~singleOnly, false);
    else
        [jsc,bg,checkBoxes, innerPnl]=Radio.Panel4(jsa, dfltIdxs, ...
            itemsPerScrollWindow, ~singleOnly, true);
    end
    if ~singleOnly && isa(jsa(1),'javax.swing.JPanel')
        priorCmps=Gui.GetJavaComponents(innerPnl);
    end
    if ~isempty(choiceTitle)
        Gui.SetTitledBorder(choiceTitle, jsc);
    end
    
    if ~isempty(disableIdxs)
        try
            nDisableIdxs=length(disableIdxs);
            for iii=1:nDisableIdxs
               di=disableIdxs(iii);
               checkBoxes.get(di-1).setEnabled(false);
            end
        catch ex
            ex.getReport
        end
    elseif disableChoices
        nCheckBoxes=checkBoxes.size;
        for iii=0:nCheckBoxes-1
            checkBoxes.get(iii).setEnabled(false);
        end
    end
else
    lst=javaObjectEDT('javax.swing.JList', jsa);
    if disableChoices
        lst.setEnabled(false);
    end
    lst.setBorder(javax.swing.BorderFactory.createEmptyBorder(4, 10, 4, 8));
    jListbox = handle(lst, 'CallbackProperties');    
    set(jListbox, 'MousePressedCallback',@myCallbackFcn);
    % Define the mouse-click callback function
    doubleClicked=false;
    scroll=javaObjectEDT('javax.swing.JScrollPane', lst);
    limit=itemsPerScrollWindow;
    lst.setVisibleRowCount(limit);
    d=lst.getPreferredSize;
    if d.width>650
        d.width=650;
        lst.setPreferredSize(d);
    end
    jsc=scroll;
end
if ~isempty(xtraCmp1) && ~strcmp(xtraWhere1,'south buttons')...
        && ~strcmp(xtraWhere1,'south west buttons')
    ly2=javaObjectEDT('java.awt.BorderLayout', 2, 8);
    pnl2=javaObjectEDT('javax.swing.JPanel', ly2);
    pnl2.add(xtraCmp1, xtraWhere1);
    if ~isempty(xtraCmp2)
        pnl2.add(xtraCmp2, xtraWhere2);
    end
    pnl.add(pnl2, 'Center');
    ly3=javaObjectEDT('java.awt.BorderLayout', 2, 10);
    pnl3=javaObjectEDT('javax.swing.JPanel', ly3);
    pnl2.add(pnl3, 'Center');
    pnl2=pnl3;
else
    pnl2=pnl;
end
allChbPnl=[];
if ~singleOnly && radioOrCheckBox
    ly3=javaObjectEDT('java.awt.BorderLayout');
    pnl3=javaObjectEDT('javax.swing.JPanel', ly3);
    noAll=isstruct(args) && isfield(args, 'noAll') && args.noAll;
    if ~noAll && nOptions>1
        allChb=javaObjectEDT('javax.swing.JCheckBox', allMsg);
        allChb.setMnemonic('a');
        if length(unique(defaults))==length(options)
            allChb.setSelected(true);
        end
        allH=handle(allChb,'CallbackProperties');
        allChbPnl=Gui.BorderPanel;
        allChbPnl.add(allChb, 'West');
        pnl3.add(allChbPnl, 'North');
        set(allH, 'ActionPerformedCallback', @(h,e)doAll(h,e));
        
    end
    pnl3.add(jsc, 'Center');
    checkBoxPnl=pnl3;
    pnl2.add(pnl3, 'Center');
else
    pnl2.add(jsc, 'Center');
end
pnl3=javaObjectEDT('javax.swing.JPanel');
if hasIcon
    pnl3.setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 28, 0, 0));
end
if ischar(guide)
    lbl=javaObjectEDT('javax.swing.JLabel', guide);
    lbl.setHorizontalAlignment(javax.swing.JLabel.CENTER);
    pnl3.add(lbl);
else
    pnl3.add(guide);
end
pnl2.add(pnl3, 'North');
if ~isempty(subTitle) && ischar(subTitle)
    lbl2=javaObjectEDT('javax.swing.JLabel', subTitle);
    lbl2.setHorizontalAlignment(javax.swing.JLabel.CENTER);
    pnl.add(lbl2, 'North');
end
if nargin<=3 || isempty(defaults)
    defaults=0;
end
if ~radioOrCheckBox
    if singleOnly
        lst.setSelectionMode(1);
    end
    if ~isempty(defaults)
        lst.setSelectedIndices(int32(defaults));
        lst.ensureIndexIsVisible(defaults(1));
    end
end

if isempty(javaWindow)
    jFrame = Gui.ParentFrame;
else
    jFrame = javaWindow;
end
mainDlg=javaObjectEDT('javax.swing.JDialog', jFrame);
dlg=handle(mainDlg, 'CallbackProperties');
set(dlg, 'WindowClosingCallback', @(h,e)windowClose());
if ~isempty(title)
    mainDlg.setTitle(title);
end
if ~noCancel
    done=javaObjectEDT('javax.swing.JButton', 'Ok');
else
    done=javaObjectEDT('javax.swing.JButton', 'Done');
end
doneH=handle(done,'CallbackProperties');
set(doneH, 'ActionPerformedCallback', @(h,e)close(true));
cancel=javaObjectEDT('javax.swing.JButton', 'Cancel');
cancel.setIcon(Gui.Icon('cancel.gif')); 
cancelH=handle(cancel,'CallbackProperties');
set(cancelH, 'ActionPerformedCallback', @(h,e)close(false));
edu.stanford.facs.swing.CpuInfo.registerEscape(mainDlg, cancel);
south=javaObjectEDT('javax.swing.JPanel', java.awt.BorderLayout);
southSouth=javaObjectEDT('javax.swing.JPanel', java.awt.BorderLayout);
southEast=javaObjectEDT('javax.swing.JPanel');
southSouth.add(southEast, 'South');
if ~isempty(rememberId)
    rememberCb=RememberedAnswers.GetCheckBox;
    southEast.add(rememberCb);
    southEast.add(Gui.Label(' '));
end
if ~noCancel
    southEast.add(cancel);
end
southEast.add(done);
south.add(southSouth, 'East');
if ~isempty(xtraCmp1) && strcmp(xtraWhere1,'south west buttons') 
    south.add(xtraCmp1, 'West');  
    if ~isempty(xtraCmp2)
        if strcmp('South', xtraWhere2)
            south.add(xtraCmp2, 'North');
        else
            pnl.add(xtraCmp2, xtraWhere2);
        end
    end
else
    if ~isempty(xtraCmp1) && strcmp(xtraWhere1,'south buttons')
        south.add(xtraCmp1, 'Center');
    end
end
mainDlg.getRootPane.setDefaultButton(done);
pnl.add(south, 'South');
pnl.setBorder(javax.swing.BorderFactory.createEmptyBorder(2, 5, 10, 5));
mainDlg.add(pnl);
sortGui=[];
if radioOrCheckBox 
    if isempty(allChbPnl)
        allChbPnl=Gui.BorderPanel;
        pnl3.add(allChbPnl, 'North');
        allMsg='All';
        allChb=[];
    end
    if (singleOnly && nOptions>7) || (~singleOnly && nOptions>2)
        sortGui=SortGui(mainDlg, allChb, allMsg, allChbPnl, options, ...
            checkBoxes, innerPnl);
        if isa(jsa(1),'javax.swing.JPanel')
            sortGui.fncRefresh=@refreshPanelOrder;
        end
        if isfield(args, 'sortProps') 
            searchDflt=isfield(args, 'sortSearch') && args.sortSearch;
            sortGui.setProperties(args.sortProps, args.sortProp, searchDflt);
        end
    end
else
    sortGui=[];
end
mainDlg.pack;
if ~isempty(jFrame)
    mainDlg.setLocationRelativeTo(jFrame);
end
Gui.LocateJava(mainDlg, javaWindow, where);
if radioOrCheckBox
    nCh=checkBoxes.size;
    for iCh=1:nCh
        cb1=checkBoxes.get(iCh-1);
        cb2= handle(cb1, 'CallbackProperties');
        set(cb2, 'ActionPerformedCallback', @(h,e)innerChbCb(h,e));
    end
end
cancelled=true;
if ~ispc
    setAlwaysOnTopTimer(mainDlg)
end
MatBasics.RunLater(@(h,e)noTip(), .25)
if ~isempty(sortGui)
    sortGui.setAllChbText;
end
mainDlg.setModal(modal);
Gui.SetJavaVisible(mainDlg);
drawnow;
conclude;

    function noTip
        BasicMap.Global.closeToolTip;
    end

    function conclude
        if ~radioOrCheckBox
            javaIdxs=lst.getSelectedIndices;
            N=length(javaIdxs);
        else
            doubleClicked=false;
            N=0;
        end
        idxs=[];
        if doubleClicked || ~cancelled
            if N==0
                if radioOrCheckBox
                    if singleOnly
                        idxs=Radio.Choice(bg);
                    else
                        idxs=getSelectedIdxs;
                    end
                else
                    %idxs=1;
                end
            else
                for i=1:N
                    idx=javaIdxs(i)+1;
                    idxs(end+1)=idx;
                    answer=StringArray.IndexOf(options, idx);
                    disp(answer);
                end
            end
        end
        if ~cancelled
            if ~isempty(property) && ~isempty(properties)
                if singleOnly
                    properties.set(property, num2str(idxs));
                else
                    properties.set(property, num2str(idxs-1));
                end
            end
            if ~isempty(rememberId)
                if rememberCb.isSelected
                    rememberAnswer(idxs);
                end
            end
        end
    end

    function windowClose
        close(false);
    end

    function doAll(h,e)
        isSelected=h.isSelected;
        N2=checkBoxes.size;
        for ii=1:N2
            cb=checkBoxes.get(ii-1);
            if ~isempty(Gui.WindowAncestor(cb))
                cb.setSelected(isSelected);
            end
        end
        innerChbCb(h,e);
    end

    function close(saved)
        if isstruct(args) && isfield(args, 'closeFnc')
            feval(args.closeFnc, saved, idxs, checkBoxes);
            return;
        end
        cancelled=~saved;
        if ~isempty(checkFnc)
            conclude;
            ok=feval(checkFnc, idxs, cancelled, mainDlg);
            if ~ok
                return;
            end
        end
        mainDlg.dispose;
    end

    function myCallbackFcn(jListbox,jEventData)
        % Determine the click type
        % (can similarly test for CTRL/ALT/SHIFT-click)
        if jEventData.getClickCount==2
            w=Gui.Wnd(lst);
            if ~isempty(w)
                doubleClicked=true;
                if modal
                    w.dispose;
                else
                    close(true);
                end
            end
        end
    end

    function rememberAnswer(idxs)
        curFig=get(0, 'currentFigure');
        [~, ~, ~, ~, quadrant]=Gui.FindScreen(curFig);
        ra=RememberedAnswers;
        if strcmpi('west', quadrant{2})
            where='north east++';
        else
            where='north west++';
        end
        Gui.LocateJava(mainDlg, javaWindow, where);
        mainDlg.setModal(false);
        mainDlg.setVisible(true);
        
        ra.remember(rememberId, rememberCb, idxs);
        raIfRemembered=ra;
    end

    function [idxs_, N]=getSelectedIdxs
        [idxs_, N]=Gui.GetSelectedChbIdxs(checkBoxes);
    end

    function innerChbCb(h, e)
        if ~isempty(sortGui)
            idxs_=sortGui.setAllChbText;
        else
            idxs_=Gui.GetSelectedChbIdxs(checkBoxes);
        end
        if ~isempty(checkBoxFnc)
            feval(checkBoxFnc, h, e, idxs_, checkBoxes);
        end
    end

    function innerPnl=refreshPanelOrder
        checkBoxPnl.remove(jsc);
        [jsc,~,~, innerPnl]=Radio.Panel2(jsa, [], ...
            itemsPerScrollWindow, true, [], priorCmps, ...
            sortGui.sortIdxs, sortGui.visibleIdxs);
        checkBoxPnl.add(jsc, 'Center');
        jsc.repaint;
    end

end
