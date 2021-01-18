%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function [answer, yes, cancelled, raIfRemembered, jd]...
    =questDlg(theMsg, Title, varargin)
if nargin<2
    Title='Confirm...';
end
rememberCb=[];
raIfRemembered=[];
rememberId=[];
msgType=javax.swing.JOptionPane.QUESTION_MESSAGE;
if nargin<3
    if nargout>2
        varargin={'Yes', 'No', 'Cancel', 'Yes'};
    else
        varargin={'Yes', 'No', 'Yes'};
    end
end 
[msgType, jsa, default]=getMsgTypeAndOptions(msgType, varargin);
[theMsg, where, property, properties, default, myIcon, javaWin,~, checkFnc,...
     modal, pauseSecs, southWestComponent, rememberId, rememberOnly]=...
     decodeMsg(theMsg, default);
qu=javax.swing.JOptionPane.YES_NO_OPTION;
pane=javaObjectEDT('javax.swing.JOptionPane', theMsg, msgType, qu);
adds=java.util.ArrayList;
if ~isempty(rememberId)
    rememberCb=RememberedAnswers.GetCheckBox;
    adds.add(rememberCb);         
end
if ~isempty(southWestComponent)
    adds.add(southWestComponent);
end
if adds.size>0
    N=jsa.length;
    if ispc
        objs=javaArray('java.lang.Object', N+1+adds.size);
        for i=1:N
            objs(i)=jsa(i);
        end
        objs(N+1)=javax.swing.JLabel('<html>&nbsp;&nbsp;&nbsp;</html>');
        for i=0:adds.size-1
            objs(end-i)=adds.get(i);
        end
    else
        objs=javaArray('java.lang.Object', N+adds.size);
        for i=0:adds.size-1
            objs(end-i)=adds.get(i);
        end
        for i=1:N
            objs(i)=jsa(i);
        end
    end
    pane.setOptions(objs);
else
    pane.setOptions(jsa);
end
if ~isempty(rememberId)
    ra=RememberedAnswers;
    idx=ra.get(rememberId);
    if idx>0
        answer=jsa(idx);
        yes=strcmpi('Yes',answer) || strcmpi('Ok', answer);
        cancelled=false;
        disp(['Remembered answer used "' char(answer) '"']);
        raIfRemembered=ra;
        return;
    end
    idx=BasicMap.Global.getNumeric(rememberId, 0);
    if idx>0 && idx<=length(jsa)
        default=jsa(idx);
    end
end
pane.setInitialValue(default);
if ~strcmp('none', myIcon)
    if msgType==0
        myIcon='error.png';
    elseif msgType==1
        myIcon='facs.gif';
    elseif msgType==2
        myIcon='warning.png';
    elseif isempty(myIcon)
        myIcon='facs.gif';
    end
else
    myIcon='blank.png';
end
if ~isempty(myIcon)
    pane.setIcon(Gui.Icon(myIcon));
end
MatBasics.RunLater(@(h,e)noTip(), .35)
jd=PopUp.Pane(pane, Title, where, javaWin, modal, pauseSecs, ...
    false, isempty(checkFnc));
answer=pane.getValue;
if isnumeric(answer)
    answer='';
end
if strcmp(answer, 'uninitializedValue')
    answer=char(default);
end
yes=strcmpi('Yes',answer) || strcmpi('Ok', answer);
cancelled=isempty(answer) || strcmpi('Cancel', answer);
if ~isempty(property)
    if ~cancelled
        properties.set(property, answer);
    end
end
if ~cancelled && ~isempty(rememberCb)
    idx=indexOf(answer);
    BasicMap.Global.set(rememberId, num2str(idx));
    if rememberCb.isSelected
        rememberAnswer;
    end
end
if ~modal && ~isempty(checkFnc)
    jd.setResizable(true);
    dlg=handle(jd, 'CallbackProperties');
    set(dlg, 'WindowClosingCallback', @(h,e)close([]));
    
    btn=javaObjectEDT('javax.swing.JButton');
    btnClass=btn.getClass;
    nChoices=length(jsa);
    for ii=1:nChoices
        answ=char(jsa(ii));
        btn2=Gui.FindFirst(jd, btnClass, answ);
        btnAls=btn2.getActionListeners;
        nFirstBtnAls=length(btnAls);
        for i=1:nFirstBtnAls
            btn2.removeActionListener(btnAls(i));
        end
        
        if ~isempty(btn2)
            set(handle(btn2, 'CallbackProperties'), ...
                'ActionPerformedCallback', @(h,e)close(h));
        end
        if strcmpi('Cancel', char(btn2.getText))
            Gui.RegisterEscape(jd.getRootPane, btn2)
        end
    end
end

    function close(h)
        if isempty(h)
            finalAnsw='';
        else
            finalAnsw=char(h.getText);
        end
        try
            ok=feval(checkFnc, jd, finalAnsw);
            if ~ok
                return;
            end
        catch ex
            ex.getReport
        end
        jd.dispose;
    end

    function idx=indexOf(str)
        for i=1:jsa.length
            if isequal(str, char(jsa(i)))
                idx=i;
                return;
            end
        end
        idx=-1;
    end
    function rememberAnswer
        if ~isempty(rememberOnly) && ~isequal(answer, rememberOnly)
                msg(['<html><center>In this particular case only the<br>' ...
                    'answer <b>' rememberOnly '</b> can be remembered!!'...
                    '</center></html>'], 5, 'south east');
                return;
        end
        curFig=get(0, 'currentFigure');
        [~, ~, ~, ~, quadrant]=Gui.FindScreen(curFig);
        ra=RememberedAnswers;
        idx=indexOf(answer);
        pane.setInitialValue(jsa(idx));
        if strcmpi('west', quadrant{2})
            where='north east++';
        else
            where='north west++';
        end
        PopUp.Pane(pane, Title, where, javaWin, false, pauseSecs);
        ra.remember(rememberId, rememberCb, idx);
        raIfRemembered=ra;
    end
    
    function noTip
    BasicMap.Global.closeToolTip;
    end
end
