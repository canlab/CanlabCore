%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

function [jd, pane]=msgBox(msgArgs, varargin)
app=BasicMap.Global;
app.closeToolTip;
msgType=javax.swing.JOptionPane.INFORMATION_MESSAGE;
title='Please note...';
where=[];
if nargin>1
    N=length(varargin);
    title=varargin{1};
    if N==2
        v=lower(varargin{2});
        if String.Contains(v,'west')||String.Contains(v, 'east') ||...
                String.Contains(v,'south') || String.Contains(v,'north')
            where=v;
        else
            msgType=getMessageType(v, msgType);
        end
    end
end
if ~isempty(where)
    [msgTxtOrObj, ~, ~, ~, ~, myIcon, javaWindow, ~, checkFnc, ...
        modal, pauseSecs]=decodeMsg(msgArgs, 0);
else
    [msgTxtOrObj, where, ~, ~, ~, myIcon, javaWindow, ~, checkFnc, ...
        modal, pauseSecs]=decodeMsg(msgArgs, 0);
end
if isstruct(msgArgs) && isfield(msgArgs, 'cancel_option') && ...
        msgArgs.cancel_option
    qu=javax.swing.JOptionPane.CANCEL_OPTION;
    pane=javaObjectEDT('javax.swing.JOptionPane', msgTxtOrObj, msgType, qu);
else
    
    pane=javaObjectEDT('javax.swing.JOptionPane', msgTxtOrObj, msgType);
end
if isempty(myIcon)
    if msgType==0
        myIcon='error.png';
    elseif msgType==1
        myIcon='facs.gif';
    elseif msgType==2
        myIcon='warning.png';
    else
        myIcon=[];
    end
end
if ~isempty(myIcon)
    pane.setIcon(Gui.Icon(myIcon));
end
suppressParent=isstruct(msgArgs) && isfield(msgArgs, 'suppressParent') ...
    && msgArgs.suppressParent;
if isstruct(msgArgs) && isfield(msgArgs, 'widthLimit') 
    jd=PopUp.Pane(pane, title, where, javaWindow, modal, pauseSecs, ...
        suppressParent, true, msgArgs.widthLimit);
else
    jd=PopUp.Pane(pane, title, where, javaWindow, modal, pauseSecs, ...
        suppressParent);
end
if ~modal && ~isempty(checkFnc)
    dlg=handle(jd, 'CallbackProperties');
    set(dlg, 'WindowClosingCallback', @(h,e)close());
    btn=javaObjectEDT('javax.swing.JButton');
    btnClass=btn.getClass;
    firstBtn=Gui.FindFirst(jd, btnClass, 'Ok');
    if isempty(firstBtn)
        firstBtn=Gui.FindFirst(jd, btnClass, PopUp.CLOSE_LABEL);
    end
    if ~isempty(firstBtn)
        set(handle(firstBtn, 'CallbackProperties'), ...
            'ActionPerformedCallback', @(h,e)close);
        if app.toolBarFactor>0
            txt=(char(firstBtn.getText));
            txt=strrep(txt, '<small>', app.smallStart);
            txt=strrep(txt, '</small>', app.smallEnd);
            firstBtn.setText(txt);
        end
    end
end

    function close
        try
            ok=feval(checkFnc, jd);
            if ~ok
                return;
            end
        catch ex
            ex.getReport
        end
        jd.dispose;
    end

end

