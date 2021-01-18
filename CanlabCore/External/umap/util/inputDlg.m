%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University
%   License: BSD 3 clause
%   


function [answer, cancelled]=inputDlg(msg,title,varargin)
if nargin==1
    title='Input required';
end
[msgType, jsa,default,~,isNumerics]=getMsgTypeAndOptions(...
    javax.swing.JOptionPane.QUESTION_MESSAGE, varargin);
[msg, where, property, properties, default, myIcon, javaWin,~,~,modal]...
    =decodeMsg(msg, default);

if ~isempty(jsa)
    inputValue=char(jsa(1));
else
    inputValue='';
end
if msgType==0
    myIcon='error.png';
elseif msgType==1
    myIcon = 'facs.gif';
elseif msgType==2
    myIcon='warning.png';
else
    myIcon='question.png';
end
pane=javaObjectEDT('javax.swing.JOptionPane', msg, msgType);
pane.setWantsInput(true);
pane.setInitialSelectionValue(inputValue);
pane.selectInitialValue();
pane.setIcon(Gui.Icon(myIcon));
pane.setOptionType(javax.swing.JOptionPane.OK_CANCEL_OPTION);
PopUp.Pane(pane, title,where, javaWin, modal);
answer=pane.getInputValue;
cancelled=strcmp(answer,'uninitializedValue');
if cancelled
    answer='';
end
end
