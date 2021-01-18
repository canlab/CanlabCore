%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function answer=ask(question, defaultAnswer, title, where, cancelToo)

if nargin<5
    cancelToo=true;
    if nargin<4
        where='center';
        if nargin<3
            title=[];
            if nargin<2
                defaultAnswer=1;
            end
        end
    end
end
if isempty(title)
    title='Please confirm....';
end
pu=PopUp(question, where, title, false, [],[], true);
pu.addYesNo(cancelToo, defaultAnswer);
pu.dlg.setAlwaysOnTop(true);
pu.dlg.setModal(true);
Gui.SetJavaVisible(pu.dlg);
answer=pu.answer;
