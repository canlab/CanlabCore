%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function [yes, cancelled, raIfRemembered]=askYesOrNo(theMsg, title, ...
    where, defaultIsYes, rememberId, property)
if nargin<4
    defaultIsYes=true;
    if nargin<3
        where=[];
    if nargin<2
        title=[];
    end
    end
end
raIfRemembered=[];
if isempty(title)
    title= 'Please confirm...';
end
if nargin>2
    if isempty(where)
        where='center';
    end
    if ~isstruct(theMsg)
        m.msg=theMsg;
        m.where=where;
        if nargin>4 
            if ~isempty(rememberId)
                m.remember=rememberId;
            end
            if nargin>5
                m.property=property;
            end
        end
        theMsg=m;
    else
        theMsg.where=where;
    end
end
if defaultIsYes
    dflt='Yes';
else
    dflt='No';
end
if nargout>1
    [~,yes,cancelled, raIfRemembered]=questDlg(theMsg, title, 'Yes', 'No', ...
        'Cancel', dflt);
else
    [~,yes,cancelled, raIfRemembered]=questDlg(theMsg, title, 'Yes',...
        'No', dflt);
end
end