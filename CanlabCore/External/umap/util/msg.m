%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function [jd, pane]=msg(msssageTextOrObj, pauseSecs, where, title, ...
    icon, checkFnc, suppressParent)
if nargin<7
    suppressParent=false;
    if nargin<6
        checkFnc=[];
        if nargin<5
            icon='facs.gif';
            if nargin<4
                title='Note...';
                if nargin<3
                    where='center';
                    if nargin<2
                        pauseSecs=15;
                    end
                end
            end
        end
    end
end
if isstruct(msssageTextOrObj)
    st=msssageTextOrObj;
    if ~isfield(st, 'pauseSecs')
        st=setfield(st, 'pauseSecs', pauseSecs);
    end
    if ~isfield(st, 'modal')
        st=setfield(st, 'modal', false);
    end
    if ~isfield(st, 'suppressParent')
        st=setfield(st, 'suppressParent', suppressParent);
    end
    if ~isfield(st, 'where')
        st=setfield(st, 'where', where);
    end
    if ~isfield(st, 'icon')
        st=setfield(st, 'icon', icon);
    end
    if ~isfield(st, 'checkFnc')
        st=setfield(st, 'checkFnc', checkFnc);
    end
    [jd, pane]=msgBox(st, title);
    return;
end

if nargin<5
    [jd, pane]=msgBox(struct('msg', msssageTextOrObj, 'modal', false, ...
        'pauseSecs', pauseSecs, 'where', where, ...
        'suppressParent', suppressParent), title);
else
    if ~isempty(icon)
        [jd, pane]=msgBox(struct('msg', msssageTextOrObj, 'modal', false, ...
            'pauseSecs', pauseSecs, 'where', where, 'icon', icon, ...
            'checkFnc', checkFnc,'suppressParent', suppressParent), title);
    else
        [jd, pane]=msgBox(struct('msg', msssageTextOrObj, 'modal', false, ...
            'pauseSecs', pauseSecs, 'where', where, 'icon', icon, ...
            'checkFnc', checkFnc, ...
            'suppressParent', suppressParent), title, 'plain');
    end
end
end
