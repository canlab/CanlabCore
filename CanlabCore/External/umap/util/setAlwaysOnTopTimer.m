%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

function setAlwaysOnTopTimer(jd, time, focus, ignoreIfPc)
if nargin<4
    ignoreIfPc=true;
    if nargin<3
        focus=[];
        if nargin<2
            time=1;
        end
    end
end
if isempty(jd)
    return;
end
jd.toFront;
if ispc && ignoreIfPc
    return;
end
%disp(['Setting window on top for ' num2str(time) ' seconds']);
javaMethodEDT( 'setAlwaysOnTop', jd, true);
tmr=timer;
tmr.StartDelay=time;
tmr.TimerFcn=@(h,e)dismiss;
start(tmr);

    function dismiss
        javaMethodEDT('setAlwaysOnTop', jd, false);
        if ~isempty(focus)
            javaMethodEDT('requestFocus',focus);
        end
    end
end
