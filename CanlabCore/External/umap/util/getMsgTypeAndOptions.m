%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function [msgType, jsa, default, defaultIdx, isNumerics]=...
    getMsgTypeAndOptions(msgType, args)
defaultIdx=1;
N=length(args);
isNumerics={};
if N>0 && iscell(args{1})
    args=args{1};
    N=length(args);
end
if N<=1
    if N<1
        msgType=javax.swing.JOptionPane.INFORMATION_MESSAGE;
        msgTypeFound=true;
    else
        [msgType, msgTypeFound]=getMessageType(args{1}, msgType);
    end
    if msgTypeFound
        args={'Yes', 'No', 'Yes'};
        N=3;
    end
    start=1;
    nOptions=N;
end
if N>1
    [msgType, msgTypeFound]=getMessageType(args{1}, msgType);
    if ~msgTypeFound
        nOptions=N-1;
        start=1;
        lastOptionBeforeDefault=start+(nOptions-1);
    else
        start=2;
        if N==2
            nOptions=1;
            lastOptionBeforeDefault=0;
        else
            nOptions=N-2;
            lastOptionBeforeDefault=start+nOptions;
        end
    end
    for i=start:lastOptionBeforeDefault
        if strcmp(args{i}, args{end})
            defaultIdx=i-start+1;
            break;
        end
    end
    if defaultIdx>nOptions
        defaultIdx=1;
    end
end
if nOptions<=0
    jsa={};
    default='';
else
    jsa=javaArray('java.lang.String', nOptions);
    for i=1:nOptions
        arg=args{i+start-1};
        if isnumeric(arg)
            isNumerics{end+1}=true;
            jsa(i)=java.lang.String(num2str(arg));
        else
            isNumerics{end+1}=false;
            jsa(i)=java.lang.String(arg);
        end
    end
    default=jsa(defaultIdx);
end
end