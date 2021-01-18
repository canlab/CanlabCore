%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function [type, found]=getMessageType(arg, dflt)
if nargin<2
        dflt=javax.swing.JOptionPane.PLAIN_MESSAGE;
end
found=true;
if strcmpi(arg, 'error')
        type=javax.swing.JOptionPane.ERROR_MESSAGE;
    elseif strcmpi(arg, 'question')
        type=javax.swing.JOptionPane.QUESTION_MESSAGE;
    elseif strcmpi(arg, 'warning') || strcmpi(arg, 'warn');
        type=javax.swing.JOptionPane.WARNING_MESSAGE;
    elseif strcmpi(arg, 'information')
        type=javax.swing.JOptionPane.INFORMATION_MESSAGE;
    elseif strcmpi(arg, 'plain')
        type=javax.swing.JOptionPane.PLAIN_MESSAGE;
else
        found=false;
        type=dflt;
end
end