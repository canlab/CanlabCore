function str = scn_get_datetime(varargin)
% :Usage:
% ::
%
%     str = scn_get_datetime
%
% pass in 'ymd' to get the string in yyyy_mm_dd-HH_MM format, so that
% alphanumeric order will correspond to chronological order
%
% Returns a string with the date and time
% Useful for annotating data and output

if strmatch('ymd', varargin)
    datestring=datestr(date, 'yyyy-mm-dd');

else
    datestring=date;
end

    cval = fix(clock); 

    str = sprintf('%s_%02d:%02d', datestring, cval(4),  cval(5));

    str(str == ':') = '_';
end
