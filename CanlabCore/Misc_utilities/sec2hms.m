function [day, hour, minute, second, display_string] = sec2hms(sec)
%SEC2HMS  Convert time from seconds (sec) into [days, hours, minutes, seconds]
%
%   [day, hour, minute, second, display_string] = SEC2HMS(SEC) converts the number of seconds in
%   SEC into days, hours, minutes and seconds.
%
% Example:
% t1 = clock;
% [day, hour, minute, second, display_string] = sec2hms(etime(clock, t1))
% disp(display_string)

day    = fix(sec/(3600*24)); % get number of days
sec    = sec - 3600*24*day;  % remove the days
hour   = fix(sec/3600);      % get number of hours
sec    = sec - 3600*hour;    % remove the hours
minute = fix(sec/60);        % get number of minutes
sec    = sec - 60*minute;    % remove the minutes
second = fix(sec);

% build display string

display_string = sprintf('%02d d %02d h %02d min %02d sec',day, hour, minute, second);

end
