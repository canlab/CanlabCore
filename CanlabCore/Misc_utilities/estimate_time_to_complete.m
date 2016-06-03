function [day, hour, minute, second, display_string] = estimate_time_to_complete(starttime, n_done, n_total, varargin)
% Estimate time to complete a job given start time and how much done so far
%
% [day, hour, minute, second, display_string] = estimate_time_to_complete(starttime, n_done, n_total, [verbose])
% 
% Example:
%
% starttime = clock;
% n_done = 10;
% n_total = 1000;
% [day, hour, minute, second, display_string] = estimate_time_to_complete(starttime, n_done, n_total, 1);
% print and erase:
% fprintf(display_string); bspace = repmat('\b', 1, length(display_string)); fprintf(bspace)
%
% 2016 Tor Wager

doverbose = 0;

if length(varargin) > 0
    doverbose = varargin{1};
end

time_elapsed = etime(clock, starttime);

time_per_unit = time_elapsed ./ n_done;

time_remaining = time_per_unit .* (n_total - n_done);

[day, hour, minute, second, display_string] = sec2hms(time_remaining);

if doverbose, fprintf(display_string); end

end


