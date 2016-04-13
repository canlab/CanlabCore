function [cnt, tot, lenmat] = cnt_runs(tc)
% Counts number of contiguous 1s in a timeseries
%
% :Usage:
% ::
%
%     function [cnt, tot, lenmat] = cnt_runs(tc)
%

tot = sum(tc);
[a1 b1] = max(diff(tc));
[a2 b2] = min(diff(tc));
cnt = 0;
T = length(tc);
lenmat = zeros(T,1);

if (b1 == 1),    
    b1 =0;    
end;

while (b2 > 1),
    if (b1 == 1 | b1>b2),
        b1 =0;
    end;

    

    if (b1 < b2),
        len = b2 - b1;
        tc(1:b2) = 0;
        cnt = cnt +1;
        lenmat(cnt) = len;

    end;

    

    [a1 b1] = max(diff(tc));
    [a2 b2] = min(diff(tc));

end;



if (tc(end) == 1),

    cnt = cnt+1;

    len = length(tc) - b1;

    lenmat(cnt) = len;

end

return
