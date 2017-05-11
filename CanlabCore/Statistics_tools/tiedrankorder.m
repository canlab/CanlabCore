function [r] = tiedrankorder(x,tiemethod,epsx)
%TR Local tiedrank function to compute results for one column

if nargin<2, tiemethod = 'average'; end
if nargin<3, epsx = zeros(size(x)); end


% Sort, then leave the NaNs (which are sorted to the end) alone
[sx, rowidx] = sort(x(:));
epsx = epsx(rowidx);
epsx = epsx(:);
numNaNs = sum(isnan(x));
xLen = numel(x) - numNaNs;


% Use ranks counting from low end
ranks = [1:xLen NaN(1,numNaNs)]';

if isa(x,'single')
   ranks = single(ranks);
end

% Adjust for ties.  Avoid using diff(sx) here in case there are infs.
ties = sx(1:xLen-1)+epsx(1:xLen-1) >= sx(2:xLen)-epsx(2:xLen);
tieloc = [find(ties); xLen+2];
maxTies = numel(tieloc);

tiecount = 1;
while (tiecount < maxTies)
    tiestart = tieloc(tiecount);
    ntied = 2;
    while(tieloc(tiecount+1) == tieloc(tiecount)+1)
        tiecount = tiecount+1;
        ntied = ntied+1;
    end
    
    switch tiemethod
        case 'average'
            % Compute mean of tied ranks
            ranks(tiestart:tiestart+ntied-1) = ...
                sum(ranks(tiestart:tiestart+ntied-1)) / ntied;
        case 'first'
            % Use min of tied ranks
            ranks(tiestart:tiestart+ntied-1) = ranks(tiestart:tiestart+ntied-1);
            
        case 'min'
            % Use min of tied ranks
            ranks(tiestart:tiestart+ntied-1) = ranks(tiestart);                
    
        case 'max'
            % Use max of tied ranks
            ranks(tiestart:tiestart+ntied-1) = ranks(tiestart+ntied-1);
    end
    tiecount = tiecount + 1;
end

% Broadcast the ranks back out, including NaN where required.
r(rowidx) = ranks;
r = reshape(r,size(x));
