function [CP2,CNTmat,TOTmat,LENmat] = zero_crossing(M,T,XX,ZIU,ZIL,Z,mu)
% Estimate out of control points and change-point using zero-crossing
% method.
%
% :Usage:
% ::
%
%     [CP2,CNTmat,TOTmat,LENmat] = zero_crossing(M,T,XX,ZIU,ZIL,Z,mu)
%
% Used in ewma5.m  See ewma5 for description of variables.


CP2 = NaN .* zeros(M,1); % Estimation via Zero-crossings
CNTmat  = zeros(M,1);
TOTmat = zeros(M,1);
LENmat = zeros(M,T);

wh = all(XX' - repmat(mean(XX'),size(XX,2),1) <= eps) | any(isnan(XX),2)';
wh = find(~wh);
    
ooc = ZIU - ZIL;    % pos if occ +, neg if ooc -, zero otherwise
for d=wh,

    % Calculate last 'zero-crossing'

    amiout = find(ooc(d,:)); 
    if isempty(amiout)
        amiout = 0;
    else
        amiout = amiout(1);
        amiout = ooc(d,amiout);  % first ooc tp, -1 or 1
    end
    
    switch amiout
        
        case 1
            dat = Z(d,:)-mu(d);
            tmp = ZIU(d,:);
            [a,b] = max((conv([1 1 1],tmp) == 3));  % must be OOC 3 consecutive points
            [a,zc] = max((dat(1:b) < 0).*(1:b));
            CP2(d) = zc;

        case -1
            dat = Z(d,:)-mu(d);
            tmp = ZIL(d,:);
            [a,b] = max((conv([1 1 1],tmp) == 3));  % must be OOC 3 consecutive points
            [a,zc] = max((dat(1:b) > 0).*(1:b));
            CP2(d) = zc;
            
        case 0
            CP2(d) = NaN;
            tmp = 0 .* ZIL(d,:);
        otherwise
            error('This should never happen'),keyboard
    end

    % Calculate number of out-of-control runs and their widths

    [cnt tot len] = cnt_runs(tmp);
    CNTmat(d) = cnt;
    TOTmat(d) = tot;
    LENmat(d,:) = len';

end;

CP2(find(CP2 == 1)) = NaN;

return
