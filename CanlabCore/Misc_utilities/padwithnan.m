function [Anew Bnew] = padwithnan(A, B, dim)
% returns the two input arrays, with the smaller padded to the size of the
% larger in the particular dimension(s) with NaNs
%
% :Usage:
% ::
%
%     [Anew Bnew] = padwithnan(A, B, dim)
%

    if(nargin ~= 3)
        error('padwithnan: incorrect number of arguments');
    end


    diff = size(A,dim) - size(B,dim);

    if diff > 0
        Bnew = pwn_expand(B, dim, diff);
        Anew = A;
    elseif diff < 0
        Anew = pwn_expand(A, dim, diff);
        Bnew = B;
    else
        Anew = A;
        Bnew = B;
    end
end

function outarray = pwn_expand(inarray, dim, diff)
    sznandims = size(inarray);
    sznandims(dim) = abs(diff);
    outarray = cat(dim, inarray, NaN(sznandims));
end
