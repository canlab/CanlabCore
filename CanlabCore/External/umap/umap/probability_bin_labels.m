%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Funded by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function y=probability_bin_labels(ptrs, labels, N)
assert(max(ptrs)<=N, 'Expect bin pointer can be higher than bin');
if size(ptrs,2)==1
    assert(size(ptrs,1)>N, 'Expect more bin pointers than bins');
elseif size(ptrs,1)==1
    assert(size(ptrs,2)>N, 'Expect more bin pointers than bins');
else
    assert(false, 'Expect 1 column for bin pointers');
end
y=zeros(1,N);
y(ptrs)=labels;
end