function [H,mi,H2] = hist2(A,B,res,varargin)
% 2-D histogram with res bins
%
% :Usage:
% ::
%
%     [H,mi,H2] = hist2(A,B,res,[plot])
%
% A and B can be 3D, as in image volumes
% mi is mutual information, a la spm_mireg.m
%
% ..
%    tor wager, may 6, 2003
% ..

A = A(:);
B = B(:);

r1 = abs((max(A)-min(A))./(res-1));
[h1,x] = histc(A,min(A):r1:max(A)-r1);

r2 = abs((max(B)-min(B))./(res-1));
[h2,y] = histc(B,min(B):r2:max(B)-r2);

% for i = 1:res
%     tmp = x == i;
%     for j = 1:res
%         H(i,j) = sum(tmp & y==j);
%     end
% end

for i = 1:res
    tmpx{i} = find(x == i);
    tmpy{i} = find(y == i);
end
for i = 1:res
    for j = 1:res
        H(i,j) = length(intersect(tmpx{i}, tmpy{j}));
    end
end

% mutual information, from that given in spm_mireg.m
H2  = H/(sum(H(:))+eps);
s1 = sum(H2,1);
s2 = sum(H2,2);
H2  = H2.*log2((H2+eps)./(s2*s1+eps));
mi = sum(H2(:));

if length(varargin) > 0, figure;imagesc(H), xlabel('Image 2'),ylabel('Image 1'),end

return
