function r=randSequence(maxN,N)
% r=randSequence(maxN,N)
% produces a random  increasing sequence 
% that includes N integers 1:maxN
% (N<maxN)

if N>=maxN, error('N is not less than maxN!!!'); end;

fact=maxN/((maxN-N)*2);
ind=[];

while length(ind)~=N,
  ind=find(round(rand(1,maxN)*fact)>0);
end;
r=1:maxN;
r=r(ind);