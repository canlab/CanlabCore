function events=balancedRnd(n,nVals,nOnes,overlapYes)

% r=balancedRnd(n,nVals,nOnes[,overlapYes])
% creates a n X nVals matrix of 0 and 1's with equal fraction 0<nOnes<=1 for all vals
% simply iterates over many rand values
% overlapYes-alows overlaps
% assuming that one needs as many null events as any other

if nargin<4, overlapYes=0; end;

n0=n;
n=floor(n/nVals)*nVals;
r=zeros(1,n);
events=zeros(nVals,n);


if ~overlapYes,
   nSub=n/nVals;
	nOnes=round(nSub*nOnes);
   nind=find(r==0);
   for i=1:nVals, 
      ind=randSequence(length(nind),nOnes);
      if ~isempty(ind),
         r(nind(ind))=i;
         nind=find(r==0);
      end;
   end;
   for i=1:nVals,
      if ~isempty(ind),
	      ind=find(r==i);
         events(i,ind)=1;
      end;
   end;
   
else % if overlaps are alowed
	nOnes=round(n*nOnes);
   for i=1:nVals, 
      ind=randSequence(n,nOnes);
      if ~isempty(ind) events(i,ind)=1; end;
   end; 
end;


if n0>n, events=[events zeros(size(events,1),n0-n)]; end;

