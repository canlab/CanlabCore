function pv = constrain_toverlap(pv,t,numframes)
% pv = constrain_toverlap(pv,t,numframes)
% constrains paramvec to be at least t elements apart

pv2 = cat(2,pv{:});

tmp=[]; 
for i=1:length(pv), tmp=[tmp i*ones(size(pv{i}))];,end

pv2 = [pv2; tmp];
pv2 = sortrows(pv2',1);

alltons = pv2(:,1);
tmp = [Inf; diff(pv2(:,1))];

wh = find(tmp < t);
        
while ~isempty(wh) & alltons(wh(1)) < numframes
    alltons(wh(1):end) = alltons(wh(1):end) + (t - tmp(wh(1)));
    tmp = [Inf; diff(alltons)];
    wh = find(tmp < t);
end

pv2 = [alltons pv2]; pv2(pv2(:,1) > numframes,:) = [];

for i=1:max(pv2(:,3))
    pv{i} = pv2(pv2(:,3)==i,1)';
end

return
