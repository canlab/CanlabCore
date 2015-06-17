function o = X2onsets(X,TR,numcond)

[delta, wd] = DX_find_delta(X);

for i = 1:size(delta,2)
    ons(:,i) = find(delta(:,i));
end

ons = ons .* TR;

ons = ons - ons(1);


o = ons';
o = o(:);
o = diff(o);
o = [0;o];
o = reshape(o,numcond,length(o)./numcond);