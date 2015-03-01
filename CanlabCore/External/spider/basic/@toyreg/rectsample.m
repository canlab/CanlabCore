function X=rectsample(N)


X=rand(N,2);
X=X+repmat([-0.5,-0.5],N,1);