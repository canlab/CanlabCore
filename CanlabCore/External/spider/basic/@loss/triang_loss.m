function dat = triang_loss(algo,dat)

% so called triangular loss function that calculates the mean error in
% degrees for a circular classification problem where 0 degree is
% equivalent to 180 degrees


[x y]=get_xy(dat);

labels = [0 1 2 3 4 3 2 1];
labels = labels *22.5 ;
L=[];

for i = (0:7)
    L = [L circshift(labels',i)];
end

if size(y,2)>1
    %convert back into angles from basis-vectors
    y = (y+1)/2;
    x = (x+1)/2;
    y = y*[ 1 2 3 4 5 6 7 8]';
    x = x*[ 1 2 3 4 5 6 7 8]';
end

loss = trace( L(x, y)) / length(y);

dat=data([get_name(dat) ' -> triang_loss='  num2str(loss,4) ],[],loss);
