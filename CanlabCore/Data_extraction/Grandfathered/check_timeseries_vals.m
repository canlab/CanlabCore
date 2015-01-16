function chk = check_timeseries_vals(V,dat,coords)
% function chk = check_timeseries_vals(V,dat,coords)
%
% Checks a random subset of up to 5 images against extracted data values
% to make sure the right data is being extracted.
%
% V is memory-mapped volumes (see spm_vol)
% dat is extracted data
% coords is voxel coordinates, n rows x 3 columns
%
% tor wager

n = min(5,length(V));

% get random set of n images, up to 5
wh = randperm(length(V));
wh = wh(1:n);

% select these rows in V and dat
% dat is [images, coordinates]
V = V(wh);
dat = dat(wh,:);

% get random set of nn coordinates, up to 5

nc = size(coords,1);
nn = min(5,nc);
whc = randperm(nc);
whc = wh(1:nn);
coords = coords(whc,:);

% select these columns in dat
dat = dat(:,whc);

v = spm_read_vols(V);
for i = 1:nn    % for each coordinate
    dat2(:,i) = squeeze(v(coords(i,1),coords(i,2),coords(i,3),:));
end

chk = dat - dat2;
chk = any(chk(:));

if chk, warning('Problem with timeseries!! Extracted data do not match expected values.');,end

return
    