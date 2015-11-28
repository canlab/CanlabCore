function [dat, dati] = union(dat1, dat2, outputname)
% Union and intersection masks for two image_vector objects
%
% :Usage:
% ::
%
%    [dat_union, dat_intersection] = union(dat1, dat2, outputname)
%
%     dat = union(dat1, dat2, outputname)
%     outputname = character array name for union image
%                   INCLUDE .img at the end.
%
% ..
%    NOTE: must now be in same space!
%    tor
% ..

isdiff = compare_space(dat1, dat2);
if isdiff == 1 || isdiff == 2
    error('Spaces are not the same!')
end


vdat1 = reconstruct_image(dat1);
vdat2 = reconstruct_image(dat2);

sz = size(vdat1);

vdat1 = vdat1(:);
vdat2 = vdat2(:);

both = vdat1 | vdat2;
oneandtwo = vdat1 & vdat2;

dat = dat1;
dat.volInfo.fname = 'REMOVED BY UNION';
dat.volInfo.descrip = 'UNION mask';
dat.volInfo.image_indx = both;
dat.volInfo.nvox = length(both);
dat.volInfo.wh_inmask = find(both);
dat.volInfo.n_inmask = sum(both);

[x, y, z] = ind2sub(sz, 1:length(vdat1));
xyz = [x' y' z'];
dat.volInfo.xyzlist = xyz(both, :);

if dat.volInfo(1).n_inmask < 50000
    dat.volInfo(1).cluster = spm_clusters(dat.volInfo(1).xyzlist')';
else
    dat.volInfo(1).cluster = ones(dat.volInfo(1).n_inmask, 1);
end

if isa(dat1, 'statistic_image')
    
    [p, ste, sig] = deal(zeros(size(both)));
    p(dat1.volInfo.wh_inmask, :) = dat1.p(:, 1);
    ste(dat1.volInfo.wh_inmask, :) = dat1.ste(:, 1);
    sig(dat1.volInfo.wh_inmask, :) = dat1.sig(:, 1);

end

if isa(dat2, 'statistic_image')
    p(:, 2) = 0;
    ste(:, 2) = 0;
    sig(:, 2) = 0;
    p(dat1.volInfo.wh_inmask, 2) = dat2.p(:, 1);
    ste(dat1.volInfo.wh_inmask, 2) = dat2.ste(:, 1);
    sig(dat1.volInfo.wh_inmask, 2) = dat2.sig(:, 1);

end

if isa(dat1, 'statistic_image') || isa(dat2, 'statistic_image')
    p = p(both, :);
    ste = ste(both, :);
    sig = sig(both, :);
    
end

% could average here...
dat.dat = single(both(both));

if isa(dat1, 'statistic_image')

    dat.p = min(p, [], 2);
    dat.ste = min(ste, [], 2);
    dat.sig = max(sig, [], 2);
end


% intersection
dati = dat;
dati.volInfo.fname = 'REMOVED BY UNION';
dati.volInfo.descrip = 'INTERSECTION mask';

dati.dat = single(oneandtwo(both));

if isa(dat1, 'statistic_image')

    dati.p = max(p, [], 2);
    dati.ste = max(ste, [], 2);
    dati.sig = min(sig, [], 2);
end

if nargin > 2
    
    dat.fullpath = fullfile(pwd, outputname);
    write(dat)
    
end


end


