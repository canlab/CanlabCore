
% input in seconds

TR = 2;     % 1 = leave in seconds; 2, downsample by 2, etc.

len = 200;  % length in s
tp = 5;    % estimates 20 time points

scanspersess = [40 40 10];     % images per session; length is num sessions

[X,delta,delta_hires,hrf] = onsets2delta(ons,TR,len);

[DX,sf] = tor_make_deconv_mtx3(delta,tp,1,0,1,0,scanspersess);


% for phys data, .01 s -> 1s
Y = RESAMPLE(X,1,100);

