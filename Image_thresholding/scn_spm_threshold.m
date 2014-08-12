
% img should be residual images, I believe.

VM = spm_vol(mask);

[FWHM,VRpv] = spm_est_smoothness(img, mask);
R = spm_resels_vol(VM, FWHM)';

volInfo = iimg_read_img(mask, 2);

u = spm_uc(.05, [1 20], 'T', R, 1, volInfo.n_inmask);

[P p Em En EN] = spm_P(c,k,Z,df,STAT,R,n,S)

k = 20
u = .01
[P Pn Em En EN] = spm_P(1, k, u, [1 20], 'T', R, 1, volInfo.n_inmask);
