function P3 = reslice_seg_anatomy()

disp(['This script reslices all segmented n*het1*img anatomicals into space of functionals.'])
disp(['Output images are thresholded (1 or 0) to be mask images.'])

mdir = spm_get(-1,'*','Enter Main Directory for Study',pwd);
mypwd = pwd;
eval(['cd ' mdir])

!ls
dwcard = input('Enter wildcard for all individual subjects:','s');

P = get_filename(dwcard,['anatomy' filesep 'n*het1*img']);
disp(['reslicing these images:']), P

d = dir(dwcard);
P2 = get_filename([d(1).name '*'],'scan1/sn*ra*0001.img');
Pt = P2(1,:);

disp(['Reslicing to space of: ' Pt])

reslice_imgs(Pt,P,1);

P3 = get_filename(dwcard,['anatomy' filesep 'rn*het1*img']);

disp('Output in:')
P3

return




