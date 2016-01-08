function Pout = scale_imgs_by_csf(hP)
% Takes a string matrix of image file names
% finds the mean and std of the CSF space
% specified in a mask (hard-coded)
% and standardizes images by these values
%
% :Usage:
% ::
%
%     Pout = scale_imgs_by_csf(hP)
%
% Writes SC* images (SCaled)
%
% assumes images are spatially normalized.
% uses a canonical CSF mask!
%
% ..
%    tor wager
% ..

mP = which('canonical_ventricles.img');	
[tmp,mP] = reslice_imgs(deblank(hP(1,:)),mP,0);

M = spm_general_hist(hP,mP,'canonical_ventricles');
disp(' ')
disp(['Means are: ' num2str(M(:,1)')])
disp(['Var    is: ' num2str(M(:,2)')])

for i = 1:size(hP,1)

	V = spm_vol(deblank(hP(i,:)));
	v = spm_read_vols(V);
	
	v = v - M(i,1);
	v = v ./ M(i,2);

	[d,f,e] = fileparts(V.fname);
	V.fname = [d filesep 'SC' f e];

	if i == 1, Pout = V.fname;, else, Pout = str2mat(Pout,V.fname);, end
	spm_write_vol(V,v);
	disp(['Written ' V.fname])

end
