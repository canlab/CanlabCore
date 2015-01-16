function tor_spm_mean_ui(P,varargin)
%tor_spm_mean_ui(Pinputnames,[outputname])
% promts for a series of images and averages them
% FORMAT spm_mean_ui
%_______________________________________________________________________
%
% spm_mean_ui simply averages a set of images to produce a mean image
% that is written as type int16 to "mean.img" (in the current directory).
%
% The images must have the same dimensions, orientations (as defined by
% the Origin header field or any associated *.mat files), and the same
% voxel sizes.
%
% This is not a "softmean" - zero voxels are treated as zero.
%_______________________________________________________________________
% @(#)spm_mean_ui.m	2.4 John Ashburner, Andrew Holmes 98/10/21
SCCSid = '2.4';

% Tor modified: select output file name
%-----------------------------------------------------------------------
if nargin > 1, 
	outname = varargin{1};
else
	outname = 'mean.img';
end



%-Say hello
%-----------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SCCSid);


%-Select images & check dimensions, orientations and voxel sizes
%-----------------------------------------------------------------------
if nargin == 0
	fprintf('\t...select files')
	P = spm_get(Inf,'.img','Select images to be averaged');
   fprintf(' ...mapping & checking files')
end

Vi = spm_vol(P);

n  = prod(size(Vi));
if n==0, fprintf('\t%s : no images selected\n\n',mfilename), return, end

if n>1 & any(any(diff(cat(1,Vi.dim),1,1),1)&[1,1,1,0])
	error('images don''t all have same dimensions'), end
if any(any(any(diff(cat(3,Vi.mat),1,3),3)))
	error('images don''t all have same orientation & voxel size'), end


%-Compute mean and write headers etc.
%-----------------------------------------------------------------------
fprintf(' ...computing')
Vo = struct(	'fname',	outname,...
		'dim',		[Vi(1).dim(1:3),4],...
		'mat',		Vi(1).mat,...
		'pinfo',	[1.0,0,0]',...
		'descrip',	'spm - mean image');

%-Adjust scalefactors by 1/n to effect mean by summing
for i=1:prod(size(Vi))
	Vi(i).pinfo(1:2,:) = Vi(i).pinfo(1:2,:)/n; end;

%-Write basic header
spm_create_image(Vo);

%-Use spm_add to do the donkey work
Vo.pinfo(1,1) = spm_add(Vi,Vo);

%-Write header (complete with scaling information)
spm_create_image(Vo);


%-End - report back
%-----------------------------------------------------------------------
fprintf(' ...done\n')
fprintf('\tMean image written to file ''%s'' in current directory\n\n',Vo.fname)
