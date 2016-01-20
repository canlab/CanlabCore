function [anybad,wh] = iimg_check_volinfo(maskInfo,imgInfo)
% Checks a series of image .mat files and dims against a reference
% (maskInfo)
%
% :Usage:
% ::
%
%     anybad = iimg_check_volinfo(maskInfo,imgInfo)
%
% :Inputs:
%
%   **maskInfo** and **volInfo :**
%        are spm-style volume info structures
%
% see spm_vol.m
 
wh = [];

n = length(imgInfo);
notok = zeros(1,n);

tol = .01;

for i=1:n
    chk = abs(maskInfo.mat - imgInfo(i).mat) > tol; 
    
    chk = any(diag(chk(1:3, 1:3)));
    
% %     chk = chk(:);
% %     chk = chk(1:end-1);             % eliminate SPM scale factor and translation
% %     chk1 = any(chk);
    
    chk2 = any(maskInfo.dim(1:3) - imgInfo(i).dim(1:3));
    
    notok(i) = chk | chk2;
end

anybad = any(notok);

if anybad
    wh = find(notok);

    disp('The following images'' mat files or dims differed from the first:')
    disp(num2str(wh));

    disp('First mat:');
    disp(maskInfo.mat);
    disp(maskInfo.dim);
    
    disp('Bad mats:');
    for i = wh
        disp(imgInfo(i).fname);
        disp(imgInfo(i).mat);
        disp(imgInfo(i).dim);
    end
end

return
