function [anybad] = check_spm_matfiles(P)
% Check a list of SPM-style images to see if they all have the same
% space, dims.
%
% :Usage:
% ::
%
%     [anybad] = check_spm_matfiles(P)
%

V=spm_vol(P);
for i=2:length(V),
    chk = V(1).mat - V(i).mat; chk = chk(:);
    chk = chk(1:end-1);             % eliminate SPM scale factor
    chk = any(chk);
    notok(i) = chk;
end

if any(notok)
    wh = find(notok);

    disp(['The following images'' mat files differend from the first:'])
    disp(num2str(wh));

    disp(['First mat:']);
    disp(V(1).mat);
    
    disp(['Bad mats:']);
    for i = wh
        disp(V(i).fname);
        disp(V(i).mat);
    end
end

anybad = any(notok);
return
