function p = add_surface(surf_mat)
    Ps = which(surf_mat); %'c:\tor_scripts\3DheadUtility\surf_single_subj_T1_gray.mat';
    if isempty(Ps), disp(['I need the file: ' surf_mat]); return; end

    %Ps = which('surf_single_subj_grayR.mat');
    %Ps = which('surf_brain_render_T1_preCarmack.mat');
    load(Ps)
    p = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5], ...
        'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',1,'SpecularExponent',200);

    set(p,'FaceAlpha',.3)
end
