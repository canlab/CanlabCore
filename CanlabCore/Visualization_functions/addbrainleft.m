function p = addbrainleft(varargin)
% :Usage:
% ::
%
%    handle = addbrainleft(enter arg to suppress lighting changes)
%
% quick function to add transparent brain surface to figure


Ps = which('surf_single_subj_T1_gray.mat'); %'c:\tor_scripts\3DheadUtility\surf_single_subj_T1_gray.mat';
Ps = which('surf_single_subj_grayL.mat');
%Ps = which('surf_brain_render_T1_preCarmack.mat'); 
load(Ps)
 p = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5], ...
  'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',1,'SpecularExponent',200);

if length(varargin) == 0
 lighting gouraud;camlight right
 axis image; myLight = camlight(0,0);set(myLight,'Tag','myLight');
 set(gcf, 'WindowButtonUpFcn', 'lightFollowView');lightFollowView
end

 set(p,'FaceAlpha',.3)
 
  drawnow
  axis vis3d
  %view(135,30)
  return
