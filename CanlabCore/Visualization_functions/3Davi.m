% Makes an avi movie file called head3d[x].avi
% This script spirals up, right, and in 72 degrees
%
% notes: my indeo5 one wouldn't work in media player
%
% start with the image in the location you want to zoom in on
% but with no zoom.

mov = avifile(['head3d2.avi'],'Quality',100,'Compression','Indeo3','Fps',10);
H = gca;

% -------------------------------------------
% * set initial view
% -------------------------------------------
set(gcf,'Color','k')
axis off
[az,el] = view;
az = az - 36; el = el - 20;
view(az,el)

% -------------------------------------------
% * rotate and add frames
% -------------------------------------------
for i = 1:36
	az = az + 1;
	el = el + 1;
	view(az,el); 
    lightFollowView
	camzoom(1.03); drawnow
	mov = addframe(mov,H);
end

mov = close(mov);
