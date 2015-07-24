% start in main expt directory, above ind img/anatomy directories
% run this script
%
% makes movies of all slices for all subjects in all_anatomy subdir

D = dir('all_anatomy/*img');

for i = 1:length(D)

img = readim2(['all_anatomy' filesep D(i).name]);

fig=figure;
set(fig,'DoubleBuffer','on');
set(gcf,'Color','k')
colormap(gray)
set(gca,'NextPlot','replace','Visible','off')
axis off

try
	
	mov = avifile([D(i).name(1:8) '.avi'],'Quality',80,'Compression','Indeo3','Fps',3)
     
	for i=1:size(img,3)
		imagesc(img(:,:,i)); axis image;
		drawnow
       	% set(h,'EraseMode','xor');
       	F = getframe(gca);
       	mov = addframe(mov,F);
	end
       
	mov = close(mov);
	
catch
	
	mov = close(mov);

end




end
