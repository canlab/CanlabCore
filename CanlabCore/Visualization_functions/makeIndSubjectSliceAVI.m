% start in main expt directory, above ind img/anatomy directories
% run this script
%
% makes movies of all slices for all subjects in EXPT.subjects (cell array)

for mysub = EXPT.subjects

img = readim2([mysub{1} filesep 'anatomy' filesep 'het1overlay']);

fig=figure;
set(fig,'DoubleBuffer','on');
set(gcf,'Color','k')
colormap(gray)
set(gca,'NextPlot','replace','Visible','off')
axis off

try
	
	mov = avifile([mysub{1} '.avi'],'Quality',80,'Compression','Indeo3','Fps',3)
     
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
