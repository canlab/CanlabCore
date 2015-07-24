%for mysub = {'RH' 'AB' 'DC' 'PD' 'PS' 'KH' 'MS' 'NS' 'RM' 'HP' 'LB'}
%for mysub = {'RH' 'TD'}

clear mysub
a = pwd;
a = a(end-1:end)
mysub{1} = ['biman5_' a '_nocompress'];
%img = readim2([mysub{1} '3dvol']);

if ~(exist('homocor.img') == 2), homocorr('t1spgr.img'), end
img = readim2(['homocor']);

fig=figure;
set(fig,'DoubleBuffer','on');
set(gcf,'Color','k')
colormap(gray)
set(gca,'NextPlot','replace','Visible','off')
axis off
set(gcf,'Color','w')

%try
	clear mov
	mov = avifile([mysub{1} '.avi'],'Quality',100,'Compression','none','Fps',3)
    mymax = mean(mean(mean(img))) + 8*mean(std(mean(img)));
	for i=1:size(img,3)
		imagesc(img(:,:,i),[0 mymax]); axis image;axis off;camzoom(1.3)
		drawnow
       	% set(h,'EraseMode','xor');
       	F = getframe(gca);
       	mov = addframe(mov,F);
	end
       
	mov = close(mov);
	
    %catch
	
	%mov = close(mov);

    %end




%end
