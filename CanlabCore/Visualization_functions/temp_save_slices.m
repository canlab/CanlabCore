img = readim2('homocor');
fig=figure;
set(fig,'DoubleBuffer','on');
set(gcf,'Color','k')
colormap(gray)
set(gca,'NextPlot','replace','Visible','off')
axis off
set(gcf,'Color','w')
mymax = mean(mean(mean(img))) + 8*mean(std(mean(img)));
	for i=1:size(img,3)
		imagesc(img(:,:,i),[0 mymax]); axis image;axis off;camzoom(1.3)
        saveas(gcf,['slice' num2str(i)],'jpg')
    end
    
    close
    