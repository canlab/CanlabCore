% Load the image
% =============================================
disp('Loading canonical T1 image.');drawnow
basename = 'single_subj_T1';
[array,hdr] = readim2(basename);
%origin = hdr.origin(1:3,1)';
clear D
for i = 1:size(array,3)
    E(:,:,i) = rot90(array(:,:,i));
end

startslice = 80;
endslice = 40;
index = 1;
figure
E(:,:,1:6) = [];

for i = startslice:-2:endslice
clf;cla
D = E;
D(:,:,i:end) = [];
disp('  Plotting surface');
p1 = patch(isosurface(D, 50),'FaceColor',[1,.75,.65],...,...
    'EdgeColor','none');
p2 = patch(isocaps(D, 50),'FaceColor','interp',...
    'EdgeColor','none');
view(300,20); axis tight; axis image; axis off
colormap(gray(100));
camlight left; camlight; camlight left; lighting gouraud
isonormals(D,p1)
if i == startslice
	myaxis = axis;
end
axis(myaxis)
drawnow

eval(['saveas(gcf,''colinhead-' num2str(index) ''',''jpg'')'])
index = index+1;
end
