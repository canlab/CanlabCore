function mea_visualise(plotmat,xaxis,yaxis,caxis)
% This program allows the visualisation of 3D images in separate subplots
%
% ..
%    chris summerfield 2003. summerfd@paradox.columbia.edu
% ..

if min(min(min(plotmat))) ~= max(max(max(plotmat))) 
    constrain_axis=1;
else
    constrain_axis=0;
end

if exist('xaxis','var')==0;
    xaxis=1:size(plotmat,3);
end

if exist('yaxis','var')==0;
    yaxis=1:size(plotmat,2);
end

if exist('caxis','var')==0;
    caxis=[min(min(min(plotmat))) max(max(max(plotmat)))];
end


numplot=size(plotmat,1);
scrsz = get(0,'ScreenSize');
figure('Position',[50 50 scrsz(3)/3 scrsz(4)/3])
for n=1:numplot;
    plotmat1=squeeze(plotmat(n,:,:));
    subplot('position',[0.1 1-((n/numplot)*0.95) 0.8 0.8/numplot])
    imagesc(xaxis,yaxis,plotmat1);
    if constrain_axis==1;
    set(gca, 'Clim',caxis);
    end
    colorbar
end;

%%%mea_visualise
