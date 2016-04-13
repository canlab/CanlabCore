function [hh1,hh2,hh3,hl,a1,a2,a3] = make_figure_into_orthviews()
% :Usage:
% ::
%
%    [hh1,hh2,hh3,hl,a1,a2,a3] = make_figure_into_orthviews
%
% Copies a surface rendering or glass brain into three separate view
% panels, one saggital, one axial, and one coronal
%
% returns handles to objects in each view and hl light handles
%
% :Examples:
% ::
%
%    [hh1,hh2,hh3,hl,a1,a2,a3] = make_figure_into_orthviews;
%    axes(a1)
%    text(-55,60,70,'L','FontSize',24);text(55,60,70,'R','FontSize',24);
%    axes(a3)
%    text(-55,60,70,'L','FontSize',24);text(55,60,70,'R','FontSize',24);
%
% ..
%    tor wager
% ..

fo = gcf;
f = figure('Color','w'); a1=subplot(1,3,1); set(gca,'FontSize',16);
a2=subplot(1,3,2); set(gca,'FontSize',16);
a3=subplot(1,3,3);set(gca,'FontSize',16);

figure(fo);
o1=get(gca,'Children');
hh1=copyobj(o1,a1);
hh2=copyobj(o1,a2);
hh3=copyobj(o1,a3);

axes(a1); , axis image,
axes(a2); view(270,0), axis image
hl = lightangle(270,0);
axes(a3); view(0,0), axis image
hl(2) = lightangle(0,0);

axes(a1)
text(-55,60,70,'L','FontSize',24);text(55,60,70,'R','FontSize',24);
axes(a3)
text(-55,60,70,'L','FontSize',24);text(55,60,70,'R','FontSize',24);

return
