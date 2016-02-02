function [hh,h2]=plot_ellipse(x,y,theta,a,b)
% PLOT_ELLIPSE
%
% :Usage:
% ::
%
%    [h,h2]=plot_ellipse(x,y,theta,a,b)
%
% This routine plots an ellipse with centre (x,y), axis lengths a,b
% with major axis at an angle of theta radians from the horizontal.
%
% ..
%    Author: P. Fieguth
%            Jan. 98
%
%    http://ocho.uwaterloo.ca/~pfieguth/Teaching/372/plot_ellipse.m
% ..


np = 100;
ang = [0:np]*2*pi/np;
pts = [x;y]*ones(size(ang)) + [cos(theta) -sin(theta); sin(theta) cos(theta)]*[cos(ang)*a; sin(ang)*b];
hh = plot( pts(1,:), pts(2,:), 'LineWidth',3,'Color','k'); hold on;
h2 = fill(pts(1,:), pts(2,:),'r','FaceAlpha',.2);

return

