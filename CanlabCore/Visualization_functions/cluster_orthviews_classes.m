function colors2 = cluster_orthviews_classes(cl,classes,overlay,myview,domontage, varargin)
% :Usage:
% ::
%
%    colors2 = cluster_orthviews_classes(cl,classes,overlay,myview,domontage, [orthview axis], [colors cell])
%
% Makes montage and cluster orthviews of classes of clusters in different
% colors (hard coded colors right now).
%
% :Inputs:
%
%   **cl:**
%        clusters structure with regions
%
%   **classes:**
%        integers with classes (i.e., networks) to be color-coded on plot
%
%   **overlay:**
%        anatomical underlay image
%
%   **myview:**
%        if entered, shows centers on plot
%
%   **domontage:**
%        if non-zero, make montages of networks
%
%   **[orthview axis]:**
%        Optional; integer for which orthviews axis to use, 1:kk
%
% :Output:
%
%   Cell of colors, in order, for use in other functions
%
% :Examples:
% ::
%
%    classes = c.ClusterSolution.classes;
%    overlay = EXPT.overlay;
%    cluster_orthviews_classes(cl,c.ClusterSolution.classes, EXPT.overlay, 'saggital', 0);
%    cluster_orthviews_classes(cl,c.ClusterSolution.classes, EXPT.overlay, 'axial', 0);
%    cluster_orthviews_classes(cl,c.ClusterSolution.classes, EXPT.overlay, 'coronal', 0);
%    cluster_orthviews_classes(cl,c.ClusterSolution.classes, EXPT.overlay, [], 0);
%
% ..
%    tor wager
%    Colors update, minor improvements, Jan 2010
% ..

whichorth = 1;
if ~isempty(varargin), whichorth = varargin{1}; end

if length(varargin) > 1
    colors2 = varargin{2};
    disp('Using input colors: '); %colors2{:}
else
    % text colors no longer needed -- can handle string colors
    %colors = {'yo' 'bo' 'go' 'ro' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};
    disp('Using default colors. ');
    colors2 = {[1 1 0] [0 0 1] [0 1 0] [1 0 0] [1 .5 0] [0 1 1] [1 0 1] [.5 1 0] [0 .5 1] [0 1 .5] [1 0 .5]};
end

while length(colors2) < max(classes)
    colors2tmp = colors2;
    for jj = 1:length(colors2tmp)
        colors2tmp{jj}(colors2tmp{jj} > .8) = colors2tmp{jj}(colors2tmp{jj} > .8) - .2;
        colors2tmp{jj}(colors2tmp{jj} < .2) = colors2tmp{jj}(colors2tmp{jj} < .2) + .2;
    end
    colors2 = [colors2 colors2tmp];
end
colors2 = colors2(1:max(classes));

% orthviews
cluster_orthviews(cl(classes==1),colors2(1),'overlay',overlay,'solid');

for i = 2:max(classes)
    cluster_orthviews(cl(classes==i),colors2(i),'add','solid', 'handle', whichorth);
end

if ~isempty(myview) && ischar(myview)
    cluster_orthviews_showcenters(cl,myview,overlay,0);
end

if domontage

% montage: build function call
str = ['montage_clusters(overlay,'];
for i = 1:max(classes)
    clc{i} = cl(classes==i);
    str = [str 'clc{' num2str(i) '},'];
end
str = [str '{'];
for i = 1:max(classes)
    %str = [str '''' colors2{i} ''''];  % no longer need solid colors
    str = [str '[' num2str(colors2{i}) '] '];
    %if i ~= max(classes), str = [str ' ']; end
end
str = [str '}, ''nooverlap'');'];
eval(str)
end

return

