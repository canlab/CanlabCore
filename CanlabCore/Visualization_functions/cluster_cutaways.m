function O = cluster_cutaways(clusters,textprefix,mycolor,whichcuts,varargin)
% :Usage:
% ::
%
%    O = cluster_cutaways(clusters,textprefix,mycolor,whichcuts,[coords for center of cuts],[revx])
%
% groups clusters to image on brains by the first letter in whichcuts
% therefore, if you enter 'y' for whichcuts, the program will select
% clusters with similar y values (w/i 15 mm) and image them on the same
% brain.
%
% you need files called brain_render_T1.img/hdr on the path
% these should be scalped anatomicals for the brain image.
%
% you would also need 'brain_render_T1.mat' for the
% transparent brain surface, if you changed this script to get
% the transparent brain surface.
%
% if it cuts from the wrong side, try 'revx' as input into tor_3d.m
% or enter anything as 2nd var arg.
%
% :Examples:
% ::
%
%    O = cluster_cutaways(clusters,'myoutname','y','yzx',[0 0 0],[revx])
%
% Uses renderCluster_ui.m
%
% ..
%    by Tor Wager, Jan 2003
% ..

mm = cat(1,clusters.mm_center);
d = tril(dist(mm(:,2)') < 15);  %distances
bestc = [0 0 0];
revx = 0;

done = zeros(1,length(clusters));

if length(varargin) > 0, bestc = varargin{1};, else bestc = [0 0 0];, end
if length(varargin) > 1, revx = varargin{2};, else revx = 0;, end

clind = 1;
while any(~done)
    
    whclust = (d(:,clind));  
    whclust(find(done)) = 0;
    whclust = find(whclust);
    done(whclust) = 1;
    
        if ~isempty(whclust)
            figure('Color','w');
            O = struct('dohead','y','dobrain','n','get_from','workspace', ...
                'clusters',clusters, ...
            'which_cl',whclust,'whichc',whichcuts,'bestCoords',bestc,'clcol',mycolor,'addtext','n', ...
            'head','brain_render_T1','revx',revx);
        
            O = renderCluster_ui(O); set(O.headp(1:2:end),'FaceColor',[.5 .5 .5]);material dull
            h = findobj('Type','Light');delete(h)
            
            if whichcuts(1) == 'y', view(171,6),[az,el] = view; hh = lightangle(az,el); %hh = lightangle(az,el);
            elseif whichcuts(1) == 'z', view(0,70);[az,el] = view; hh = lightangle(az,el);
            elseif whichcuts(1) == 'x', view(90,10);[az,el] = view; hh = lightangle(az,el);
            end
            
            %h = findobj('Type','Light');delete(h),h(1) = camlight(90,0); h(2) = camlight(270,0); h(3) = camlight(270,0);
            
            if clusters(clind).numVox < 20 & clusters(clind).numVox > 1,
                %hold on, plot3(clusters(1).XYZmm(1,:),clusters(1).XYZmm(2,:),clusters(1).XYZmm(3,:),'o','Color',mycolor,'MarkerFaceColor',mycolor,'MarkerSize',10)
            end
            
            
            saveas(gcf,[textprefix '_' whichcuts '_group' num2str(clind)],'tif');
            saveas(gcf,[textprefix '_' whichcuts '_group' num2str(clind)],'fig');
            close
        end
            
        clind = clind+1;
 end
 
 
 return
 
 
