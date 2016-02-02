function M1 = glassbrain_avi(fps,len,clusters)
% :Usage:
% ::
%
%    M = glassbrain_avi(fps,len,clusters)
%
% works with cluster_surf.m
% assumes a gray brain where no activation, or creates one
%
% SEE ALSO cluster_surf_movie.m, which works well too.
%
% ..
%    tor wager 
% ..

dofadein = 0;
dozionimg = 0;
if isempty(clusters), doclustersurf = 0;, else, doclustersurf = 1;,end


if dozionimg
% ---------------------------------------
% Special Script
% ---------------------------------------
hh = open('Brook Zion 2.jpg');
hh = hh.Brook_Zion_2;
for i = 1:3, hh2(:,:,i) = flipud(hh(:,:,i));,end
hh = hh2;
f1 = figure('Color',[0 0 .7]); image(-size(hh,2)./2,-size(hh,1)./2,hh); hold on; axis image; axis auto
set(gca,'Color',[0 0 1],'YDir','Normal'); axis off
Ps = which('surf_single_subj_T1_gray.mat'); %'c:\tor_scripts\3DheadUtility\surf_single_subj_T1_gray.mat';
load(Ps)
 h = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5], ...
  'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',0.01,'SpecularExponent',200);

%camzoom(3)
camzoom(1.2)
drawnow

 O = struct('name','glassZion.avi','fps',5,'length',10,'timeres',1, 'H',h, ...
    'azOffset',0,'elOffset',0,'zoom',2,'add2movie',[],'closemovie',1);

    finalTrans = .5;
    O.timecourse{1} = [0:(finalTrans/(O.fps*O.length)):finalTrans] + .01;
    
    M2 = make3Davi(O);  % close this one
    O.name = 'glassZion_act.avi'; O.closemovie = 0;
    camzoom(1/3); set(h,'FaceAlpha',.01);
    M1 = make3Davi(O);  % leave this one open
    
    axis vis3d
    

% ---------------------------------------
% First part: transparent brain fade-in
% ---------------------------------------
elseif dofadein
    f1 = figure('Color','w'); h = addbrain;
    view(0,90); axis off;

    set(h,'Tag','cluster1');

    O = struct('name','glassbrainX.avi','fps',fps,'length',len(1),'timeres',1, 'H',h, ...
    'azOffset',0,'elOffset',0,'zoom',1,'add2movie',[],'closemovie',0);

    finalTrans = .5;
    O.timecourse{1} = [0:(finalTrans/(O.fps*O.length)):finalTrans];
    M1 = make3Davi(O);

else
   
    f1 = gcf;
    h3 = get(gca,'Children');for i = 1:length(h3),h4=get(h3(i));if isfield(h4,'FaceVertexCData'),h5 = h3(i);,end,end
    h = h5; %gco;
    tmp = get(h,'FaceVertexCData'); if size(tmp,2) ~= 3 | ~(size(tmp,1) > 5000), error('Make surface object current object'),end
    
    O = struct('name','glassbrainX.avi','fps',fps,'length',len(1),'timeres',1, 'H',h, ...
    'azOffset',0,'elOffset',0,'zoom',1,'add2movie',[],'closemovie',0);
    M1 = avifile('glassbrainX.avi','Quality',75,'Compression','none','Fps',O.fps);

    
end



% ---------------------------------------
% Second part: change cluster surface
% ---------------------------------------
if doclustersurf
    cl = cat(2,clusters.Z); refZ = [0 max(cl) min(cl(cl < 0)) 0];, 

    %for j = 1:length(clusters)
    %      clusters(j).origZ = clusters(j).Z;
    %  end
    %cluster_surf(clusters,h,refZ,10,{[1 0 0]},'colorscale','heatmap');
    
    cluster_surf(clusters,h,10,{[1 0 0]},'colorscale','heatmap');
    figure(f1)
end
    
% ---------------------------------------
% Second part: cluster activations
% ---------------------------------------
axis vis3d    
mov = avifile('glassbrainX_act1.avi','Quality',75,'Compression','none','Fps',O.fps);
    
% add frames

tmp = get(h,'FaceVertexCData'); 
t1 = find(~(all(tmp == .5,2))); 
ind = 1;
    
for i = [10:-1:1 1:10 10:-1:1 1:10]

    if ind > 10, ff = get(h,'FaceAlpha'); set(h,'FaceAlpha',min(1,ff+.1));, end 
    fprintf(1,'%3.0f . ',i)
    M1 = get_frame(M1,i,tmp,h,t1);
    mov = addframe(mov,gca);
    ind = ind + 1;
    
end

mov = close(mov);


% ---------------------------------------
% Third part: tip down
% ---------------------------------------
[az,el] = view; ind = 1;
ii = repmat([10:-1:1 1:10 10:-1:1 1:10 10:-1:1 1:10],1,10);

mov = avifile('glassbrainX_act2.avi','Quality',75,'Compression','none','Fps',O.fps);

for i = 5:5:90

    view(az-2*i,el - i)
    
    fprintf(1,'%3.0f . ',i)
    M1 = get_frame(M1,ii(ind),tmp,h,t1);
    mov = addframe(mov,gca);
    
    ind = ind + 1;
    lightFollowView
end

mov = close(mov);

% ---------------------------------------
% Fourth part: spin
% ---------------------------------------
[az,el] = view; ind = 1;
ii = repmat([10:-1:1 1:10 10:-1:1 1:10 10:-1:1 1:10],1,10);

mov = avifile('glassbrainX_act3.avi','Quality',75,'Compression','none','Fps',O.fps);

for i = 0:10:360

    if i < 70
        view(az-i,el+i./2)
    elseif i == 70
        [dummy,el] = view;
        view(az-i,el)
    else
        view(az-i,el);
    end
    
    fprintf(1,'%3.0f . ',i)
    M1 = get_frame(M1,ii(ind),tmp,h,t1);
    mov = addframe(mov,gca);
    
    ind = ind + 1;
    lightFollowView
end

mov = close(mov);

M1 = close(M1);

    
return



% ---------------------------------------
% Alternative: Random blobs
% ---------------------------------------
n = noisevector(size(tmp,1),[1:-.1:0],1);
gray = .5 * ones(size(tmp));
set(h,'FaceVertexCData',gray);
t2 = gray;

ind = 1;
ni = find(n > 1);
 bins = round(1:length(ni)./30:length(ni));
 cols = rand(length(bins),3);
cols(:,3) = 0;

for i = 1:length(bins)-1
    whh{ind} = ni(ni > bins(i) & ni < bins(i+1)) ;  
    t2(whh{ind},:) = repmat(cols(ind,:),length(whh{ind}),1);
    ind = ind+1;
end
set(h,'FaceVertexCData',t2);




function M1 = get_frame(M1,i,tmp,h,t1)

    fprintf(1,'%3.0f . ',i)
    t2 = tmp; t2(t1,:) = 2.5 .* t2(t1,:) ./ (i .^ .7); % scale down

    t2(t1,:) = t2(t1,:) + .05 * randn(length(t2(t1,1)),3); % add a little noise
        t2(t2 > 1) = 1;
        t2(t1,3) = 0;        % no blue
        
    t2(sum(t2,2) < .2,:) = .5;                  % replace low values with gray
    set(h,'FaceVertexCData',t2);
    drawnow
    
    H = gca;
    %mov = addframe(mov,H);
    M1 = addframe(M1,H);
    
 return
    
    
