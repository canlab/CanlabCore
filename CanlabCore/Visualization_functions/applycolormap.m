function h = applycolormap(h,mapname)
% ..
%    by Tor Wager
% ..

mp2 = colormap; % save current colormap
mp = colormap(mapname);

if isfield(get(h),'FaceVertexCData')    % for a patch object with faces/vertices
    dat = get(h,'FaceVertexCData');
    
    if size(dat,2) > 1, error('FaceVertexCData is already RGB color?');, end
    
    % scale dat to indices
    minc = 1 - min(dat);
    maxc = max(dat) ./ (size(mp,1) - 1);
    dat = dat ./ maxc;
    dat = dat + minc;
    
    d2 = round(dat);
    d2(d2 == 0) = 1;
    d2(d2 > size(mp,1)) = size(mp,1);


    for i = 1:length(dat)
        rgb(i,:) = mp(d2(i),:);
    end

    set(h,'FaceVertexCData',rgb)



    
else    % for 2-d image color data
    dat = get(h,'CData');
    
    if length(size(dat)) > 2, error('CData is already RGB color?');, end
    
    % scale dat to indices
    minc = 1 - min(min(dat));
    maxc = max(max(dat)) ./ (size(mp,1) - 1);
    dat = dat ./ maxc;
    dat = dat + minc;
    
    d2 = round(dat);
    d2(d2 == 0) = 1;
    d2(d2 > size(mp,1)) = size(mp,1);

    for i = 1:size(dat,1)
        for j = 1:size(dat,2)
            val = d2(i,j);
            rgb(i,j,:) = mp(val,:);
        end
    end

    set(h,'CData',rgb)

end


colormap(mp2)   % reapply 

return
