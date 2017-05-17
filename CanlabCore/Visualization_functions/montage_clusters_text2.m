function montage_clusters_text2(cl)
% ..
%    Tor Wager
% ..

rc = ceil(sqrt(length(cl))); % how many slices, which ones


for i = 1:length(cl)
    
    subplot(rc,rc,i);
    montage_clusters_maxslice([],cl(i),{'r'});
    
end


return

