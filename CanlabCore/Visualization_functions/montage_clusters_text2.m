function montage_clusters_text2(cl)
% montage_clusters_text2(cl)
% tor wager


% how many slices, which ones
rc = ceil(sqrt(length(cl)));


for i = 1:length(cl)
    
    subplot(rc,rc,i);
    montage_clusters_maxslice([],cl(i),{'r'});
    
end


return

