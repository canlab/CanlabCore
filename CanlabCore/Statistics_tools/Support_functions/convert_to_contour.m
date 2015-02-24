function contourmap=convert_to_contour(groupspace,clus,varargin);

scalefact=100;
if length(varargin)>0;
    varargin{1}=scalefact;
end

fgs=round(groupspace*scalefact);
contourmap=zeros(scalefact,scalefact);

for c=1:max(clus)    
    eachgroupspace=groupspace(find(clus==c),:);
        for n=1:size(eachgroupspace,1);
            contourmap(fgs(n,1),fgs(n,2))=1;
        end
     clear eachgroupspace;
end




