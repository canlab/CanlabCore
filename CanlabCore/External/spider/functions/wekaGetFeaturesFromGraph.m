function Feats=wekaGetFeaturesFromGraph(a)
% this function retrieves the used features from a weka tree "Drawable"
% object
if(~isempty(strfind(class(a),'weka.classifiers.trees')))
    s=char(a.graph);
else
    s=char(a);    
end
Feats=[];
I=strfind(s,'inp');
for i=I
    f=sscanf(s(i+3:end),'%d');
    Feats=[Feats,f];
end

Feats=unique(Feats);