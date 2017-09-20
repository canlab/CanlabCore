function layer2colors=riverplot_recolor_layer2(sim_matrix,layer1colors)


% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
%    Created July 2016 by Tor Wager
%    
%    8/21/2017 Stephan Geuter
%    changed normalization to used absolute similarity values (l21). Before, negative
%    values could distort the normalization factor resulting RGB values
%    outside the [0-1] interval. 
% ..



[n2, n1]=size(sim_matrix);
for i=1:n2
    temp_color=[0 0 0];
    for j=1:n1
        temp_color= temp_color + layer1colors{j}*abs(sim_matrix(i,j))/sum(abs(sim_matrix(i,:)));
    end
    
    if any(isnan(temp_color));
       temp_color=[.5 .5 .5]; 
    end
    layer2colors{i}=temp_color;
end