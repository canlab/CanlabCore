function layer2colors=riverplot_recolor_layer2(sim_matrix,layer1colors)

[n2, n1]=size(sim_matrix);
for i=1:n2
    temp_color=[0 0 0];
    for j=1:n1
    temp_color= temp_color + layer1colors{j}*sim_matrix(i,j)/sum(sim_matrix(i,:));
    end
    layer2colors{i}=temp_color;
end