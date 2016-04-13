function lineH = cl_line_plots(clusters)

% :Usage:
% ::
%
%    cl_line_plots(clusters)
%
% purpose:  ???
%
% :Input:
%
%   **clusters**
%        
%       cluster object
%
% :Output:
%
%   **lineH**
% 
%       handle
%
%



    index = 1;
    tabtext = [];
    
    for i = 1:length(clusters)
        
        cl = clusters(i);
        
        cl_center = mean(cl.XYZmm');
        
        line_out = cl_center + sign(cl_center) .* 125;
        line_out(abs(line_out) > 125) = sign(abs(line_out) > 125) * 125;
        %line_out = cl_center + (sign(cl_center) .* ([100 100 100] ./ abs(cl_center)))
        %line_out(abs(line_out) > 125) = sign(line_out(abs(line_out) > 125)) .* 125;
        %line_out
        
        text_pos = line_out + sign(cl_center) * 2;
        
        xyz = [cl_center' line_out'];
        x = xyz(1,:);
        y = xyz(2,:);
        z = xyz(3,:);
        
        hold on
        
        %lineH(index) = plot3(x,y,z,'k','LineWidth',1.5);
        index = index + 1;
        
        linetext = [num2str(i) ' (' num2str(cl.numVox) ' voxels): max t = ' num2str(max(clusters(i).Z)) ' mni: ' num2str(round(cl_center))];
        tabtext = str2mat(tabtext,linetext);
        
        %lineH(index) = text(text_pos(1),text_pos(2),text_pos(3),num2str(i));
        lineH(index) = text(cl_center(1),cl_center(2),cl_center(3),num2str(i),'FontSize',12);
        index = index + 1;
    end
    
    tabtext
return
    