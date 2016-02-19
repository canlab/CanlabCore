function ov_mat = cl_overlap(in1,in2)

% :Usage:
% ::
%
%    cl_overlap(in1,in2)
%
% purpose:  compute overlap of two clusters
%
% :Input:
%
%   **in1**
%        
%       cluster object
%
%   **in1**
%        
%       cluster object
%
%
% :Output:
%
%   **ov_mat**
% 
%       matrix of overlaps
%
%

for i = 1:length(in1)
    
    for j = 1:length(in2)
        
        a = in1(i);
        b = in2(j);

        if size(a.XYZ,2) < size(b.XYZ,2)
            a.XYZ(:,end+1:(end+size(b.XYZ,2) - size(a.XYZ,2))) = NaN;
        elseif size(b.XYZ,2) < size(a.XYZ,2)
            b.XYZ(:,end+1:(end+size(a.XYZ,2) - size(b.XYZ,2))) = NaN;
        end
    
        ov_mat(i,j) = sum(sum(a.XYZ == b.XYZ) == 3);

    end
    
end

return