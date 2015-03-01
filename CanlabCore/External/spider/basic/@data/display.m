function [s,it]=display(d,tab,it)


if nargin>1
    [s,it]=display(d.algorithm,tab,it);
    return;
end

disp(d.algorithm.name);
s=(['X = ' num2str(size(d.X,1)) 'x' num2str(size(d.X,2))  ]);
% if size(d.X,1)<=10 & size(d.X,2)<=10  & size(d.Y,1)<=10 & size(d.Y,2)<=10       
%     if size(d.X,2)==1 & size(d.Y,2)==1 
%         if ~isstruct(d.X(1,1))
%             disp('     X     Y'); disp([d.X d.Y]);
%         end
%     else
%         if ~isempty(d.X) & prod(size(d.X))>1  disp('X:'); disp(d.X); end
%         if ~isempty(d.Y) & prod(size(d.X))>1  disp('Y:'); disp(d.Y); end
%     end
%  else
     disp(['data dimensions:'])
    disp([s '  Y = ' num2str(size(d.Y,1)) 'x' num2str(size(d.Y,2)) ]);
% end 

