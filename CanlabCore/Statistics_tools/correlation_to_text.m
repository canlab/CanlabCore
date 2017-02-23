function [str,sigmat] = correlation_to_text(a,cm,varargin)
%[str,sigmat] = correlation_to_text(matrix,value to mark with *,[cell array%of names])
%
% input correlation matrix and critical value
% or raw data and empty cm value
% (determines threshold based on n)
%
% see also print_correlation.m

if length(varargin) > 0, nms = varargin{1};,else, nms = {[]};,end

if isempty(cm)
    % choose default value
    [rci,sig,z,p,cm] = r2z(.5,size(a,1),.05);
    a = corrcoef(a);
end
    
% names

str=sprintf('Crit=%3.2f\t',cm);
for j=1:size(a,2),
    
    if length(nms) < j, nms{j} = ['V' num2str(j)];,end
    
    str=[str sprintf('%s\t',nms{j})];,
end
str=[str sprintf('\n')];,
    

% table
sigmat=a*0;
for i = 1:size(a,1),
    
    str=[str sprintf('%s\t',nms{i})];,
    
    for j=1:size(a,2),
        
        if abs(a(i,j))>cm & a(i,j) ~= 1,t='*';,
            sigmat(i,j)=1;
        else,t='';,
            
        end,
        
        str=[str sprintf('%3.3f%s\t',a(i,j),t)];,
    end,
    str=[str sprintf('\n')];,
end

disp(str)

return
