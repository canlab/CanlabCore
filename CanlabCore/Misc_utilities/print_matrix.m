function print_matrix(x,varargin)
% prints matrix values as tab delimited, 2 decimal places
%
% :Usage:
% ::
%
%     print_matrix(x,[col names cell array], [row names cell], [format string])
%
% :Examples:
% ::
%
%    t = [1 2; 3 4; 5 6];
%    print_matrix(t,{'col1' 'col2'},{'row1' 'row2' 'row3');
%
%    print_matrix(rand(5), [], [], '%3.2f');
%    print_matrix(rand(5), {'A' 'B' 'C' 'D' 'E'}, {'A' 'B' 'C' 'D' 'E'}, '%3.2f');
%    print_matrix(rand(5), {'A' 'B' 'C' 'D' 'E'}, {'A' 'B' 'C' 'D' 'E'}, '%d');
%
% ..
%    tor wager
% ..

fmtstring = '%3.4f';
if length(varargin) > 2 && ~isempty(varargin{3})
    fmtstring = varargin{3};
end

% set up
s = size(x,2);
str = [repmat([fmtstring '\t'],1,s) '\n'];
colnames = ''; % default added by Luka, 5/2013

if length(varargin) > 0 && ~isempty(varargin{1})
    colnames = varargin{1};
    
    if ~isempty(colnames) && iscolumn(colnames), colnames = colnames'; end

end

if length(varargin) > 1 && ~isempty(varargin{2})

    rownames = varargin{2};
    if ~isempty(rownames) && iscolumn(rownames), rownames = rownames'; end

    colnames = [{' '} colnames];
end 


% Print names
if length(varargin) > 0 && ~isempty(varargin{1})
    for i = 1:length(colnames)
        fprintf(1,'%s\t',colnames{i}), 
    end
    fprintf(1,'\n')
end

% print matrix
if length(varargin) > 1 && ~isempty(varargin{2})
    for i = 1:size(x,1)
        fprintf(1,'%s\t',rownames{i})
        fprintf(1,str,x(i,:));
    end
else
    disp(sprintf(str,x'))
end




    


end
