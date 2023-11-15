function print_matrix(x,varargin)
% prints matrix values as tab delimited, 2 decimal places
%
% :Usage:
% ::
%
%     print_matrix(x,[col names cell array], [left row names cell], [format string], [right row names cell])
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
% 2023/10/23 Michael Sun
% - Now includes the option to include a 'right_rowlabel', for example, to
% include a column for significance stars.


fmtstring = '%3.4f';
if length(varargin) > 2 && ~isempty(varargin{3})
    fmtstring = varargin{3};
end

% set up
s = size(x,2);
str = [repmat([fmtstring '\t'],1,s)];
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

if length(varargin) > 3 && ~isempty(varargin{4})

    right_rownames = varargin{4};
    if ~isempty(right_rownames) && iscolumn(right_rownames), right_rownames = right_rownames'; end

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
        if exist('right_rownames', 'var')
            fprintf(1, '%s', right_rownames{i})
        end
        fprintf(1,'\n')
    end

else
    disp(sprintf(str,x'))
end




    


end
