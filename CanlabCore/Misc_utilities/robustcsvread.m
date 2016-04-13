function MM=robustcsvread(filename, varargin)
% ROBUSTCSVREAD reads in CSV files with different number of columns
% on different lines
%
% This returns a struct, with one field per column of the csv file.
% Each field is a cell array whose length = rows in the csv file.  Column
% names are assumed to be in the first row.
% If column names are invalid struct field names, edits them by replacing
% funky characters with an underscore, or if first char is a number, I
% prepend aa_ to the field name.
%
% :Inputs:
%
%   **varargin:**
%        cols:  how many cols to read in, by defaults reads them all
%
%        rows_to_skip:  how many rows to skip
%
%        delim:  cell delimiter
%
%        missing:  followed by cell array, first cell is val for missing,
%                   second cell is what to replace with
%
% ..
%    extended by Yoni Ashar, 10/2012
%
%    original code off the fileexchange, by
%    robbins@bloomberg.net
%    michael.robbins@us.cibc.com
% ..

cols = NaN; 
n_h = 0;
delim = ',';
missing_data = 0;
for i=1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'cols'
            cols = varargin{i+1};
            case 'rows_to_skip'
                n_h = varargin{i+1};
            case 'delim'
                delim = varargin{i+1};
            case 'missing'
                missing_data = 1;
                missing = varargin{i+1}{1};
                replace_missing = varargin{i+1}{2};
        end
    end
end


if ~exist(filename, 'file'), MM=[]; warning('File does not exist'); fprintf('%s\n', filename), return, end

fid=fopen(filename,'r');
slurp=fscanf(fid,'%c');
fclose(fid);
M=strread(slurp,'%s','delimiter','\n');
for i=1:length(M)
    temp=strread(M{i},'%s','delimiter',delim);
    for j=1:length(temp)
        if missing_data && strcmp(temp{j}, missing), temp{j} = replace_missing; end
        MM{i,j}=temp{j};        
    end;
end;


% Yoni's extension is below
data = MM;
clear MM;

for i=1:n_h
    data(1,:) = [];
end

%one line silly hack for a project Yoni needed... please forgive.
if isequal(data{1,2},'"ANSWER - ENTER A'), data{1,2} = 'A'; end

for i=1:min(cols,size(data,2))
    data{1,i} = strtrim(data{1,i});
    % replace bad characters with _
    data{1,i} = regexprep(data{1,i}, '\W', '_');
    % if col name starts with a number, preprend an 'aa'
    if isstrprop(data{1,i}(1), 'digit'), data{1,i} = ['aa_' data{1,i}]; end
    
    cmd = ['MM.' data{1,i} ' = data(2:end,' num2str(i) ');'];
    eval(cmd);
end

