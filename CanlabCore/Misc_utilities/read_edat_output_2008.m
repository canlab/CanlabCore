% DATA = read_edat_output_2008(fname, varargin)
%
% Function that creates a structure DATA containing columns
% of the edat file output (saved in text tab delimited "excel" format)
% 
% For this code to work on a Mac, you must: 1) export .edat2 file as an Excel file, 
% then 2) open this file in Excel on a Mac and save as a .csv, 3) read that
% .csv file
%
% Tor Wager, Oct 2008
% 
% Examples:
% -----------------------------------------
% fname = 'myfile.txt'; 
% DATA = read_edat_output_2008(fname)
%
% Defaults: 
% These are the default formats this function expects:
% tab delimited, 1 header row, then row of column names, then data
%
% You can override some of them by using the following --
% E.g., for zero header rows and comma delimited data:
% DATA = read_edat_output_2008(fname, 'nheaderrows', 0, 'mydelimiter', ',')
%
% You can force the number of columns to be a certain value by doing the
% following:
% DATA = read_edat_output_2008(fname, 'nheaderrows', 1, 'numc', 103);
%
% This could be useful if your last row contains empty cells at the end,
% which will mess up the automatic calculation of number of columns.


function DATA = read_edat_output_2008(fname, varargin)
    
% ------------------------------------------------------------------------
% read the database initially to get all column names, etc.
% ------------------------------------------------------------------------
DATA = [];

% assume you have one header row
nheaderrows = 1;
mydelimiter = '\t';
numc = [];

for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % reserved keywords
                case 'nheaderrows', nheaderrows = varargin{i+1};
                case 'mydelimiter', mydelimiter = varargin{i+1};

                case 'numc', numc = varargin{i+1};
                        
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
end
   
col1 = textread(fname, '%s%*[^\n]', 'delimiter', mydelimiter, 'headerlines', nheaderrows, 'bufsize', 100000);

nrows = length(col1);


%numc = 148;
[d] = textread(fname,'%s','headerlines', nheaderrows, 'delimiter', mydelimiter, 'bufsize', 100000);

for i = 1:length(d), if isempty(d{i}), d{i} = 'NaN'; end, end

if isempty(numc)
    % determine from data
    numc = length(d) ./ nrows;
else
    % force numc to be input value
    disp(['Forcing number of columns to be ' num2str(numc)])
    totalel = nrows * numc;
    if length(d) < totalel
        numtoadd = totalel - length(d);
        for i = 1:numtoadd
            d{end+1} = 'NaN';
        end
    else
        d = d(1:totalel);
    end
end

% % if strcmp(d{1},'This file contained edited data.'),
% %     % add an extra input row!
% %     disp('edited data: extra row.')
% %     addone = 1;
% %     [d] = textread(fname,'%s','headerlines',2,'delimiter','\t');
% % end

if mod(length(d), nrows) == 0
    fprintf('This data matrix appears to be OK, with %3.0f rows (including col. names) x %3.0f cols\n', nrows, numc)
else
    fprintf('This data matrix IS NOT OK\n')
    fprintf('i think it has %3.0f header rows, %3.2f data rows (including col. names) and %3.2f cols\n', nheaderrows, nrows, numc)
    
    disp('Here are the first 10 rows of the first column: ')
    disp(char(col1{1:10}))
    
    disp('One potential cause is if you have empty spaces in the last row in your file.')
    disp('then matlab will think it''s reached the end of the file and you will have too few rows.')
    disp('Try entering ''numc'' as an input with the number of columns, or delete the last row.')
    
    return
end

names = d(1:numc)'; disp(['Last name column read is ' names{end}])

disp(' ')
disp('Your column names:')
fprintf('%s ', names{:});
disp(' ')

% ------------------------------------------------------------------------
% replace illegal characters that will cause Matlab problems
% ------------------------------------------------------------------------

for i = 1:length(names), names{i}(findstr('.',names{i})) = '_'; end
for i = 1:length(names), names{i}(findstr('[',names{i})) = ''; end
for i = 1:length(names), names{i}(findstr(']',names{i})) = ''; end
for i = 1:length(names), names{i}(findstr(':',names{i})) = ''; end
for i = 1:length(names), names{i}(findstr('=',names{i})) = ''; end


% ------------------------------------------------------------------------
% Get data matrix in appropriate format
% ------------------------------------------------------------------------

% This only works if there are no blanks
%dmat = reshape(d, numc, nrows)';

% The code below works with blanks

% dmat(1, :) = [];  % get rid of names

% adjust for header rows
nrows = nrows - 1;

fmtstring = '%s%*[^\n]';

for i = 1:numc
    %mycol = dmat(:, i);
    
    % read in the nth column, skipping header rows and column names
    mycol = textread(fname, fmtstring, 'delimiter', mydelimiter, 'headerlines', nheaderrows + 1, 'bufsize', 100000);
    
    fmtstring = ['%*s' fmtstring];  % add this so we skip this col next time
    
    % Get the first non-empty entry and determine if it's numeric or text
    isvalidnumber(i) = 0;
    for j = 1:length(mycol)
        empt = isempty(mycol{j});
        
        if ~empt
            isvalidnumber(i) = ~isempty(str2num(mycol{j}));
            break
        end
    end

    
    if isvalidnumber(i)
        % This column has numbers
        for j = 1:nrows
            try
                mynum = str2num(mycol{j});  % sometimes, with 17:40:12 notation, returns more than 1 number
                if isempty(mynum), mynum = NaN; end
                
            catch
                disp('We seem to be mixing up text and numbers in columns.  Are your rows and cols correct?')
                disp('Here''s the current column (1st 15 rows)')
                mycol(1:15)
                
                disp('Here''s the first 10 rows and cols:')
                dmat(1:10, 1:10)
                
                keyboard
            end
            
            DATA.(names{i})(j, 1) = mynum(1);  
        end
    else
        % This column is text
        DATA.(names{i}) = mycol;
    end
    
end


end



% % 
% % 
% % 
% % % ------------------------------------------------------------------------
% % % create formatting string with number of columns and output var names
% % % ------------------------------------------------------------------------
% % 
% % fmt = repmat('%s',1,numc);
% % fmt = [fmt '%*[^\n]'];
% % 
% % outs = ['[' names{1}];
% % for i = 2:numc, outs = [outs ',' names{i}]; end
% % outs = [outs ']'];
% % 
% % % works, but requires an extra step
% % %outs = ['[out{1}'];
% % %for i = 2:numc, outs = [outs ',out{' num2str(i) '}'];, end
% % %outs = [outs ']'];
% % 
% % % ------------------------------------------------------------------------
% % % read the database again with proper formatting and output names
% % % ------------------------------------------------------------------------
% % 
% % if addone
% %     str = [outs ' = textread(''' fname ''',fmt,''headerlines'',3,''delimiter'',''\t'');'];
% % else
% %     str = [outs ' = textread(''' fname ''',fmt,''headerlines'',2,''delimiter'',''\t'');'];
% % end
% % 
% % eval(str)
% % 
% % warning off
% % ww = whos('*RT*'); ww = ww(end).name; eval(['ww = isempty(str2num(' ww '{1}));'])
% % warning on
% % if ww,
% %     str = [outs ' = textread(''' fname ''',fmt,''headerlines'',3,''delimiter'',''\t'');'];
% %     eval(str)
% % end
% % 
% % % not necessary if all filenames are OK
% % %for i = 1:length(outs),
% % %    eval([names{i} ' = out{i};'])
% % %end
% % 
% % 
% % % ------------------------------------------------------------------------
% % % convert blanks to NaN's, to avoid losing placeholders, and
% % % convert columns with numeric information to numeric vectors
% % % ------------------------------------------------------------------------
% % for i = 1:numc
% %     
% %     % The replacing blanks part
% %     eval(['wh = strmatch('' '',str2mat(' names{i} '{:}));'])
% %     eval([names{i} '(wh) = {''NaN''};'])
% %     
% %     % the conversion to numbers part
% %     eval(['a = str2num(str2mat(' names{i} '{:}));'])
% %     if isempty(a) | sum(isnan(a)) == length(a)
% %         % leave it alone; it's text
% %     else
% %         eval([names{i} ' = a;']), 
% %     end
% % end

