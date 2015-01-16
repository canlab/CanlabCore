function outcell = parse_char_to_cell(invar,sepval)
% function outcell = parse_char_to_cell(invar,sepval)
% take a row of characters and separate into cells, breaking at either
% spaces or tabs
%
% tor wager
%
% Example:
% copy from Excel as row, then parse:
% disorder = ['SAD	PTSD	PTSD	PTSD	PTSD	PTSD	SP	SAD	PTSD	PTSD	PTSD	PTSD	PTSD	SAD	SAD	PTSD	PTSD	SP	SP	SP	PTSD	PTSD	PTSD	PTSD	SAD	SAD	SAD	SP	SAD	SAD	SAD	SP	SP	SP	SAD	SP	PTSD	SP	PTSD	PTSD'];
% disorder = parse_char_to_cell(disorder, 'tab');

switch sepval
    case 'space'
whsep = find(invar == ' ');
    case 'tab'
whsep = find(invar == sprintf('\t'));
    otherwise 
        error('Enter space or tab')
end

st = [1 whsep+1];
en = [st(2:end) - 1 length(invar)+1];
n = length(st);

for i = 1:n
    outcell{i} = invar(st(i):en(i)-1);
end

end

    