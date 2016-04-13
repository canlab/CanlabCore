function AUs = CERTreader(fname)
% Given an output file from CERT, returns a struct with the data
%
% :Input:
%
%   **fname:**
%        filename
%
% :Output:
%
%   struct with the data


r=robustcsvread(fname, 'rows_to_skip', 1, 'delim', '\t');

fs = fields(r);
for i=1:length(fs)
    r.(fs{i}) = str2num(char(r.(fs{i})));
end

AUs=r;
