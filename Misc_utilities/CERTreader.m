% given an output file from CERT, returns a struct with the data
%
% input:
%   filename
%
% output:
%   struct with the data

function AUs = CERTreader(fname)

r=robustcsvread(fname, 'rows_to_skip', 1, 'delim', '\t');

fs = fields(r);
for i=1:length(fs)
    r.(fs{i}) = str2num(char(r.(fs{i})));
end

AUs=r;