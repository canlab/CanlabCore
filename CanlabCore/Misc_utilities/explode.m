function [split,numpieces]=explode(string,delimiters)
%EXPLODE    Splits string into pieces.
%   EXPLODE(STRING,DELIMITERS) returns a cell array with the pieces
%   of STRING found between any of the characters in DELIMITERS.
%
%   [SPLIT,NUMPIECES] = EXPLODE(STRING,DELIMITERS) also returns the
%   number of pieces found in STRING.
%
%   Input arguments:
%      STRING - the string to split (string)
%      DELIMITERS - the delimiter characters (string)
%   Output arguments:
%      SPLIT - the split string (cell array), each cell is a piece
%      NUMPIECES - the number of pieces found (integer)
%
%   Example:
%      STRING = 'ab_c,d,e fgh'
%      DELIMITERS = '_,'
%      [SPLIT,NUMPIECES] = EXPLODE(STRING,DELIMITERS)
%      SPLIT = 'ab'    'c'    'd'    'e fgh'
%      NUMPIECES = 4
%
%   See also IMPLODE, STRTOK
%
%   Created: Sara Silva (sara@itqb.unl.pt) - 2002.04.30
%
%   Modified: Matthew Davidson (matthew@psych.columbia.edu) - 2005.10.11

if isempty(string) % empty string, return empty and 0 pieces
    split{1}='';
    numpieces=0;

elseif isempty(delimiters) % no delimiters, return whole string in 1 piece
    split{1}=string;
    numpieces=1;

else % non-empty string and delimiters, the correct case

    remainder=string;
    i=0;
    
    % If the first letter is in the delimiters, add a blank piece at the
    % beginning
    firstletteridx = strfind(delimiters,string(1));
    if ~isempty(firstletteridx)
        i=i+1;
        split{i}='';
        remainder = string(2:end);
    end

    while ~isempty(remainder)
        [piece,remainder]=strtok(remainder,delimiters);
        i=i+1;
        split{i}=piece;
    end
    numpieces=i;

end
