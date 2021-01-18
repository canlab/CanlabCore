%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Funded by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function [commaFormattedString] = CommaFormat(value)
  % Split into integer part and fractional part.
  [integerPart, decimalPart]=strtok(num2str(value),'.'); 
  % Reverse the integer-part string.
  integerPart=integerPart(end:-1:1); 
  % Insert commas every third entry.
  integerPart=[sscanf(integerPart,'%c',[3,inf])' ... 
      repmat(',',ceil(length(integerPart)/3),1)]'; 
  integerPart=integerPart(:)'; 
  % Strip off any trailing commas.
  integerPart=deblank(integerPart(1:(end-1)));
  % Piece the integer part and fractional part back together again.
  commaFormattedString = [integerPart(end:-1:1) decimalPart];
  return; % CommaFormat
end