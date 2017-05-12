function [indic, nms, condf] = string2indicator(str,varargin)
%[indic, nms, condf] = string2indicator(str,varargin)
% 
% Takes a cell vector of string labels and returns an indicator matrix
% Optional argument is a cell vector of strings for which values to match.
%
% Examples:
%[indic,nms] = string2indicator(CL{1}(1).valence);
%[indic,nms] = string2indicator(CL{1}(1).valence,{'neg' 'pos'});

if length(varargin) > 0
    nms = varargin{1};
else
    nms = unique(str);
end

if ~iscell(nms), for i = 1:length(nms), nms2{i} = nms(i);,end, nms = nms2;,end
        
nms(strcmp(nms,'NaN')) = [];

indic = zeros(length(str),1);

for k = 1:length(nms)
    indic(:,k) = strcmp(str,nms{k});
end

if size(nms,1) > size(nms,2)
    nms = nms';
end

if nargout > 2
% make condition function
	condf = indic2condf(indic);
end

return
