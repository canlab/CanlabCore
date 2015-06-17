function out = sampleInSeconds(stimList,ISI,varargin)
% out = sampleInSeconds(stimList,ISI,varargin)
% input: stimList, output: stimlist sampled in .1 seconds
%                          OR sampled at your specified frequency

scale = ceil(ISI*10);

if nargin > 2
    scale = round(ISI/varargin{1});
end

numstim = size(stimList,1);

out = zeros(numstim*scale,1);

for i = 0:numstim-1

    out(i*scale+1,1) = stimList(i+1,1);

end

return

