%% Circular Graph Examples
% Copyright 2014 The MathWorks, Inc.

%% 1. Adjacency matrix of 1s and 0s
% Create an example adjacency matrix made up of ones and zeros.
rng(0);
x = rand(50);
thresh = 0.93;
x(x >  thresh) = 1;
x(x <= thresh) = 0;

%%
% Call CIRCULARGRAPH with only the adjacency matrix as an argument.
circularGraph(x);

%%
% Click on a node to make the connections that emanate from it more visible
% or less visible. Click on the 'Show All' button to make all nodes and
% their connections visible. Click on the 'Hide All' button to make all
% nodes and their connections less visible.

%% 2. Supply custom properties
% Create an example adjacency matrix made up of various values and supply
% custom properties.
rng(0);
x = rand(20);
thresh = 0.93;
x(x >  thresh) = 1;
x(x <= thresh) = 0;
for i = 1:numel(x)
  if x(i) > 0
    x(i) = rand(1,1);
  end
end

%%
% Create custom node labels
myLabel = cell(length(x));
for i = 1:length(x)
  myLabel{i} = num2str(round(1000000*rand(1,1)));
end

%%
% Create custom colormap
figure;
myColorMap = lines(length(x));

circularGraph(x,'Colormap',myColorMap,'Label',myLabel);