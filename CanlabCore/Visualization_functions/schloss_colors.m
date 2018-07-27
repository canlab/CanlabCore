function colors = schloss_colors(varargin)
% colors = schloss_colors(varargin)
%
% Color map generated from Karen Schloss's website:
% http://vrl.cs.brown.edu/color
%
% Gramazio, Connor C., David H. Laidlaw, and Karen B. Schloss. 2017. ?Colorgorical: Creating Discriminable and Preferable Color Palettes for Information Visualization.? IEEE Transactions on Visualization and Computer Graphics 23 (1): 521?30.
%
% Goals: Perceptual discrimination, color preference
%
% n = 8
% colors = schloss_colors(n);

colors = { 0.6275 0.9098 0.3569
0.1098 0.3725 0.1176
0.5529 0.8235 0.8471
0.1725 0.4431 0.5804
0.1098 0.8784 0.6980
0.1216 0.2353 0.6510
0.7961 0.4000 0.9176
0.5294 0.0824 0.3137
0.7765 0.7529 0.9961
0.4980 0.1255 0.6745
0.9255 0.4706 0.6588
0.4392 0.4784 0.8941
0.3686 0.2314 0.3412
0.9647 0.1412 0.5608
0.9961 0.0863 0.9569
0.1451 0.1412 0.9765
0.6196 0.4706 0.5843
0.8902 0.8431 0.4118
0.4784 0.1882 0.0118
0.9490 0.6627 0.4000
0.8157 0.1216 0.0941
0.2275 0.6510 0.0353
0.4941 0.6078 0.2392
0.2980 0.9529 0.1725
0.3059 0.2824 0.0353
0.6314 0.8706 0.9412
0.0824 0.2863 0.4588
0.5373 0.9216 0.4824
0.1059 0.3176 0.1137
0.7882 0.8667 0.5294
0.2667 0.2627 0.7059
0.6314 0.5608 0.9725
0.1451 0.5686 0.6157};

for i = 1:size(colors, 1)
    
    catcolors(1, i) = {cat(2, colors{i, :})};
    
end

colors = catcolors;

if ~isempty(varargin)
    
    n = varargin{1};
    
    % add if needed
    while n > length(colors)
        % add
        
        newcolors = colors;
        
        % could rotate color map here - to do in the future!
        
        colors = [colors newcolors];
        
    end
    
    wh = round(linspace(1, length(colors), n));
    
    colors = colors(wh);
    
end


