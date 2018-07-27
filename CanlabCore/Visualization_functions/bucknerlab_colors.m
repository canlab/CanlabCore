function [colorcell, colormap] = bucknerlab_colors(varargin)

rgb = [120  18 134
255   0   0
70 130 180
42 204 164
74 155  60
0 118  14
196  58 250
255 152 213
220 248 164
122 135  50
119 140 176
230 148  34
135  50  74
12  48 255
0   0 130
255 255   0
205  62  78];

colormap = rgb ./ 255;

try % diff versions of matlab are different?
    
    colorcell = mat2cell(colormap, ones(length(rgb), 1), 3);
    
catch
    
    colorcell = mat2cell(colormap, ones(length(rgb), 1));
end

end

