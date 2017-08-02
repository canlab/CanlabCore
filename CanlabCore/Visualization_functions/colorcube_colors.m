function colors = colorcube_colors(varargin)

wh = [10 11 13 17 18 21 22 23 28 30 32 24 36];

nofigs = isempty(findobj(0, 'Type', 'Figure'));

cm = colormap('colorcube');
cm = [cm(wh, :); cm];

% Maybe different versions of mat2cell? Not so stable...
nc = ones(size(cm, 1), 1);
try
    cm = mat2cell(cm, nc, 3)';
catch
    cm = mat2cell(cm, ones(nc, 1))';
end

colors = cm;

if nofigs
    % we created one for colormap, so close
    close
end

if ~isempty(varargin)
    n = varargin{1};
    
    while length(colors) < n, colors = [colors colors]; end
    
    colors = colors(1:n);
end

end