function texth = add_text_label(cl, mytext, mycolor, varargin)
%
% texth = add_text_label(cl, mytext, mycolor, varargin)
% 
% For making movies and adding text easily

yplus = 0;
xplus = 0;
zplus = 0;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}

            % functional commands
            case 'yplus', yplus = varargin{i+1};
            case 'xplus', xplus = varargin{i+1};
            case 'zplus', zplus = varargin{i+1};
    
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end
    
x = cl(1).mm_center(1) - 10 + xplus;
y = cl(1).mm_center(2) - 20 + yplus;
z = cl(1).mm_center(3) - 0 + zplus;

texth = text(x, y, z, mytext,'Color',mycolor,'FontWeight','b','FontSize',36);

end
