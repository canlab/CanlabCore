function cl = region2struct(cl)
% Convert a region object to a simple structure, primarily for
% compatibility with other, older CANlab tools.
%
% :See also: cluster2region, for the reverse transformation

warning off
for i = 1:length(cl)
    cl2(i) = struct(cl(i));
end
warning on

cl = cl2;

end

