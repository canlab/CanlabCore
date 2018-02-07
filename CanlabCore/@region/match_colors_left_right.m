function [all_colors, leftmatched, rightmatched, midline, leftunmatched, rightunmatched] = match_colors_left_right(r, varargin)
% Given a region object and optional colorfun, generate a list of colors
% for each region, assigning the same color to symmetric regions in left
% and right hemispheres.
%
% [all_colors, leftmatched, rightmatched, midline, leftunmatched, rightunmatched] = match_colors_left_right(r, varargin)
%
% - input custom color function
% [all_colors, leftmatched, rightmatched, midline, leftunmatched, rightunmatched] = match_colors_left_right(r, @(n) custom_colors([1 0 .5], [.5 0 1], n));

colorfun = @scn_standard_colors;

if length(varargin) > 0, colorfun = varargin{1}; end

% Get left, right, midline
% -------------------------------------------------------------------------

xyz = cat(1, r.mm_center);

lr = sign(xyz(:, 1));

wh = lr < 1 & lr > -1;
midline = r(wh);

isleft = lr <= -1;
left = r(isleft);

isright = lr >= 1;
right = r(isright);

whleft = find(isleft);  % index into full list
whright = find(isright);  % index into full list

% -------------------------------------------------------------------------

if ~any(whleft) | ~any(whright)
    % cannot match regions
    
    all_colors = colorfun(length(r));
    return
    
end

% Match regions
% -------------------------------------------------------------------------

for i = 1:length(whleft)
    
    rr = r(whleft(i));
    
    %d = dist(rr.mm_center(2:3), xyz(isright, 2:3)')';  % index into rightregions; neural network toolbox now!
    d = distance_euclid(rr.mm_center(2:3), xyz(isright, 2:3));  % index into right regions
    
    [mind, whmin] = min(d); % closest one on right
    rmatch = right(whmin);
    
    % rr must be closest one to right-hand match too
    %d = dist(rmatch.mm_center(2:3), xyz(isleft, 2:3)')';  % index into right regions
    d = distance_euclid(rmatch.mm_center(2:3), xyz(isleft, 2:3));  % index into right regions

    [mind2, whmin2] = min(d); % closest one on left

    if whmin2 == i
        whrightmatch(i, 1) = whmin;  % index into right regions of closest matching on right        
    else
        whrightmatch(i, 1) = 0;
    end
end

% Clusters
leftmatched = left(whrightmatch > 0);
rightmatched = right(whrightmatch(whrightmatch > 0));

leftunmatched = left(~whrightmatch);
rightunmatched = right;
rightunmatched(whrightmatch(whrightmatch > 0)) = [];

%%

% index into original vectors
index_leftmatch = whleft(whrightmatch > 0);

index_rightmatch = whright(whrightmatch(whrightmatch > 0));

index_midline = find(lr < 1 & lr > -1);

index_leftunmatched = whleft(~whrightmatch);

z = ones(size(right));
z(whrightmatch(whrightmatch > 0)) = 0;
index_rightunmatched = whright(find(z));

nmatch = length(index_leftmatch);
nextra = length(index_midline) + length(index_leftunmatched) + length(index_rightunmatched);

c = colorfun(nextra + nmatch);

all_colors(index_leftmatch) = c(1:nmatch);
all_colors(index_rightmatch) = c(1:nmatch);
c(1:nmatch) = [];

n = length(index_midline);
all_colors(index_midline) = c(1:n);
c(1:n) = [];

n = length(index_leftunmatched);
all_colors(index_leftunmatched) = c(1:n);
c(1:n) = [];

n = length(index_rightunmatched);
all_colors(index_rightunmatched) = c(1:n);
c(1:n) = [];

end % function


%%
% whright = find(isright);  % index into full list
% 
% for i = 1:length(whright)
%     
%     rr = r(whright(i));
%     
%     d = dist(rr.mm_center(2:3), xyz(isleft, 2:3)')';  % index into left regions
%     
%     [mind, whmin] = min(d); % closest one on left
%     lmatch = left(whmin);
%     
%     % rr must be closest one to right-hand match too
%     d = dist(lmatch.mm_center(2:3), xyz(isright, 2:3)')';  % index into left regions
% 
%     [mind2, whmin2] = min(d); % closest one on right
%     rmatch = right(whmin2);
%     
%     if whmin2 == i
%         whleftmatch(i, 1) = whmin; 
%     else
%         whleftmatch(i, 1) = 0;
%     end
% end
% 
% rightmatched2 = right(whleftmatch > 0);
% leftmatched2 = left(whleftmatch(whleftmatch > 0));
% 
% %%
% unique([whleftmatch(whleftmatch > 0); find(whrightmatch > 0)])
% 
% unique([whrightmatch(whrightmatch > 0); find(whleftmatch > 0)])

%cluster_orthviews(midline, 'unique');

