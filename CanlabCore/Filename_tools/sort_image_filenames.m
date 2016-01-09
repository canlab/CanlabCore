function [P, resort] = sort_image_filenames(P)
% :Usage:
% ::
%
%     [P, indices] = sort_image_filenames(P)
%
% Not all image name listing functions return imgs in the correct numbered
% order!
%
% This function resorts a string matrix of image file names by image number, 
% in ascending order
% At most, filename can have one number in it, or error is returned
%
% P can be a string matrix of a cell of string matrices
%
% ..
%    Tor Wager, April 2005
% ..

if iscell(P) % if cell, call this function recursively
    for i = 1:length(P)
        fprintf(1,'%3.0f',i);
        P{i} = sort_image_filenames(P{i});
    end
    fprintf(1,'\n')
    return
end

if size(P,1) == 1, return, end

nums = [];
for j = 1:size(P,1)
    [dd,ff,ee]=fileparts(P(j,:));

    % find the last real number in the name
    tmp = nums_from_text(ff);
    tmp(tmp ~= real(tmp)) = [];
    nums(j) = tmp(end);
end

[n2,resort] = sort(nums);

if any(n2-nums), 
    fprintf(1,'Resorting.')
    
    P = P(resort,:);
end



return
