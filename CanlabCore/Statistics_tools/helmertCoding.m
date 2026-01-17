% This code was drafted by Chatgpt in response to the following prompt:
%   "I want to fit a model while controlling for a categorical variable 
%   using helmert codes. How do I do this in matlab? dummyvarcoding does 
%   not include this option so I need to create my own codes."
% Code was verified and modified by BP, 12/07/2024
function helmertMatrix = helmertCoding(categories)
    % Convert to unique levels
    uniqueCategories = unique(categories);
    k = length(uniqueCategories);
    
    % Initialize Helmert matrix
    helmertMatrix = zeros(length(categories), k - 1);
    
    % Generate Helmert contrasts
    for i = 1:k-1
        scale = sum(uniqueCategories >= uniqueCategories(i));
        contrast = (categories == uniqueCategories(i)) - ...
                   (categories > uniqueCategories(i));
        contrast(contrast > 0) = contrast(contrast > 0) * (scale - 1);
        helmertMatrix(:, i) = contrast/scale;
    end
end