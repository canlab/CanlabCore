function com = center_of_mass(XYZ,Z)
% This function returns the center of mass of a cluster of voxels or mm coordinates
% defined as the nearest in-list coordinate to the average of the 
% coordinate values weighted by the Z-score
%
% :Usage:
% ::
%
%     com = center_of_mass(XYZ,Z)
%
% assigns a rank to each coordinate based on Z scores
% and includes 
%
% enter a 3 x n list of XYZ coordinates
% returns 1 x 3 center of mass
%
% ..
%   by Tor Wager
%
%   Functions called:
%   C:\matlabR12\toolbox\matlab\elmat\repmat.m
%
%   Edits: Oct 2012, to fix bug if Z values sum to exactly zero
%   Now scales weights by abs value, instead of original value (Z)
% ..

if size(Z,1) > size(Z,2), Z = Z'; end

if all(Z) < 10*eps
    disp('Center_of_mass warning: All Z values are zero');
    Z = ones(size(Z));
end

% Z = Z ./ sum(Z);

Z = repmat(Z,3,1);

if any(isnan(Z(:)) | isinf(Z(:))), Z = 1; end

XYZw = XYZ .* abs(Z); % Wani modified Z into abs(Z), 11/12/12

%center = sum(XYZw,2) ./ sum(Z,2);
center = sum(XYZw,2) ./ sum(abs(Z),2);

center = repmat(center,1,size(XYZ,2));

whch = find(sum((XYZ - center).^2) == min(sum((XYZ - center).^2)));

com = XYZ(:,whch)';
com = com(1,:);

return
