function conj = conjunction(si1, si2, varargin)
% Returns the conjunction of two statistic_images.  considers positive and
% negative activations separately.
%
% :Inputs:
%
%   Two thresholded statistic images.  Optional 3rd argument:  -1 to
% get only negative conjunction, or 1 to get only positive conjunction
%
% :Output:
%
%   A statistic_image with all voxels suprathreshold (in the same direction) in both input
% images.  Voxel values are set to 1 and -1, to indicate direction.
%
% ..
%    Yoni Ashar, Sept. 2015
% ..

direction = 0; % default: both pos and neg
if nargin > 2
    direction = varargin{1};
end

if compare_space(si1, si2), error('Objects are not in same space'); end

% split first si into pos/neg significant voxels
si1_pos = si1;
posinds = si1.dat > 0;
si1_pos.sig(~posinds) = 0; % unset significance of neg voxels

si1_neg = si1;
neginds = si1.dat < 0;
si1_neg.sig(~neginds) = 0; % unset significance of pos voxels

% repeat for 2nd si
si2_pos = si2;
posinds = si2.dat > 0;
si2_pos.sig(~posinds) = 0;

si2_neg = si2;
neginds = si2.dat < 0;
si2_neg.sig(~neginds) = 0;

% the pos conjuction
pos_conj = si1_pos;
pos_conj.dat( pos_conj.dat > 0) = 1; % set all pos values to 1
pos_conj.sig = si1_pos.sig & si2_pos.sig;

% the neg conjuction
neg_conj = si1_neg;
neg_conj.dat( neg_conj.dat < 0) = -1; % set all neg values to -1
neg_conj.sig = si1_neg.sig & si2_neg.sig;

if direction==0 % both pos and neg
    conj = si1; % to load image space, etc.
    conj.dat = pos_conj.dat + neg_conj.dat;
    conj.sig = pos_conj.sig + neg_conj.sig;

elseif direction < 0 % only neg
    conj = neg_conj;

else % only pos
    conj = pos_conj;
end
        
    
