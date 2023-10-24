function conj = conjunction(si1, si2, type, varargin)
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
% 
% ..
%    Zizhuang Miao, Oct. 2023
% .. 
% 
% :Inputs:
%
%   Added another input argument "type", allowed values are "indicator" or
% "values". "indicator" will return an output with voxel values set to 1 
% and -1 (as previously designed); "values" will return the average voxel
% values of the two input maps. Default to "values" to prevent crashing
% previous scripts using this function.
% 

direction = 0; % default: both pos and neg
if nargin == 2
    type = "values"; % default: average of two maps
end
if nargin > 3
    direction = varargin{1};
end

if compare_space(si1, si2), error('Objects are not in same space'); end

if ~ismember(type, {'indicator', 'values'}), error('The input type is not allowed'); end

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
pos_conj.dat(pos_conj.dat > 0) = 1; % set all pos values to 1
pos_conj.dat(pos_conj.dat < 0) = 0; % set all neg values to 0
pos_conj.sig = si1_pos.sig & si2_pos.sig;

% the neg conjuction
neg_conj = si1_neg;
neg_conj.dat(neg_conj.dat < 0) = -1; % set all neg values to -1
neg_conj.dat(neg_conj.dat > 0) = 0; % set all pos values to 0
neg_conj.sig = si1_neg.sig & si2_neg.sig;

if strcmp(type, 'indicator')    % output a [1 -1] map
    if direction==0 % both pos and neg
        conj = si1; % to load image space, etc.
        conj.dat = pos_conj.dat + neg_conj.dat;
        conj.sig = pos_conj.sig + neg_conj.sig;
    
    elseif direction < 0 % only neg
        conj = neg_conj;
    
    else % only pos
        conj = pos_conj;
    end
elseif strcmp(type, 'values')    % output a value map
    conj = si1; 
    conj.dat = (si1.dat + si2.dat)/2;    % averaging the two maps
    if direction==0 
        conj.sig = pos_conj.sig + neg_conj.sig;
    
    elseif direction < 0 % only neg
        conj.sig = neg_conj.sig;
    
    else % only pos
        conj.sig = pos_conj.sig;
    end
end
    
