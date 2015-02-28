function weights = get_w(a) 

%
%               w = get_w() 
% 
%   Returns the weight vector from the svm 
%   AE: i have changed this procedure to include non linear feature
%   selection with RFE. When non-linear, it gives  a score vector for each
%   feature (used in RFE) which is not the margin but is sorted as the margin

weights=a.w;
