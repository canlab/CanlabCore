function [V, S, varexp, w, Yhat] = plssquash(X, Y, varargin)
% Decomposes data X into K components that are ordered in their 
% covariance with Y, designed to predict orthogonal parts of Y
%
% :Usage:
% ::
%
%     [V, S, varexp, w, Yhat] = plssquash(X, Y, varargin)
%
% :Optional Inputs:
%   - case {'noplot'}, turn off plotting
%   - case 'ndims', save only first ndims (K) vectors
%
% :Outputs:
%
%   **V:**
%        'eigenvetors', or weights, on data (columns)
%
%   **S:**
%        score matrix, N x K
%
%   **varexp:**
%        sqrt(r-square) with first k components predicting Y
%
%   **w:**
%        V*b, integrated weights.  for predicting new data, pred = X*w
%
%   **Yhat:**
%        X*V*b, or S*b
%
% ..
%    Tor Wager, 9/12/09
% ..

doplot = 1;
ndims = length(Y) - 1; 

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'noplot'}, doplot = 0;
            case 'ndims', ndims = varargin{i + 1};
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end
    
% univariate covariance-based weights
% imperfect prediction
X = X - repmat(mean(X), size(X, 1), 1); %scale(X);
Y = scale(Y);

%
% neither V nor S appears to have orthogonal columns, though many are
% later columns, as df in Y is approached, appear to be highly colinear
% this may be because all variance in Y is basically explained...

clear V  % V is voxel weights for each component 1:K
clear S
rY = Y; % Y values to successively predict; intialize to Y

for i = 1:ndims  %ceil(length(Y)./2)

    % V = voxel weights, based on univariate relationship with rY
    V(:,i) = (X' * rY); %.^ .5;
    V(:,i) = V(:,i) ./ norm(V(:,i));

    % score matrix, N x K
    % these will ultimately be predictors for Y
    % chosen to maximize predictive power with few orthogonal components
    S = X*V;

    rY = Y - S * pinv(S) * Y;

    b = pinv(S) * Y;
    Yhat = S * b;

    rsq(i) = 1 - var(rY) ./ var(Y);
    
    if doplot
        create_figure('Fit'); plot(Yhat, Y, 'kx'); refline; 

        rsq_adj(i) = 1 - (1-rsq(i)) * ((size(Y,1)-1) ./ (size(Y,1)-i-1));

        title(sprintf('cum. r-squared = %3.1f, adj = %3.1f', 100*rsq(i), 100*rsq_adj(i)));
        drawnow; pause(.3)
    end
    
end

varexp = sqrt(rsq);

b = pinv(S) * Y;
Yhat = S * b;
w = V * b;  % final weight vector

%prediction = mean(training) + data*S*b
%Yhat = X * V * b;
%figure; plot(Yhat, Y, 'kx')

end
