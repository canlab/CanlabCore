function vif = getvif(model, no_add_intercept, varargin)

% function vif = getvif(model, [no_add_intercept], varargin)
%
% first varagin must be no_add_intercept.
% pass in 'wh' followed by vector to return the vifs for only those columns,
% though vif will still be calculated on entire model

if nargin < 2, no_add_intercept = 0; end

ind=strcmp('wh',varargin);
if ~isempty(ind), 
    wh=varargin{ind+1}; 
    if size(wh,1) > size(wh,2), wh=wh';end %transpose if needed
else
    wh = 1:size(model,2);    
end

vif = [];
for i = wh

    if no_add_intercept
        X = model;
    else
        X = [model ones(size(model,1),1)];
    end
    
    % check that model has an intercept, or that all preds are
    % mean-centered.  if both of these conditions are false, Rsquared is uninterpretable see [http://statistiksoftware.blogspot.com/2013/01/why-we-need-intercept.html] 
    % and thus VIF is uninterpretable.  Yoni, May 2015
    hasint = any(var(X)==0);
    totallymc = sum(mean(X))==0;
    if ~hasint & ~totallymc
        warning('Model has no intercept, and model is not mean-centered; VIF is not interpertable');
        vif = NaN;
        return
    end

     y = X(:,i);

     X(:,i) = [];

     b = X\y;
     fits = X * b;

     rsquare = var(fits) / var(y);

     % sometimes rounding error makes rsquare>1
     if rsquare >= 1, rsquare = .9999999; end

     vif(end+1) = 1 / (1 - rsquare);

end  
