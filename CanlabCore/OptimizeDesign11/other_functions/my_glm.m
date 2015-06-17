function [t,contrastse] = my_glm(x,y,c)

    %function t = my_glm(x,y,c)

    % t is the resulting T score of the comparison

    % x is the design matrix

    % c is a vector of contrasts

    % y is the vector of data

    % Modified 2/13/00 by Tor Wager
    % Modified 3/17/01 by Tor to fix contrastse.
    % tests only the first contrast!
    % Modified 5/24/01 by Tor Wager

    x = [x ones(size(x,1) , 1)];

    c = [c; zeros(1,size(c,2))];

    if size(c, 1) > size(x, 2)
        c = c(1:size(x, 2), :);
    end

    % testing:
    %x =[x(:,1) ones(size(x,1),1)];
    %c = [1;0];

    p = size(x,2);

    n = size(x,1);



    % this stuff looks wrong to me...
    % xtx_inv = pinv(x'*x);

    % tested 5/24/01:  pinv(a) = inv(a'*a)*a'
    % beta_est = xtx_inv*x'*y;	   % this would be ok, if xtx_inv were right.


    % ...replaced it with this:

    xtx_inv = inv(x'*x);
    beta_est = x \ y;

    RSS = y - x*beta_est;

    RSS = sum(RSS.^2);



    var_est = RSS/(n-p);

    v = var_est*diag((c' * xtx_inv * c));	% variance of beta or contrast

    t = (beta_est' * c) ./ (sqrt(v)');


    contrastse = sqrt(v);

    % se of contrast


    % testing stuff
    % ------------------------------------------------------
    % 05/24/01
    % this is wrong because:
    % 	doesn't account for smoothing and autocorrelation, and it has to.
    %	positive bias in t-scores for null data!
    %
    % I think the variance formulas are right.
    % What's effective df?
    % I think it's b/c predictors are correlated, and no account is made of that.
    % If predictors are correlated, these results are wrong!
    % ALSO: smoothing positively biases the t-values, not negatively...why?
    % ------------------------------------------------------
    %figure;subplot(3,1,1);  %imagesc(x)
    %plot(x*beta_est,'LineWidth',2);hold on;plot(y,'r')
    %subplot(3,1,2);plot(x(:,1)*beta_est(1));hold on;plot(y,'r')
    %subplot(3,1,3);plot(x(:,2)*beta_est(2));hold on;plot(y,'r')
    %beta_est
    %beta_est2 = xtx_inv*x'*y
    %corrcoef(x)
    %xtx = x'*x
    %xtxi = inv(xtx)
    %xtx_inv
    %cxtxicprime = (c' * xtx_inv * c)
    %var_est
    %t
    %input('Press return to continue...')
    %close





    return