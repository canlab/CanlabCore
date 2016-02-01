function goodness = classify_naive_bayes_objfun(dat, Xi, volInfo, a, s, g, k, t, h)
% :Usage:
% ::
%
%     goodness = classify_naive_bayes_objfun(dat, Xi, volInfo, 0, .9, 1, .05, .1, .1);
%
% :Examples:
% ::
%
%    objfun = @(t, h) classify_naive_bayes_objfun(dat, Xi, volInfo, 0, .9, 1, .05, t, h);
%    goodness = objfun(.1, .1)

bayes_model = classify_naive_bayes('setup', dat, Xi, a, s, g, k);

doplot = 0;
reducedY = bayes_meta_feature_abstract(dat, t, h, bayes_model, volInfo, doplot);

xval = classify_naive_bayes('xval', reducedY, Xi, 0, 0, g, k );

goodness = mean(xval.prop_correct_by_class);

end

