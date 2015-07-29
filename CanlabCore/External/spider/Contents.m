% THE SPIDER Version 1.71 (24/7/2006)  -- Current Repository
%
% Basic library objects. 
%   data         - Storing input data and output results 
%   data_global  - Implementation of data object that limits memory overhead
%   algorithm    - Generic algorithm object
%   group        - Groups sets of objects together (algorithms or data) 
%   loss         - Evaluates loss functions
%   get_mean     - Takes mean loss over groups of algs
%   chain        - Builds chains of objects: output of one to input of another
%   param        - To train and test different hyperparameters of an object
%   cv           - Cross validation using objects given data
%   kernel       - Evaluates and caches kernel functions
%   distance     - Evaluates and caches distance functions
%
% Statistical Tests objects.
%   wilcoxon     - Wilcoxon test of statistical significance of results
%   corrt_test   - Corrected resampled t-test - for dependent trials
%
% Dataset objects.
%   spiral       - Spiral dataset generator.
%   toy          - Generator of dataset with only a few relevant features
%   toy2d        - Simple 2d Gaussian problem generator
%   toyreg       - Linear Regression with o outputs and n inputs 
%
% Pre-Processing objects
%   normalize    - Simple normalization of data
%   map          - General user specified mapping function of data
%
% Density Estimation objects.
%   parzen       - Parzen's windows kernel density estimator
%   indep        - Density estimator which assumes feature independence
%   bayes        - Classifer based on density estimation for each class
%   gauss        - Normal distribution density estimator
%                    
% Pattern Recognition objects.
%   svm          - Support Vector Machine (svm)
%   c45          - C4.5 for binary or multi-class 
%   knn          - k-nearest neighbours
%   platt        - Conditional Probability estimation for margin classifiers
%   mksvm        - Multi-Kernel LP-SVM
%   anorm        - Minimize the a-norm in alpha space using kernels
%   lgcz         - Local and Global Consistent Learner 
%   bagging	     - Bagging Classifier
%   adaboost     - ADABoost method
%   hmm          - Hidden Markov Model 
%   loom         - Leave One Out Machine 
%   l1           - Minimize l1 norm of w for a linear separator 
%   kde          - Kernel Dependency Estimation: general input/output machine
%   dualperceptron       - Kernel Perceptron
%   ord_reg_perceptron   - Ordinal Regression Perceptron (Shen et al.)
%   splitting_perceptron - Splitting Perceptron (Shen et al.)
%   budget_perceptron    - Sparse, online Pereceptron (Crammer et al.)
%   randomforest - Random Forest Decision Trees         WEKA-Required
%   j48          - J48 Decision Trees for binary        WEKA-Required
%
% Multi-Class and Multi-label objects. 
%   one_vs_rest  - Voting method of one against the rest (also for multi-label)
%   one_vs_one   - Voting method of one against one
%   mc_svm       - Multi-class Support Vector Machine by J.Weston
%   c45          - C4.5 for binary or multi-class 
%   knn          - k-nearest neighbours
%               
% Feature Selection objects.
%   feat_sel     - Generic object for feature selection + classifier
%   r2w2_sel     - SVM Bound-based feature selection
%   rfe          - Recursive Feature Elimination (also for the non-linear case)
%   l0           - Dual zero-norm minimization (Weston, Elisseeff)
%   fsv          - Primal zero-norm based feature selection (Mangasarian)
%   fisher       - Fisher criterion feature selection
%   mars         - selection algorithm of Friedman (greedy selection)
%   clustub      - Multi-class feature selection using spectral clustering
%   mutinf       - Mutual Information for feature selection.
%      
% Regression objects.
%   svr          - Support Vector Regression
%   gproc        - Gaussian Process Regression 
%   relvm_r      - Relevance vector machine 
%   multi_rr     - (possibly multi-dimensional) ridge regression   
%   mrs          - Multivariate Regression via Stiefel Constraints      
%   knn          - k-nearest neighbours
%   multi_reg    - meta method for independent multiple output regression
%   kmp          - kernel matching pursuit
%   kpls         - kernel partial least squares
%   lms          - least mean squared regression [now obselete due to multi_rr]
%   rbfnet       - Radial Basis Function Network (with moving centers)
%   reptree      - Reduced Error Pruning Tree       WEKA-Required
%
% Model Selection objects.
%   gridsel      - select parameters from a grid of values 
%   r2w2_sel     - Selecting SVM parameters by generalization bound
%   bayessel     - Bayessian parameter selection 
%
% Unsupervised objects.
%   one_class_svm - One class SVM
%   kmeans       - K means clustering
%   kvq          - Kernel Vector Quantization
%   kpca         - Kernel Principal Components Analysis
%   ppca         - Probabilistic Principal Component Analysis
%   nmf          - Non-negative Matrix factorization
%   spectral     - Spectral clustering
%   mrank        - Manifold ranking
%   ppca         - Probabilistic PCA
%
% Reduced Set and Pre-Image objects.
%   pmg_mds      - Calculate Pre-Images based on multi-dimensional scaling
%   pmg_rr       - Calculate Pre-Images based on learning and ridge regression
%   rsc_burges   - Bottom Up Reduced Set; calculates reduced set based on gradient descent
%   rsc_fp       - Bottom Up Reduced Set; calculates reduced set for rbf with fixed-point iteration schemes
%   rsc_mds      - Top Down Reduced Set; calculates reduced set with multi-dimensional scaling
%   rsc_learn    - Top Down Reduced Set; calculates reduced set with ridge regression
%   rss_l1       - Reduced Set Selection via L1 penalization
%   rss_l0       - Reduced Set Selection via L0 penalization
%   rss_mp       - Reduced Set Selection via matching pursuit
%
% See also help on: demos, train, test.
% Most of these algorithms are available as extra download on the website:
% http://www.kyb.tuebingen.mpg.de/bs/people/spider/index.html 

