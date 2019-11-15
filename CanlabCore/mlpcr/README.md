# Multilevel Principle Component Regression (MLPCR) Toolbox
## Created by Bogdan Petre, 2019
### Bogdan.Petre.GR@dartmouth.edu

mlpcr.m
The top level script. This script has hyperparameters which need to be fit 
somehow. User must decide which objective function to optimize and how to 
optimize it to obtain hyperparameter fits. Bayesian hyperparameter optimization 
works well for this (see bayesopt() documentation in Matlab r2016b or later). 
mlpcr_full() implements hyperparameter optimization and MLPCR model fitting 
using some sensible default choices. See help documentation for example use.
Usage is the same for mlpcr_full, except you prefix a couple extra arguments
related to hyperparameter optimization.

mlpcr_full.m 
Fits an MLPCR model using bayesian hyperparameter optimization with 
a MSE objective function for determining optimal hyperparameters to use. The 
user must specify which mlpcr() options are hyperparameters, but either pca 
dimensions alone or dimensions and covariance patterns might be suitable 
choices. mlpcr_full() provides a quick and easy way to implement mlpcr(), but 
the user should think carefully about whether or not MSE is an appropriate 
objective function for their task. Other choices might be to only consider 
within subject error, and discard between subject error, or to consider between 
subject error but discard intercepts depending on what the analyst believes is 
a meaningful measure in their outcome data. Consider updating mlpcr_cv_pred.m 
STATS object to return any objective function metrics that you find to be 
useful.

mlpcr_cv_pred.m 
Performs CV according to user specified CV folds and returns out of fold 
predictions. With an appropriately designed wrapper function defining the
objective function to optimize, mlpcr_cv_pred() could be useful for custom
implementations of hyperparameter optimization schemes if the defaults 
(bayesopt() with a MSE objective function) implemented in mlpcr_full()) are not
desirable.

multithreadWorkers.m
Automatically asigns threads to workers. Useful when fewer workers are requested 
than threads. For instance, if doing 5-fold cross validation on a 16 core 
machine multithreadWorkers() will assign 3 threads to each worker (by default 
workers are single threaded in matlab). This can speed up some operations, 
especially matrix math, and affects pca performance. multithreadWorkers() is 
potentially useful in many applications, not just when using the mlpcr toolbox.
