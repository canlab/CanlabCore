% -----------------------------------------------------------------------------------------------------------------------------
% Author: Guilherme V. Rocha
%         Department of Statistics
%         University of California, Berkeley
%         gvrocha@stat.berkeley.edu, gvrocha@gmail.com
% 2006/09
% -----------------------------------------------------------------------------------------------------------------------------
% This script provides an example of using the lasso routines for matlab
% The data set used is the same as the one used in Efron et. al (2004) ``Least angle regression'' paper
% The ``stress'' data sets are modified versions of this data set (adding noise variables and LD columns to test the algorithm)

clear variables;

addpath('./../');
addpath('./../../');

load matlab_diabetes_data.txt;
x = matlab_diabetes_data(:,1:10);
y = matlab_diabetes_data(:,11);

options.trace = 1;

res_ls       = lasso(y, x, options);
sizes_ls     = abs(res_ls.nbeta)*ones(size(res_ls.nbeta,2),1);
rel_sizes_ls = sizes_ls/max(sizes_ls);

cv_res = lasso_cv(y, x, 10);

load stress_diabetes_data_v1.txt;
x = stress_diabetes_data_v1(:,1:end-1);
y = stress_diabetes_data_v1(:,end);

stress_ls           = lasso(y, x, options);
stress_sizes_ls     = abs(stress_ls.nbeta)*ones(size(stress_ls.nbeta,2),1);
stress_rel_sizes_ls = stress_sizes_ls/max(stress_sizes_ls);

load stress_diabetes_data_v2.txt;
x = stress_diabetes_data_v2(:,1:end-1);
y = stress_diabetes_data_v2(:,end);

stress_ls_v2           = lasso(y, x, options);
stress_sizes_ls_v2     = abs(stress_ls_v2.nbeta)*ones(size(stress_ls_v2.nbeta,2),1);
stress_rel_sizes_ls_v2 = stress_sizes_ls_v2/max(stress_sizes_ls_v2);

axis_limits = [0 1 -max([abs(res_ls.nbeta(end,:)) abs(stress_ls.nbetas(end,:)) abs(stress_ls_v2.nbetas(end,:))]) max([abs(res_ls.nbetas(end,:)) abs(stress_ls.nbetas(end,:)) abs(stress_ls_v2.nbetas(end,:))])];
figure(1);plot(rel_sizes_ls, res_ls.nbeta);axis(axis_limits);
figure(2);plot(rel_sizes_ls, res_ls.nbeta);axis(axis_limits);
figure(3);plot(stress_rel_sizes_ls_v2, stress_ls_v2.nbeta);axis(axis_limits);
