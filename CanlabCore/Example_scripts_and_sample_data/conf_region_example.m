% This script contains an example dataset for multivariate confidence
% intervals from Johnson & Wichern, p. 237, ex. 5.3, and demonstrates the
% use of the function conf_region to estimate and plot the confidence
% region.

% Johnson & Wichern, p. 237, ex. 5.3
% multivariate confindence regions
% data is microwave radiation, door closed (x1) and open (x2)
% 4th root of measured radiation is the var. of interest.

x = [    0.1500    0.3000
    0.0900    0.0900
    0.1800    0.3000
    0.1000    0.1000
    0.0500    0.1000
    0.1200    0.1200
    0.0800    0.0900
    0.0500    0.1000
    0.0800    0.0900
    0.1000    0.1000
    0.0700    0.0700
    0.0200    0.0500
    0.0100    0.0100
    0.1000    0.4500
    0.1000    0.1200
    0.1000    0.2000
    0.0200    0.0400
    0.1000    0.1000
    0.0100    0.0100
    0.4000    0.6000
    0.1000    0.1200
    0.0500    0.1000
    0.0300    0.0500
    0.0500    0.0500
    0.1500    0.1500
    0.1000    0.3000
    0.1500    0.1500
    0.0900    0.0900
    0.0800    0.0900
    0.1800    0.2800
    0.1000    0.1000
    0.2000    0.1000
    0.1100    0.1000
    0.3000    0.3000
    0.0200    0.1200
    0.2000    0.2500
    0.2000    0.2000
    0.3000    0.4000
    0.3000    0.3300
    0.4000    0.3200
    0.3000    0.1200
    0.0500    0.1200];

X = x.^.25;

figure;
[ci,ax,S,e,lam,Fm,Fc] = conf_region(X, 1);