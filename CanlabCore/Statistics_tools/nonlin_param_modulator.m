function [h_mean h_by_mod d_mean d_by_mod intcpt auc_mean_trial auc_by_modulator errval] =  nonlin_param_modulator(y, ons, modulator_centered, tr, xvals)
% :Usage:
% ::
%
%     [h_mean h_by_mod d_mean d_by_mod intcpt auc_mean_trial auc_by_modulator errval] =  nonlin_param_modulator(y, ons, modulator_centered, tr, xvals)
%
% USE THIS with nonlin_param_mod_brain.m
%
% This function takes data (y) and other things (onsets, etc.) and
% produces parameter estimates for amplitude, duration, amp*modulator,
% dur*modulator
% 
% See rt_fit_brain.m for more info, and for examples creating fit
% plots, etc.
%
% :Examples:
% ::
%
%    y = rand(360, 1);
%    [h_mean h_by_mod d_mean d_by_mod intcpt auc_mean_trial auc_by_modulator errval] =  nonlin_param_modulator(y, ons, scale(RTs, 1), 2);
%    nonlin_parammod_predfun(y, ons, scale(RTs), [h_mean h_by_mod d_mean d_by_mod intcpt], 'plot');
%    hold on; plot(y, 'k')

%    % THIS stuff is the same for each voxel
%   RTcenter = RTs - mean(RTs);
%   xvals = (1:length(y))';
%
%  % fitting function: times is a dummy var to get this to work with nonlin_fit
%   fhan = @(p, times) nonlin_parammod_predfun(y,ons,RTcenter,p);
%
%   fitting_fun = @(y) nonlin_fit(y, xvals, 'link', fhan, 'start',[1 1 1 1 mean(y)]);
%
% ..
%    Tor Wager, May 2008
% ..


    if nargin < 5 || isempty(xvals)
        xvals = (1:length(y))';
    end
    
    if nargin < 4 || isempty(tr)
        disp('Using default TR of 1 s');
        tr = 1;
    end
    
    
       % fitting function: times is a dummy var to
    % get this to work with nonlin_fit
    fhan = @(p, times) nonlin_parammod_predfun(y, ons, modulator_centered, p, 'tr', tr);

    % This gives you param estimates, given a data vector
    %fitting_fun = @(y) nonlin_fit(y, xvals, 'link', fhan, 'start',[1 0 1 0 mean(y)]);

    % This stuff is different for each voxel
    [p, errval] = nonlin_fit(y, xvals, 'link', fhan, 'start',[1 0 1 0 mean(y)]);  %fitting_fun(y);

    h_mean = p(1);
    h_by_mod = p(2);
    d_mean = p(3);
    d_by_mod = p(4);

    intcpt = p(5);

    auc_mean_trial = p(1) * p(3);  % height * width = area of average response

    % with unit param. modulator increase, what is change in AUC?
    % Height(PM = 1) * Width(PM = 1) - auc_mean_trial

    auc_by_modulator = (p(1) + p(2) .* 1) * (p(3) + p(4) .* 1) - auc_mean_trial;


end



% 
% true_p = [1 1.5 3 1.5 0];          % true parameters
% y = fhan(true_p);                  % generate simulated "true" signal
% [p, errval, fit] = fitting_fun(y);
% true_p
% p
% create_figure('test'); plot(y); hold on; plot(fit, 'r');

%p = lsqcurvefit(fhan,[1 1 1 1 .5], xvals, y) 

% p = params
% 1) mag. scale (amplitude)
% 2) mag. X RT slope (linear effect of RT)
% 3) duration intercept (neural epoch duration)
% 4) duration x RT slope (linear effect of RT on duration)
% 5) overall fitted response intercept



