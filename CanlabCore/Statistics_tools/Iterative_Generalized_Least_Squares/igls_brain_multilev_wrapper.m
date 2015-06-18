% [a, b, c1, c, ab, aste, bste, c1ste, cste, abste, ap, bp, c1p, cp, abp, ...
% aind, bind, c1ind, cind, abind, aiste, biste, c1iste, ciste, abiste] = ...
% igls_brain_multilev_wrapper(dmpfc, hr, pag, 'boot');
%
% Wrapper function for mediation_brain_multilev
% Returns separate outputs for each variable that deserves a separate image name
%
% Tor Wager, Oct 2007
%
% Outputs:
% % names = {'intercept_b.img' 'slope_b.img' 'intercept_rfxvar_b.img' 'slope_rfxvar_b.img' ...
% %     'intercept_indiv.img' 'slope_indiv.img' ...
% %     'sigma.img' 'isconverged.img' 'intercept_t.img' 'slope_t.img'  ...
% %     'intercept_rfxvar_t.img' 'slope_rfxvar_t.img' 'intercept_p.img' 'slope_p.img'  ...
% %     'intercept_rfxvar_p.img' 'slope_rfxvar_p.img' };
%
% All this function does is run igls.m on whatever you put in, and then
% separates the output from a structure into separate output variables,
% which will then be written as separate image files.

function [intercept_b, slope_b, intercept_rfxvar_b, slope_rfxvar_b, intercept_indiv, slope_indiv  ...
        sigma, isconverged, intercept_t, slope_t, intercept_rfxvar_t, slope_rfxvar_t, ...
        intercept_p, slope_p, intercept_rfxvar_p, slope_rfxvar_p, intercept_LRTrfxvar_p, slope_LRTrfxvar_p, ...
        cov_int_b cov_int_t cov_int_p cov_slope_b cov_slope_t cov_slope_p] = ...
        igls_brain_multilev_wrapper(X, Y, varargin)

    % main function call
    out = igls(Y, X, varargin{:});
    
    % parse outputs into separate variables, for writing to separate images
    intercept_b = out.beta(1);
    slope_b = out.beta(2);
    
    intercept_rfxvar_b = out.betastar(1);
    slope_rfxvar_b = out.betastar(2);
    
    intercept_indiv = out.beta_indiv(1, :)';
    
    slope_indiv = zeros(size(intercept_indiv));
    %slope_indiv = out.beta_indiv(2, :)';
    % NOTE! NOT WRITTEN CORRECTLY!  FIX ME
    
    sigma = out.Sigma;
    
    isconverged = out.isconverged;
    
    intercept_t = out.t(1);
    slope_t = out.t(2);
    
    % NOTE: the t-values / p-values for rfx are  different from the likelihood ratio test (LRT)!
    % the LRT p is the preferred one
    intercept_rfxvar_t = out.t_randvariance(1);
    slope_rfxvar_t = out.t_randvariance(2);
    
    intercept_p = out.p(1);
    slope_p = out.p(2);
    
    intercept_rfxvar_p = out.p_randvariance(1);
    slope_rfxvar_p = out.p_randvariance(2);
    
    intercept_LRTrfxvar_p = out.pLRT_randvariance(1);
    slope_LRTrfxvar_p = out.pLRT_randvariance(2);

    % Outputs for the btwn-subjects covariate
    % igls Works now for only ONE covariate!
    % Col 3 is cov*1st level intercept, col 4 is cov*1st level slope
    if length(out.beta) > 2
        cov_int_b = out.beta(3);
        cov_int_t = out.t(3);
        cov_int_p = out.p(3);
        
        cov_slope_b = out.beta(4);
        cov_slope_t = out.t(4);
        cov_slope_p = out.p(4);

    else
        cov_int_b = 0;
        cov_int_t = 0;
        cov_int_p = 1;
        
        cov_slope_b = 0;
        cov_slope_t = 0;
        cov_slope_p = 1;
    end
    
end

% % 
% % 
% %             analysisname: 'IGLS: Iterative generalized least squares analysis'
% %                     beta: [4x1 double]
% %                 betastar: [2x1 double]
% %               beta_indiv: [0.2120 0.4492 0.4888 -0.2822 0.3759 -0.1493 1.1521 0.4335 -0.1680 0.1057]
% %               beta_names: {'Intcpt.'  'Slope1'  'Covariate_intcpt'  'Covariate_slope'}
% %                 Cov_beta: [4x4 double]
% %             Cov_betastar: [2x2 double]
% %                    Sigma: 0.2360
% %                      Phi: []
% %                     type: 'i'
% %                  arorder: 0
% %              isconverged: 1
% %                  num_obs: 100
% %                      sub: 10
% %                 num_iter: 5
% %                  epsilon: 0.0100
% %               iterations: 1
% %             elapsed_time: 0.6617
% %             inputOptions: [1x1 struct]
% %                      ste: [4x1 double]
% %                        t: [4x1 double]
% %                  df_beta: 9
% %                        p: [4x1 double]
% %                  p_tails: 'two-tailed'
% %           t_randvariance: [2x1 double]
% %              df_betastar: 9
% %           p_randvariance: [2x1 double]
% %     p_randvariance_tails: 'one-tailed'
% %                      LRT: [2x1 double]
% %        pLRT_randvariance: [2x1 double]
% %        