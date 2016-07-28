function A = sdt_A(H, F)
% Calculate signal detection theory A statistic: 
% a non-parametric estimate of sensitivity in ROC analysis
% per Zhang and Mueller (2005)
%
% :Usage:
% ::
%
%     A = sdt_A(H, F)
%
% :Inputs:
%
%   **H:**
%        hit rate
%
%   **F:**
%        false alarm rate
%
% :Outputs:
%
%   **A:**
%        A measure of sensitivity of a detector, essentially the area under
%        the ROC curve for either "normal" (H>=F) or "reverse skill" (F>H)
%        case (ranges from 0.5 to 1).
%
% Meaning of values where F > H ("reverse skill"): 
%
% The meaningful measure of sensitivity in this case is A(F,H) 
% as opposed to A(H,F).  Subject is assumed to still be sensitive, 
% but to have inverted their interpretation (or "reverse skill"). 
% 
% To distinguish from the normal case, such results are returned
% as (1 - sdt_A(F,H)), so the full range from 0 to 1 is used.
% For comparisons of sensitivity use abs(sdt_A(H,F) - 0.5).
%
% :Example:
%
% sdt_A(0.6,0.8) = 0.325. This indicates the same sensitivity as 
% sdt_A(0.8,0.6) = 0.675, but with reversed skill. 
%
% ..
%    Author: Joe Wielgosz 8/13/2009
% ..


if F <= H % Normal case

    if (F <= 0.5) && (0.5 <= H)
        A = 0.75 + (H-F)/4 - F*(1-H);
    elseif (F <= H) && (H <= 0.5)
        A = 0.75 + (H-F)/4 - F/(4*H);
    else % (0.5 < F) && (F <= H)
        A = 0.75 + (H-F)/4 - (1-H)/(4*(1-F));
    end
    
else % "Reverse skill" - more false positives than hits
    % Using Jason Buhle's approach
    A = 1 - sdt_A(F, H);
%    A = NaN;
end
    
end
