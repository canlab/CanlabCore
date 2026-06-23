function CorrLeng=AutoCorrLength(ts,T)
% function CorrLeng=AutoCorrLength(TS,T)
% Calculates the Correlation Lengths. i.e. measure how bad the things are
% in terms of autocorrelation. Adapted from:
% 
%       Straatsma, T. P., Berendsen, H. J. C., & Stam, A. J. (2016). 
%       Estimation of statistical errors in molecular simulation calculations,
%       8976(June). http://doi.org/10.1080/00268978600100071
%       
%       Consideration: In the original definition, the correlation length
%       is supposed to consider all lags (see Eq. 4), however, thi function 
%       offer curbing as lots of very far lags, in case of BOLD time
%       series, are merely crappy estimation of zero. 
%
%%%REFERENCES:
%  Afyouni, Soroosh, Stephen M. Smith, and Thomas E. Nichols. 
% "Effective Degrees of Freedom of the Pearson's Correlation Coefficient 
%  under Serial Correlation." bioRxiv (2018): 453795.
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2017
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

CorrLeng = sum(AC_fft(ts,T).^2);

end
