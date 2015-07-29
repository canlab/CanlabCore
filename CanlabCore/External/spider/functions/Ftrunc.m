function x = Ftrunc( x, prec)
  
% function x = Ftrunc( x, prec)
% 
% rounds x to prec significant digits
  
  f = 10.^fix( prec - log10( abs( x)));
  x = round( x.*f) ./ f;
  
