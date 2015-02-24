function Z = squareform(Y)
%SQUAREFORM Square matrix formatted distance. 
%   Z = squareform(Y) converts the output of PDIST function into a 
%   square format, so that Z(i,j) denotes the distance between the
%   i and j objects in the original data.
%
%   See also PDIST.

%   ZP You, 3-10-98
%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 1.7 $

[m, n] = size(Y);

if m ~= 1
   error('the input has to be a row vector');
end

m = (1 + sqrt(1+8*n))/2;

if m ~= fix(m)
   error('the size of the input is not correct');
end

Z = zeros(m);
I = ones(n,1);
J = [1 (m-1):-1:2]';
I(cumsum(J)) = (2:m);
Z(cumsum(I)) = Y';
Z = Z + Z';


