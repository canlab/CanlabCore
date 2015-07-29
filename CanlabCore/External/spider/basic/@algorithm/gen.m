function dat =  gen(algo)
%
%           data =  gen(algorithm)
%
%  Data will be generated from a generative model. This should be
%  equivalent to 
%                    test(algo)
%  (no data is supplied to the test function) 
  
  dat = generate(algo);