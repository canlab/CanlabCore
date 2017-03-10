
function a = map(hyper) 

%=======================================================================
%  Mapping object
%=======================================================================  
% A=MAP(H) returns a map object initialized with hyperparameters H. 
%
%  MAP is a general object for mapping data, it takes a supplied
%  mapping function which can transform the given input data. The 
%  function can use four hyperparmeters p1,p2,p3,p4.
%
%  Hyperparameters, and their defaults
% 
%   func='d.X=tanh((d.X+a.p1)*a.p2)' -- function we wish to use.
%                                     This can be changed at will,
%                                     it can access the object 'd' (data)
%                                     and 'a' (mapping algorithm), e.g
%                                     d.X, d.Y, a.p1, a.p2, etc. and 
%                                     usual matlab functions. It should 
%                                     store the new values in d.X. 
%                                      The default is a sigmoid mapping.
%   p1=1,p2=0,p3=[],p4=[]          -- user parameters
%=======================================================================
% Reference : 
% Author    : 
% Link      : 
%=======================================================================
    
  a.func='x=tanh((d.X+a.p2)*a.p1)';        
  a.p1=1; a.p2=0; a.p3=[]; a.p4=[];
  
  p=algorithm('map');
  a= class(a,'map',p);
  
  if nargin==1
    eval_hyper;
  end  



