
function a = toy2d(dist,hyper)

%====================================================================
% TOY2D toy data generation object
%====================================================================
% A=TOY(H) returns a toy object initialized with hyperparameters H.
%
% This generates 2d toy data problem of two Gaussians or from a 2d image.
%
% Hyperparameters, and their defaults
%  l=20          -- examples
%  seed=-1       -- random seed used to generate it, if -1 do not
%                   set seed
%  dist 		 -- name of distribution. can be the string: 'gaussians' [default]
%			        or must specify the name of a 'jpg' file (without extension) containing the colors
% 					red = class 1 , blue = class 2 , black
% Model
%  w             -- stores label model
%
% Methods:
%  generate,train,test
%
% Example:
% Possible distributions: gaussians,2circles,chess,rectangle_uneven,simple,uneven_gauss
% [r a]=train(svm('C=10'),toy2d); plot(a);
% pause
% clf,[r a]=train(svm({kernel('rbf',0.5),'C=10000'}),toy2d('2circles','l=200')); plot(a);

a.l=20;
a.n=2;
a.seed=-1;

if nargin >1
    a.dist=dist;
else
    a.dist='gaussians';
end

p=algorithm('toy2d');
a= class(a,'toy2d',p);

if nargin==2
    eval_hyper;
else
    if nargin==1
        hyper=dist;
        eval_hyper;
    end
end
