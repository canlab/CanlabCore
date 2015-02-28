function a = kmeans(hyper) 

%====================================================================
% KMEANS k-means clustering object
%====================================================================
% A=KMEANS(H) returns a k-means clustering object
% initialized with hyperparameters H.
%
% Hyperparameters and their defaults:
%   k=2              -- the number of clusters
%   a.child=distance -- distance measure to use (a distance object)
%   max_loops=1000   -- maximum number of iterations of training%
% Model:
%   mu               -- cluster centers%   y                -- cluster assignment of each training example%
% Methods:
%   training: cluster a dataset
%   testing:  assign points to clusters according to the nearest
%             cluster center
%   distortion: get value of distortion%
% Example: 
%  d=gen(spiral({'m=200','n=0.5','noise=0.35'}));
%  [r,a] =train(kmeans,d)
%  I=find(r.X==1);clf;hold on;
%  plot(d.X(I,1),d.X(I,2),'r.');
%  I=find(r.X==2);
%  plot(d.X(I,1),d.X(I,2),'b.');
% ['compare with spectral clustering ']%
%====================================================================
% Reference : chapter 10 (Richard O. Duda and Peter E. Hart) Unsupervised learning and clustering
% Author    : Richard O. Duda , Peter E. Hart
% Link      : http://www.amazon.com/exec/obidos/tg/detail/-/0471056693/002-6279399-2828812?v=glance
%====================================================================
  
  a.k=2;
  a.child=distance;
  a.mu=[]; a.y=[]; 
  a.distortion=0;  a.max_loops=1000;

  p = algorithm('kmeans');
  a = class(a,'kmeans',p);
  
  if nargin==1,
    eval_hyper;
  end;
