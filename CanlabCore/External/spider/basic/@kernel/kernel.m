function a = kernel(kType,param)   

%===========================================================================
% Kernel object - for calculating inner products in feature spaces
%===========================================================================
% Attributes: 
%  ker='linear'      -- type of kernel (linear, poly rbf, custom, SEE BELOW)
%  kerparam=[]       -- parameters of the kernel  
%  kercaching=0      -- yes if caching the kernel (only for data_global)
%  calc_on_output=0  -- calc kernel on outputs (Ys), rather than inputs (Xs)
%  output_distance=0 -- output associated distance rather than kernel
%  dat=[]            -- storage of data (for use when training/testing kernels)
% 
%
% Methods:
%  calc(k,d1,[d2]):     calc inner products in feature space between data d1 
%                       and d2 (or itself if d2 not specified) using kernel k 
%  get_norm(k,d):       calc norm of data in feature space
%  train,test
%   
% KERNEL               PARAMETERS &  DESCRIPTION
%  linear                             k(x,y)=x.y
%  poly                poly degree d, k(x,y)=(x.y+1)^d
%  rbf                 sigma,         k(x,y)=exp(-|x-y|^2/(2*sigma^2))
%  Gaussian            sigma,         k(x,y)=1/(2*pi^(N/2)*sqrt(sigma)) exp(..)
%  kmgraph             gamma,         marginalized kernel for graphs
%  spikernel                          kernel for spike trains
% 
%  weighted_linear     scale fact. w, k(x,y)=sum w_i^2 x_i y_i
%  weighted_poly       scale fact. w, k(x,y)=(sum w_i^2 x_i y_i+1)^d
%
%  rbf_of_dist         rbf kernel applied to an input distance matrix
%  poly_of_ker         poly kernel applied to an input kernel matrix
%
%  template_kernel     example to help make your own kernel 
%  custom              takes values from indices of matrix (kerparam)
%  custom_fast         takes values from global variable K 
%
% Example: - calc(kernel('rbf',2),data(rand(5)))
%          - calc(kernel('localfunctioninsearchpath',data(rand(5)))
%
% [Note: you can also access kernel params using "rbf",
% "poly" or "weighted_linear" for poly, rbf and weighted_linear kernels, which
% is implemented using alias member of @algorithm, aliasing the kerparam
% e.g a=svm(kernel('rbf')); a.rbf=4; %% set width of rbf kernel
% see help algorithm for more information]  
%================================================================================
% Reference : Learning with Kernels (Bernhard Schölkopf and Alexander J. Smola)
% Author    : Bernhard Schölkopf , Alexander J. Smola
% Link      : http://www.amazon.com/exec/obidos/tg/detail/-/0262194759/qid=1080825189/sr=1-1/ref=sr_1_1/002-6279399-2828812?v=glance&s=books
%================================================================================

a.calc_on_output=0;   %% work as an input kernel as default
a.output_distance=0;
a.callback=[];

a.kercaching = 0;
if nargin==0
    % If there is no args, linear is the default kernel 
    a.ker='linear'; 
    a.kerparam=[];
else
    if nargin==1,
        % if there is only one parameter, it may be because
        % this is a cell...
        if iscell(kType),
            a.ker = kType{1};
            for i=2:length(kType),
                a.kerparam{i-1}=kType{i};
            end
        else
            % ...or a custom kernel!
            a.ker=kType;
            a.kerparam=[];
        end
    else
        % There is more than one arg
        a.ker = kType;
        a.kerparam = param;
    end
end

a.dat=[];
algo=algorithm('kernel');
a= class(a,'kernel',algo);

if strcmp(a.ker,'rbf_of_kernel') 
    a.algorithm.alias={'rbf','kerparam'};
end;

if strcmp(a.ker,'rbf') 
    a.algorithm.alias={'rbf','kerparam'};
end;

if strcmp(a.ker,'poly') 
    a.algorithm.alias={'poly','kerparam'};
end;
if strcmp(a.ker,'weighted_linear') 
    a.algorithm.alias={'weighted_linear','kerparam'};
end;  




