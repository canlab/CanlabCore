function a = distance(c,param)   

%=============================================================================
% Distance object - for calculating inner products in feature spaces
%=============================================================================
% Attributes: 
%  dist=[];          -- type of distance (SEE BELOW)
%  child=kernel;     -- if nonempty, generate distance from this kernel 
%                       rather than from type described in dist
%  distparam=[]      -- parameters of the distance  
%  kercaching=0      -- yes if caching the distance (only for data_global)
%  calc_on_output=0  -- calc distance on outputs (Ys), rather than inputs (Xs)
%  output_kernel=0   -- output associated kernel (if any) rather than distance
%  dat=[]            -- storage of data (for use when training/testing dists)
% 
%
%
% Methods:
%  calc(d,d1,[d2]):    calc inner products in feature space between data d1 
%                      and d2 (or itself if d2 not specified) using distance d 
%  train,test
%   
% DISTANCE               PARAMETERS & DESCRIPTION
%  euclid                d(x,y)=||x-y||_2 (same as using child=kernel)
%  norm                  d(x,y)=||x-y||_p, p can be inf, -inf see help norm 
%
%  template_distance   example to help make your own distance 
%  custom              takes values from indices of matrix (distparam)
%  custom_fast         takes values from global variable K 
%
% Example: d=data(rand(5));
%          calc(distance,d)                     - use 
%          calc(distance(kernel('rbf',2)),d)
%          calc(distance('norm',1),d)
%          calc(distance('euclid'),d)
% 
%=============================================================================
% Reference : 
% Author    : 
% Link      : 
%=============================================================================




a.calc_on_output=0;   %% work as an input distance as default
a.output_distance=0;
a.child=[]; 
a.kercaching = 0;
a.dist=[]; 
a.distparam=[];

if nargin==0
    % If there is no args, linear is the default distance 
    a.dist=[]; 
    a.distparam=[];
    a.child=kernel; 
else
    if nargin==1,
        % if there is only one parameter, it may be because
        % this is a cell...
        if iscell(c),
            a.dist = c{1};
            for i=2:length(c),
                a.distparam{i-1}=c{i};
            end
        else
            % ...or a custom distance!
            if( isstr(c))
                a.dist=c;
                a.distparam=[];
            else
                % we assume c is a kernel
                a.child=c;
            end
        end
    else
        % There is more than one arg
        a.dist = c;
        a.distparam = param;
    end
end

a.dat=[];
p=algorithm('distance');
a= class(a,'distance',p);

if strcmp(a.dist,'weighted_linear') a.algorithm.alias={'weighted_linear','distparam'};end;  




