function dw = wekaNumericalData(varargin),
% DataWeka = WEKANUMERICALDATA(d,[xnames]);
% 
% Transform a spider data object into a weka instances object (with
% numerical target)
%
% By default everything is numerical except the target which can
% be nominal if the number of values is limited.
%
% Note that the dimension of the output is equal to 1. The first dimension
% of the output will be considered when it's greater than 1.
%
% Parameters
%
%   d           -- spider data object
%   xnames      -- optional list of string for the names of the input
%               features
%   ignoretarget -- a binary flag that indicates if the Y field of data
%                   object should be included.
%
% Output
% 
%   inst        -- weka instances object where the target has been set.
%               

xnames=[];
if (length(varargin)>=1)
    d = varargin{1};
end
if (length(varargin)>=2)
    xnames = varargin{2};
end
if(length(varargin)>=3)
    ignoretarget=varargin{3};
else
    ignoretarget=0;
end

[n(1),n(2),n(3)] = get_dim(d);

%% make sure the nb of target will not exceed 1

if (n(3)>1),
    error('The number of output dimensions in WEKA must not exceed 1.')   
end;

%% compute the attribute names for the input
if (isempty(xnames)),
    for i=1:n(2),
        xnames{i} = ['inp ' num2str(i)];
    end;
end;

%% =======================================================
%% creates the list of attributes, adds one for the target
attributes = weka.core.FastVector(n(2)+n(3));
for i = 1:n(2),
    %% add numeric feature
    attributes.addElement(weka.core.Attribute(xnames{i}));
end

%% =======================================================
if (n(3)>0 & ignoretarget==0)
    attributes.addElement(weka.core.Attribute('target'));
end

%% creates the instances object
dw = weka.core.Instances(get_name(d),attributes,n(1));

%% set the class (target) index
%% note that indexation starts at 0 (hence the n(1) and not n(1)+1
% BUG : Was  dw.setClassIndex(n(1)); must be n(2)
if (n(3)>0)
    dw.setClassIndex(n(2));
end

%% adds to the instances the list of input-output pairs contained in d
%% The last index is the target. 
x=get_x(d);

if (n(3)>0 & ignoretarget==0)
    y=get_y(d);
    for i = 1:n(1)
        dw.add(weka.core.Instance(1.0, [x(i,:) y(i,1)]));
    end
else
    for i = 1:n(1)
        dw.add(weka.core.Instance(1.0, x(i,:)));
    end
end

