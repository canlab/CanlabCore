function dw = wekaCategoricalData(varargin),
% DataWeka = WEKACATEGORICALDATA(d,[xnames,ynames]);
% 
% Transform a spider data object into a weka instances object (with
% categorical/nominal target)
%
% By default everything is numerical except the target which can
% be nominal if the number of values is limited.
%
% Parameters
%
%   d           -- spider data object
%   xnames      -- optional list of string for the names of the input
%               features
%   ynames      -- optional list of string for the names of the target
%               values. The classes in spider are supposed to be indexed
%               from 1 to Q (nb of classes). This string array refers to
%               the same ordering
%
% Output
% 
%   inst        -- weka instances object where the target has been set.
%               Note that the target is nominal here.

xnames=[];
ynames=[];
if (length(varargin)>=1)
    d = varargin{1};
end
if (length(varargin)>=2)
    xnames = varargin{2};
end
if (length(varargin)>=3)
    ynames = varargin{3};
end

[n(1),n(2),n(3)] = get_dim(d);
if (n(3)==1),
    n(3)=n(3)+1;
end;
%% compute the attribute names for the input
if (isempty(xnames)),
    for i=1:n(2),
        xnames{i} = ['inp ' num2str(i)];
    end;
end;

%% compute the list of values for the output
if (isempty(ynames)),   
   for i=1:n(3),
       ynames{i} = ['out ' num2str(i)];
   end
end

%% creates the FastVector for the class attribute
classValues = weka.core.FastVector(length(ynames));
for i = 1:length(ynames),
   classValues.addElement(ynames{i}); 
end

%% transform the output y into an index of the position of the y value
%% into ynames
y = get_y(d);
if (size(y,2)==1),
    %% two class problem
    ytmp = y;
    %% -1 class corresponds to the first element of ynames (indexation in
    %% weka starts at zero)
    ytmp(y==-1)=0;
    %% 1 class is the second element of ynames
    ytmp(y==1)=1;
else
    ytmp = ones(size(y,1),1);
    for i = 1:size(y,2),
        %% indexation in weka starts at zero
       ytmp(y(:,i)==1)=i-1; 
    end
end

%% creates the list of attributes, adds one for the target
attributes = weka.core.FastVector(n(2)+1); 
for i = 1:n(2),
    %% add numeric feature
    attributes.addElement(weka.core.Attribute(xnames{i}));
end
attributes.addElement(weka.core.Attribute('target',classValues));
    
%% creates the instances object
dw = weka.core.Instances(get_name(d),attributes,n(1));

%% set the class (target) index
%% note that indexation starts at 0 (hence the n(2) and not n(2)+1
dw.setClassIndex(n(2));

%% adds to the instances the list of input-output pairs contained in d
x=get_x(d);
for i = 1:n(1),
   dw.add(weka.core.Instance(1.0, [x(i,:) ytmp(i)])); 
end