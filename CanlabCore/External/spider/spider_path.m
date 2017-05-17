function s = spider_path(set)
% This function returns the path where spider is installed. 
% If an arugment is provided than it adds all spider subdirectory paths
% to the Matlab search path.

subdirs=spider_subdirs;

% if nargin==0 %% can just return s and do not set paths
if nargin>0
    for i=1:length(subdirs)
        t = cellstr(subdirs{i});
        path(path,t{:});
    end
s=subdirs{1};


end

s=subdirs{1};
k = findstr(s,'spider');
s=s(1:k+5);


if nargin > 0
    path(path,s);
else
    subdirs=spider_subdirs;
% if nargin 
    addpath(subdirs{:})
    addpath(s);

end



% show all default values for algorithms e.g C=Inf in svm, otherwise suppress
global display_tree_show_defaults;
display_tree_show_defaults=0;

% doesn't draw whole tree, only up to depth given
global display_tree_depth;
display_tree_depth=100;

% displays indexes e.g [1 3 4] if set to 1,
% displays last element of array e.g 4 of [1 3 4] if set to 2,
% or else displays a single integer index system for whole tree
%              (breadth-wise, leftmost is smallest value)
global display_tree_array_indexing;
display_tree_array_indexing=2;

% place index display at front or end of each line of output
global display_tree_index_at_front;
display_tree_index_at_front=1;

% draw brackets in display for object sets
global display_tree_show_brackets;
display_tree_show_brackets=1;

% number of characters of tab spacing for nested objects in tree
global display_tree_tab_step;
display_tree_tab_step=3;

% whether to allow recursive search for variables in statements
% such as:  b=algs({svm svm}); b.C=2
% (what you lose is it is much slower...)

global recursive_subsasgn_off;
recursive_subsasgn_off=0;
% else


if isunix
     s = [ s '/'];
 else
     s = [ s '\'];
 end


s(s=='\')='/';

%% Take into account weka.jar
wekajar=[matlabroot filesep 'java' filesep 'jar' filesep 'weka.jar'];
if( length(dir(wekajar))>0)
    javaaddpath(wekajar);
end

