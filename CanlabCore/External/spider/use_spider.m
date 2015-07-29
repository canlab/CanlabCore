function pOut = use_spider(add)
% USE_SPIDER   get ready to use the Spider toolbox
% 
% To install the Spider toolbox, add the Spider root directory permanently
% to your MATLAB search path using, for example, PATHTOOL.
% 
% To begin using the Spider, call USE_SPIDER. If you use the Spider all
% the time, it would be a good idea to call USE_SPIDER in your startup.m
% file.
% 


%% ................DEPRECATED.........
%% gb
% USE_SPIDER
% USE_SPIDER(1)
%     Adds the Spider subdirectories to the MATLAB search path, and also
%     sets up the default global options.
% 
% USE_SPIDER(0)
%     Removes the Spider subdirectories from the MATLAB search path.
% 
% P = USE_SPIDER
%     Does not add or remove paths, but returns a cell array of path
%     strings that would be added. They can then be added or removed
%     manually with ADDPATH(P{:}) or RMPATH(P{:}). Global options are not
%     set up.

if nargin < 1, add = 1; end

if nargin<1
disp(' ');
disp('SPIDER : a machine learning toolbox for Matlab(R).');
disp(' ');
disp('This program is free software; you can redistribute it and/or');
disp('modify it under the terms of the GNU General Public License');
disp('as published by the Free Software Foundation; either version 2');
disp('of the License, or (at your option) any later version.');
disp(' ');
disp('This program is distributed in the hope that it will be useful,');
disp('but WITHOUT ANY WARRANTY; without even the implied warranty of');
disp('MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the');
disp('GNU General Public License for more details:');
disp('http://www.gnu.org/copyleft/gpl.html');
disp(' ');
disp(' ');
end


if datenum(version('-date')) >= datenum('May 6 2004')
 % look at this lovely mutual incompatibility between R14+ and previous releases
 % well done MathWorks, you've made our lives needlessly difficult AGAIN!
 d = dbstack('-completenames');
 rootdir = fileparts(d(1).file);
else
 d = dbstack; 
 rootdir = fileparts(d(1).name);
end   

subdirs=spider_subdirs;
% if nargin 
addpath(subdirs{:})
addpath(rootdir);

s=spider_path(1);

if(length(strfind(javaclasspath,'weka'))>0)
    disp('WEKA support enabled!');
    disp(' ');
end



% show all default values for algorithms e.g C=Inf in svm, otherwise suppress  
global display_tree_show_defaults;     
if isempty(display_tree_show_defaults), display_tree_show_defaults=0; end

% doesn't draw whole tree, only up to depth given  
global display_tree_depth;              
if isempty(display_tree_depth), display_tree_depth=100; end

% displays indexes e.g [1 3 4] if set to 1,  
% displays last element of array e.g 4 of [1 3 4] if set to 2,  
% or else displays a single integer index system for whole tree   
%              (breadth-wise, leftmost is smallest value)  
global display_tree_array_indexing;    
if isempty(display_tree_array_indexing), display_tree_array_indexing=2; end

% place index display at front or end of each line of output  
global display_tree_index_at_front;    
if isempty(display_tree_index_at_front), display_tree_index_at_front=1; end

% draw brackets in display for object sets  
global display_tree_show_brackets;     
if isempty(display_tree_show_brackets), display_tree_show_brackets=1; end

% number of characters of tab spacing for nested objects in tree   
global display_tree_tab_step;          
if isempty(display_tree_tab_step), display_tree_tab_step=3; end

% whether to allow recursive search for variables in statements  
% such as:  b=algs({svm svm}); b.C=2  
% (what you lose is it is much slower...)  

global recursive_subsasgn_off;  
if isempty(recursive_subsasgn_off), recursive_subsasgn_off=0; end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = updir(d)
d = fileparts(d);
if d(end)==filesep, d(end) = []; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
