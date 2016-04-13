function publish_obj_oriented_help()

outputbase = '/Users/tor/Dropbox/psyc7215_class_files/Part2_Machine_Learning/Toolboxes/help'; 
if ~exist(outputbase, 'dir'), mkdir(outputbase); end


for objnames = {'fmri_data' 'statistic_image' 'image_vector' 'fmridisplay'}

objname = objnames{1};

outputdir = fullfile(outputbase, objname);
if ~exist(outputdir, 'dir'), mkdir(outputdir); end

% Dynamically generate script for each object
% This allows us to use markup to create headers, etc.
%
% output in main help directory

%scriptname = 'print_obj_oriented_help.m';
scriptname = generate_script(outputdir, objname);
addpath(outputdir)

p = struct('useNewFigure', false, 'maxHeight', 1500, 'maxWidth', 1200, ...
    'outputDir', outputdir, 'showCode', false);

publish(scriptname, p)

end

end % main function

function scriptname = generate_script(outputdir, objname)
% generates script for publish

scriptname = fullfile(outputdir, ['index.m']);  
[fid, message] = fopen(scriptname, 'w');

if fid < 0, disp(message), end

z = '_________________________________________';

maxlines = 10; % lines of help to print

% Start dynamically generating script
% Get object and method names in script
    str = ['objname = ''' objname ''';'];
    fprintf(fid, '%s\n', str);
   
    str = ['m = methods(objname);'];
    fprintf(fid, '%s\n', str);

        str = ['maxlen = 200; % characters'];
    fprintf(fid, '%s\n', str);
 
    
% Do object constructor method first

fprintf(fid, '%%%% %s Object Class\n', objname); 

h = help(objname);
%len = min(maxlen, length(h));

% Separate comment lines for Description markup

helptext = get_first_help_lines(objname, maxlines);

for i = 1:length(helptext)
    fprintf(fid, '%% %s\n', helptext{i});
end

fprintf(fid, '%% \n%% METHODS\n%% ___________________________________\n\n');

m = methods(objname);

for i = 1:length(m)
    %%% Method
    % New method begins here
    
    fprintf(fid, '%%%% %s\n', m{i}); 

    
    % Get help (in script)
    

    str = ['mname = [objname ''.' m{i} '''];'];
    fprintf(fid, '%s\n', str);
    
    str = ['h = help(mname);'];
    fprintf(fid, '%s\n', str);
    
    str = ['len = min(maxlen, length(h));'];
    fprintf(fid, '%s\n', str);
        
    % get and print method name in bold
    mname = [objname '.' m{i}];
    
    fprintf(fid, '\n%% *%s*\n', mname);
    
    % print help
    fprintf(fid, 'fprintf(''%%s\\n'', h(1:len));\n');
        
end

fclose(fid);

% for development, run this to print script:
% eval(['!more ' scriptname])



end % function

