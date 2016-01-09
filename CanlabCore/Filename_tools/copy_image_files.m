function output_names = copy_image_files(image_list, to_dir, varargin)
% Copies a set of image files from one directory to another, creating the
% directory if needed.
%
% :Usage:
% ::
%
%     function output_names = copy_image_files(image_list, to_dir, [method='copy' or 'move'], varargin)
%
% MAC OSX only!!
%
% For .img files, will look for a paired .hdr and copy it automatically.
% Works for non-image files (e.g., .txt) as well.
%
% Other variable args:
%
% Will take flags:
%
% 'append image number' -> will append a number to each file corresponding
% to the order listed, e.g., 1 for the first, 2 for the 2nd, etc.
%
% 'append string' -> followed by string to append
%
% with both flags and 'run' for append string, a matrix of
% [run1/vols.img; run2/vols.img] would become [vols_run0001.img vols_run0002.img]

method = 'copy';
donumberappend = 0;
appstr = [];

if ~isempty(varargin), method = varargin{1}; end

for i = 2:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % reserved keywords
            case 'append image number', donumberappend = 1;
            case 'append string', appstr = varargin{i + 1}; varargin{i + 1} = [];
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

switch method
    case {'copy', 'test'}
        strbase = ['!cp ']; % copy, for mac
        
    case 'move'
        strbase = ['!mv ']; % copy, for mac
        
    otherwise error('Enter ''copy'', ''move'', or ''test'' for method')
end

if strcmp(method, 'test')
    disp('TEST ONLY -- NOTHING WILL BE COPIED');
end






% get/make output directory name
if ~exist(to_dir, 'dir')
    fprintf('Creating directory: %s\n', to_dir);
    mkdir(to_dir)
end

output_names = [];

for i = 1:size(image_list, 1)  % for each
    
    numappstr = '';
    if donumberappend, numappstr = sprintf('%04d', i); end
    
    img_name = deblank(image_list(i, :));
    
    % fix SPM-specific: remove ,# at end if added by SPM
    if img_name(end-1) == ','
        img_name = img_name(1:end-2);
    end
    
    if ~exist(img_name, 'file')
        warning('copy_image_files:image does not exist', 'Check locations?');
    end
    
    [dd, ff, ee] = fileparts(img_name);
    
    % add escape characters before spaces
    dd = add_escape_chars(dd);
    
    
    inname1 = [ff ee];
    outname1 = [ff appstr numappstr ee];
    
    
    str1 = [strbase fullfile(dd, inname1) ' ' fullfile(to_dir, outname1)];
    disp(str1)
    
    if ~strcmp(method, 'test')
        eval(str1)
    end
    
    if strcmp(ee, '.img')  % if Analyze...
        inname2 = [ff '.hdr'];
        inname3 = [ff '.mat'];
        outname2 = [ff  appstr numappstr '.hdr'];  % look for header as well
        outname3 = [ff  appstr numappstr '.mat'];  % look for mat file as well (SPM)
        
        if ~exist(fullfile(dd, inname2), 'file')
            warning('copy_image_files:missing .hdr', '.hdr file missing!!!');
        end
        
        str2 = [strbase fullfile(dd, inname2) ' ' fullfile(to_dir, outname2)];
        disp(str2)
        
        str3 = [strbase fullfile(dd, inname3) ' ' fullfile(to_dir, outname3)];
        if exist(fullfile(dd, inname3), 'file'), disp(str3), end
        
        %should use this instead
        %[ok, msg, msgid] = copyfile(subject_from_dir_final, subject_to_dir);
        
        
        if ~strcmp(method, 'test')
            eval(str2)
            if exist(fullfile(dd, inname3), 'file'), eval(str3), end
        end
        
    end
    
    output_names = strvcat(output_names, fullfile(to_dir, outname1));
    
end


disp('Done.')
disp(' ');

end


function newff = add_escape_chars(ff)

wh = find(deblank(ff) == ' ');
wh = [1 wh length(deblank(ff)) + 1];
echar = '\';

for i = 1:length(wh) - 1
    
    newff{1, i} = ff(wh(i):wh(i+1) - 1);
    
    if length(wh) > 2 && i < length(wh) - 2
        newff{1, i} = [newff{1, i} echar];
    end
    
end

newff = cat(2, newff{:});

end



