function new_var = strrep_recurse(old_var, old_string, new_string, depth)
% Recursively traverses depth-first through an entire variable, replacing
% old_string with new_string everywhere it goes
%
% :Usage:
% ::
%
%     new_var = strrep_recurse(old_var, old_string, new_string)

    MAX_DEPTH = 200;
    new_var = old_var;

    if(~exist('depth', 'var') || isempty(depth))
        depth = 1;
    elseif(depth > MAX_DEPTH)
        error('Exceeding max structure depth of %d. If this is legitimate, please edit this function.\n', depth);
    end

    % Below code does not work as Matlab R2007B has some massive bug in the class() function.
    % It works when executed manually, but not as part of a func... wtf?
    %         field_class = class(new_var.(fnames{i}));
    %         switch(field_class)
    %             case 'struct'
    %                 new_var.(fnames{i}) = struct_strrep(new_var.(fnames{i}), old_string, new_string);
    %             case 'cell'
    %                 if(iscellstr(new_var.(fnames{i})))
    %                     new_var.(fnames{i}) = strrep(new_var.(fnames{i}), old_string, new_string);
    %                 end
    %             case 'char'
    %                 new_var.(fnames{i}) = char(strrep(cellstr(new_var.(fnames{i})), old_string, new_string));
    %         end

    % Not using switch due to some weird bug in class() function in 2007B
    if(isstruct(new_var))
        for i=1:length(new_var)
            fnames = fieldnames(new_var(i));
            for j=1:length(fnames)
                new_var(i).(fnames{j}) = strrep_recurse(new_var(i).(fnames{j}), old_string, new_string, depth + 1);
            end
        end
    elseif(iscellstr(new_var))
        new_var = strrep(new_var, old_string, new_string);
    elseif(ischar(new_var))
        new_var = char(strrep(cellstr(new_var), old_string, new_string));
    elseif(iscell(new_var))
        for i=1:length(new_var)
            new_var{i} = strrep_recurse(new_var{i}, old_string, new_string, depth + 1);
        end
    end
end
