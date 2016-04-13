function new_struct = struct_strrep(old_struct, old_string, new_string, depth)
% Traverses depth-first through an entire structure, replacing 
% old_string with new_string everywhere it goes
%
% :Usage:
% ::
%
%     new_struct = struct_strrep(old_struct, old_string, new_string)

    MAX_STRUCT_DEPTH = 200;
    new_struct = old_struct;

    if(~exist('depth', 'var') || isempty(depth))
        depth = 1;
    elseif(depth > MAX_STRUCT_DEPTH)
        error('Exceeding max structure depth of %d. If this is legitimate, please edit this function.\n', depth);
    end

    fnames = fieldnames(new_struct);
    for i = 1:length(fnames)
        % Below code does not work as Matlab R2007B has some massive bug in the class() function.
        % It works when executed manualy, but not as part of a func... wtf?
        %         field_class = class(new_struct.(fnames{i}));
        %         switch(field_class)
        %             case 'struct'
        %                 new_struct.(fnames{i}) = struct_strrep(new_struct.(fnames{i}), old_string, new_string);
        %             case 'cell'
        %                 if(iscellstr(new_struct.(fnames{i})))
        %                     new_struct.(fnames{i}) = strrep(new_struct.(fnames{i}), old_string, new_string);
        %                 end
        %             case 'char'
        %                 new_struct.(fnames{i}) = char(strrep(cellstr(new_struct.(fnames{i})), old_string, new_string));
        %         end

        if(isstruct(new_struct.(fnames{i})))
            for j=1:length(new_struct.(fnames{i}))
                new_struct.(fnames{i})(j) = struct_strrep(new_struct.(fnames{i})(j), old_string, new_string, depth + 1);
            end
        elseif(iscellstr(new_struct.(fnames{i})))
            new_struct.(fnames{i}) = strrep(new_struct.(fnames{i}), old_string, new_string);
        elseif(ischar(new_struct.(fnames{i})))
            new_struct.(fnames{i}) = char(strrep(cellstr(new_struct.(fnames{i})), old_string, new_string));
        elseif(iscell(new_struct.(fnames{i})))
            for j=1:length(new_struct.(fnames{i}))
                new_struct.(fnames{i}){j} = struct_strrep(new_struct.(fnames{i}){j}, old_string, new_string, depth + 1);
            end
        end
    end
end
