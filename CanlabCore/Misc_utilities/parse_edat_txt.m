function [edat_struct, edat_cells] = parse_edat_txt(fname)
% Reads EPrime .txt output (equivalent to EDAT) directly into matlab cell arrays/structures
%
% Inputs:
%
%   **fname:**
%        name of file to parse
%
% :Outputs:
%
%   **edat_struct:**
%        Structure containing the following fields:
%
%        header: 1-element structure whose fields are header items from EDAT
%
%        run: 1-element structure whose fields are run-specific items from EDAT
%
%        trials: n-element structure whose fields are trial-specific items from EDAT
%                where n is the number of trials
%
%   **edat_cells:**
%        Structure containing the following fields:
%
%        header_cols: 1-row cell array whose fields are header column names from EDAT
%
%        run_cols: 1-row cell array whose cells are run-specific column names from EDAT
%
%        trials_cols: 1-row cell array whose cells are trial-specific column names from EDAT
% 
%        header: 1-row cell array whose cells are header items from EDAT
%
%        run: 1-row cell array whose cells are run-specific items from EDAT
%
%        trials: n-row cell array whose cells are trial-specific items from EDAT
%                where n is the number of trials
%
% ..
%    EPrime helper function
%    2/22/10 Joe Wielgosz
% ..



% Two passes are necessary:
%  1) find all fields used in the file
%  2) actually read the data
for pass = 1:2
    
    if (isunix)
        warning off
        f = fopen(fname, 'r', 'l', 'UTF-16'); %this encoding worked for opening a file in Unix that was created by Windows
        warning on
    else
        f = fopen(fname);
    end
            
    line_idx = 1;
    row_idx = [0 0 0];
    start_row = false;
    end_row = false;

    % Modes during parse:
    % 1: header
    % 2: run
    % 3: trials
    for mode = 1:3
        if pass == 1
            template_struct{mode} = struct;
        else % pass == 2
            template_struct{mode} = orderfields(blank_struct(template_struct{mode}, ''));

            switch mode
                case {1}
                    edat_cells.header_cols = fieldnames(template_struct{mode})';
                case {2}
                    edat_cells.run_cols = fieldnames(template_struct{mode})';
                case {3}
                    edat_cells.trials_cols = fieldnames(template_struct{mode})';
            end
        end
    end

    while feof(f) == 0
        
        tline = fgetl(f);
        %if (skip), tline = tline(3:end); end
        
        %Check for start/end of block in file
        switch tline
            case {'*** Header Start ***'}
                mode = 1;
                start_row = true;
            case {'*** LogFrame Start ***'}
                mode = 2;
                start_row = true;
            case {'	*** LogFrame Start ***'}
                mode = 3;
                start_row = true;
            case {'*** Header End ***', '*** LogFrame End ***', '	*** LogFrame End ***'}
                end_row = true;
        end


        if start_row % Start of block in file

            row_idx(mode) = row_idx(mode)+1;
            col_idx = 1;
            start_row = false;
            row_struct = struct;

        elseif end_row % End of block in file

            if pass == 1
                
                % Include any field which has appeared before
                template_struct{mode} = combine_structs(template_struct{mode}, row_struct);

            else % pass == 2
                
                full_row_struct = orderfields(combine_structs(template_struct{mode}, row_struct));

                switch mode
                    case {1}
                        edat_struct.header(row_idx(mode)) = full_row_struct;
                        edat_cells.header(row_idx(mode),:) = struct2cell(full_row_struct);
                    case {2}
                        edat_struct.run(row_idx(mode)) = full_row_struct;
                        edat_cells.run(row_idx(mode),:) = struct2cell(full_row_struct);
                    case {3}
                        edat_struct.trials(row_idx(mode)) = full_row_struct;
                        edat_cells.trials(row_idx(mode),:) = struct2cell(full_row_struct);
                end
            end
            
            end_row = false;

        else % Middle of block, or junk

            pair = regexp(tline, '^\t*(?<col>[^:\t]+): (?<val>.+)$', 'names');
            
            if ~isempty(pair) % Ignore junk
                
                clean_col_name = regexprep(pair.col, '\.', '_');
                row_struct.(clean_col_name) = pair.val;
                col_idx = col_idx+1;
                
            end
        end
        
        line_idx = line_idx+1;

    end % while
    
    fclose(f);
    
end %for 

end



