function struct2yaml(yamlfilename, DB, yamlfilemethod, dbmethod)
% :Usage:
% ::
%
%     struct2yaml(yamlfilename, DB, yamlfilemethod, dbmethod)
%
% :Inputs:
%
%   **yamlfilemethod:**
%        'new' or 'add' (append)
%
%   **dbmethod:**
%        how the canlab database will handle the record.
%        'add', 'replace', or 'keep_existing'
%
% translate structure into YAML format text file
% this will be interpretable by the canlab database
%
% :Examples:
% ::
%
%    yamlfilename = 'YAML_tmp.yaml';
%
%    DB.study = 'NSF';     % string; study code letters
%    DB.subject = '001';   % string; subject ID number
%    DB.occasion = '21';   % string; occasion ID; unique to subj*session
%    DB.unique_id = [DB.study '_' DB.subject '_' DB.occasion];
%    DB.mean_spikes_per_image = mean(cat(2, spikesperimg{:}));
%
%    struct2yaml(yamlfilename, DB, 'add', 'replace');
%
% ..
%    See canlab_database WIKI page for more details (internal access only.)
%
%    you could upload your file to the canlab repository like this:
%    rsync -v /Users/tor/Documents/Tor_Documents/Coursework_and_Teaching/PSYC_7215_Fall_2010/Sample_data_NSF_study/SubjectData/denoised_canlab/nsf2/YAML_tmp.yaml canlab.colorado.edu:/Volumes/RAID1/labdata/qc_yaml_repository/
% ..

switch dbmethod           % database behavior - instructions to canlab database
    case 'add'              % add to existing .yaml file.
    case 'replace'        % overwrite entire existing record
    case 'keep_existing'    % keep existing record in database if it exists already
        
    otherwise
        error('dbmethod must be add, replace, or keep_existing');
end
DB.method = dbmethod;

switch yamlfilemethod
    case 'new'
        fopenstring = 'w+';
        fid = fopen(yamlfilename,fopenstring);
        fprintf(fid, 'record:\n');
        
    case 'add'
        fopenstring = 'a+';
        fid = fopen(yamlfilename,fopenstring);
        
    otherwise
        error('yamlfilemethod must be add, replace, or keep_existing');
end

N = fieldnames(DB);

for i = 1:length(N)
    fprintf(fid, '  %s: ', N{i});
    
    
    switch class(DB.(N{i}))
        case {'double' 'single'}
            fprintf(fid, '%3.6f\n', DB.(N{i}));
            
        case {'uint8', 'uint16', 'logical'}
            fprintf(fid, '%3.0f\n', DB.(N{i}));
    
        case {'char'}
            fprintf(fid, '%s\n', DB.(N{i}));
    
        otherwise
            warning('struct2yaml:unknownClass', 'Unknown class for field %s. Skipping.', N{i});
            
    end
    
end

status = fclose(fid);

switch status
    case 0
        fprintf('yaml file %s created successfully.\n', yamlfilename);
    case -1
        fprintf('yaml file %s DID NOT write successfully!\n', yamlfilename);
    otherwise
        disp('Unknown fclose status. ???');
end



end
