function queryoutput = canlab_clone_github_repository(varargin)
% Clone a named CANlab repository from Github, and add with subfolders to
% the Matlab path.  Default repository is canlabCore.  Optional inputs can
% save the new path.
%
% See https://github.com/canlab for repository names. You can also download
% or clone repositories from this site.
%
% This function will clone repositories into the current working directory,
% so run this function from your local repository parent directory.
%
% :Usage:
% ::
%
%     repo_info_structure = canlab_clone_github_repository([optional inputs])
% ..
%     Author and copyright information:
%
%     Copyright (C) July 2018, Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Optional Inputs:
%   **'reponame', 'repo', 'repository':**
%        Followed by name of repository to download.
%
%   **'groupurl':**
%        Followed by URL for user/group in Github. You could use this to
%        download repositories from other groups outside CANlab.
%
%   **'savepath':**
%        Save the Matlab path with the new repository folders added.
%        This is not done by default to avoid unintentional changes to your
%        permanent path.
%
%   **'dryrun':**
%        Skip actual clone operation - no files copied.
%
%   **'noverbose':**
%        Do not print full output from clone operation.
%
% :Outputs:
%
%   **'repo_info_structure':**
%        Query output structure from Github API query, decoded from JSON
%
%   Also: Downloads and adds toolbox folders to your Matlab path.
%
% :Examples:
% ::
% % clone canlabCore and do not save path:
%    canlab_clone_github_repository
%
% % dry run - retrieve info for CANlab_help_examples: 
%    canlab_clone_github_repository('repo', 'CANlab_help_examples', 'dryrun'); 
%
% % save Matlab path:
%   canlab_clone_github_repository('repo', 'CANlab_help_examples', 'savepath'); 
%
% :See also:
%   - canlab_toolbox_setup
%

% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
%    July 2018: Created by Tor Wager
% ..

% Notes from web:
% clone all twitter repositories using API
%curl -s https://api.github.com/orgs/twitter/repos?per_page=200 | ruby -rubygems -e 'require "json"; JSON.load(STDIN.read).each { |repo| %x[git clone #{repo["ssh_url"]} ]}'

% -------------------------------------------------------------------------
% DEFAULTS AND INPUTS
% -------------------------------------------------------------------------


reponame = 'canlabCore';
groupurl = 'https://api.github.com/repos/canlab/';
dosavepath = false;
oktogo = true;  % or false for dry run
doverbose = true;

% optional inputs with default values
for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            
            case {'reponame', 'repo', 'repository'}, reponame = varargin{i+1}; varargin{i+1} = [];
            case 'groupurl', groupurl = varargin{i+1}; varargin{i+1} = [];
            case 'savepath', dosavepath = true;
            case 'dryrun', oktogo = false;
            case 'noverbose', doverbose = false;
                
            otherwise, warning(['Unknown input string option:' varargin{i}]);
        end
    end
end

% -------------------------------------------------------------------------
% CLONE REPO
% -------------------------------------------------------------------------


repofullname = [groupurl reponame];

[isok, queryoutput] = system(['curl ' repofullname]);
if isok == 0
    fprintf('curl retrieved %s from Github\n', repofullname);
else
    disp(queryoutput);
    error('Could not retrieve repository info from Github - missing or private?');
end

if ~exist('jsondecode', 'builtin')
    error('Cannot find Matlab jsondecode. Old Matlab version? Skipping repository import')
end

try
    queryoutput = jsondecode(queryoutput);
catch
    disp('jsondecode error - skipping repository import');
    return
end

if isfield(queryoutput, 'clone_url') % only for public repos, without authenticating
    
    cloneurl = queryoutput.clone_url;
    clone_into_dir = pwd;
    
else
    disp('Cannot clone repo - maybe private?');
    return
end


%oktogo = input(sprintf('Clone %s into local dir %s?\nType 1 for yes, 0 for no: ', reponame, clone_into_dir));

if oktogo
    
    fprintf('Cloning %s into local dir %s\n', reponame, clone_into_dir);
    
    [isok2, queryoutput2] = system(['git clone ' cloneurl]);
    
    if doverbose
        disp(queryoutput2);
    end
    
    if isok2 == 0
        
        fprintf('git successfully cloned %s from Github\n', cloneurl);
        
        % -------------------------------------------------------------------------
        % Add to path
        % -------------------------------------------------------------------------

        add_to_path_dir = fullfile(clone_into_dir, reponame);
        g = genpath(add_to_path_dir);
        addpath(g);
        
        if dosavepath
            savepath
            fprintf('Added %s\nto Matlab path with subfolders and saved path.\n', add_to_path_dir);
            
        else
            fprintf('Added %s\nto Matlab path with subfolders. Type savepath in Matlab to save this.\n', add_to_path_dir);
        end
        
    else
        
        error('Could not retrieve repository info.');
    end
    
else
    
    disp('Dry run: Skipping clone operation.');
    return
    
end % oktogo

end % function
