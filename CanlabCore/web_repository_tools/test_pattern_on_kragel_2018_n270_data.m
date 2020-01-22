function test_pattern_on_kragel_2018_n270_data(varargin)
%
% test a multivariate signature pattern stored in an fmri_data object on Kragel et al. 2018 datasets, N = 270
%
% -   requires kragel_2018_nat_neurosci_270_subjects_test_images.mat
%     if not found, will attempt to download from Neurovault using
%     retrieve_neurovault_collection(). See load_image_set('kragel18_alldata')
%
% -   Dataset is from Kragel et al. 2018, Nature Neuroscience. 
%     Uses load_image_set to load it:
%     test_images = load_image_set('kragel18_alldata', 'noverbose');
%
% -   If not found, load_image_set uses the utility retrieve_neurovault_collection 
%     to pull the images from Neurovault.org if you don't have them.
%     [files_on_disk, url_on_neurovault, mycollection, myimages] = retrieve_neurovault_collection(504);
%     image_obj = fmri_data(files_on_disk);
%     It also adds meta-data to label the images.
%
% :Usage:
% ::
%
%     [list outputs here] = test_pattern_on_kragel_2018_n270_data(obj)
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2019 Tor Wager
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
% :Inputs:
%
%   **obj:**
%        An fmri_data object storing a single image of multivariate pattern weights
%
% :Optional Inputs:
%   **'noverbose':**
%        Turn off verbose output
%
%   **'noplots':**
%       Turn off plots
%
%   **'threshold_type'
%       Followed by keyphrase, e.g., 'Optimal balanced error rate'
%       Default is 'Optimal overall accuracy'
%
%
% :Outputs:
%
%   **out1:**
%        description of out1
%
%   **out2:**
%        description of out2
%
% :Examples:
% ::
%
%    % give examples of code here
%    param1 = abc();
%    param2 = xyz();
%    [out1,out2] = func_call(param1, param2)
%
% :References:
%   CITATION(s) HERE
%
% :See also:
%   - list other functions related to this one, and alternatives*
%


% ..
%    Programmers' notes:
%    List dates and changes here, and author of changes
% ..

% ----------------------------------------------------------------------
% Defaults
% ----------------------------------------------------------------------
% Must define ALL inputs here. All these can be changed by passing in arguments.

obj = [];
num_cols = 18;      % this is for the test dataset
doverbose = true;
doplot = true;
fontsize = 14;
threshold_type = 'Optimal balanced error rate'; % 'Optimal overall accuracy'; %

colors = [repmat({[1 .2 0]}, 1, 6) repmat({[.4 .5 .2]}, 1, 6) repmat({[.1 0 .9]}, 1, 6)];

sim_string = 'cosine_similarity';  % or 'dotproduct' or 'correlation'

canlab_parse_inputs_subfcn();

% ----------------------------------------------------------------------
% Display helper functions: Called by later scripts
% ----------------------------------------------------------------------

dashes = '----------------------------------------------';
printstr = @(dashes) disp(dashes);
printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);

disp('Process: test_pattern_on_kragel_2018_n270_data');
printstr(dashes)

disp(' ');

% ----------------------------------------------------------------------
% Load and prep test data
% ----------------------------------------------------------------------
test_images = load_image_set('kragel18_alldata', 'noverbose');

test_images.Y = test_images.dat_descrip.Studynumber;

vector_data = apply_mask(test_images, obj, 'pattern_expression', sim_string);

for s = 1:num_cols
    pattern_response{s} = vector_data(test_images.Y==s);
end

% clear exp
labels = test_images.additional_info;

% ----------------------------------------------------------------------
% ROC and stats
% ----------------------------------------------------------------------

[ispain, iscog, isemo] = deal(false(270, 1));

ispain(1:6*15) = true;
iscog(6*15 + 1 : (12*15)) = true;
isemo(12*15+1:18*15) = true;

% Boostrap confidence intervals for effect size
fprintf('Bootstrapping CIs\n')

bci_pain = bootci(5000, @cohens_d_2sample, vector_data, ispain);

bci_cog = bootci(5000, @cohens_d_2sample, vector_data, iscog);

bci_emo = bootci(5000, @cohens_d_2sample, vector_data, isemo);

% ----------------------------------------------------------------------
% ROC calculation and classification output
% ----------------------------------------------------------------------

printhdr('Pattern classification results')
disp('Applied pattern to test data from 18 studies, n = 270 participants')
disp('From Kragel et al. 2018, Nature Neuroscience')
disp(' ')

disp('Classification of pain vs. other')

roc_pain = roc_plot(vector_data,ispain, 'noplot', 'threshold_type', threshold_type);

fprintf('Threshold for cosine_sim = %3.4f, Cohen''s d(pain vs no) = %3.2f, Bootstrapped CI is [%3.2f %3.2f]\n', ...
    roc_pain.class_threshold, cohens_d_2sample(vector_data, ispain), bci_pain(1), bci_pain(2));

disp(' ')

disp('Classification of cognitive control vs. other')

roc_cog = roc_plot(vector_data,iscog, 'noplot', 'threshold_type', threshold_type);

fprintf('Threshold for cosine_sim = %3.4f, Cohen''s d(pain vs no) = %3.2f, Bootstrapped CI is [%3.2f %3.2f]\n', ...
    roc_cog.class_threshold, cohens_d_2sample(vector_data, iscog), bci_cog(1), bci_cog(2));

disp('Classification of emotion vs. other')

roc_emo = roc_plot(vector_data,isemo, 'noplot', 'threshold_type', threshold_type);

fprintf('Threshold for cosine_sim = %3.4f, Cohen''s d(pain vs no) = %3.2f, Bootstrapped CI is [%3.2f %3.2f]\n', ...
    roc_emo.class_threshold, cohens_d_2sample(vector_data, isemo), bci_emo(1), bci_emo(2));
disp(' ')

% ----------------------------------------------------------------------
% Plot
% ----------------------------------------------------------------------
if doplot
    
    create_figure('plot', 1, 3);
    
    clear hh
    
    subplot(1,3,3)
    
    roc = roc_plot(vector_data,ispain, 'threshold_type', threshold_type, 'nooutput');
    set(roc.line_handle, 'Color', 'r');
    hh(1) = roc.line_handle(2);
    
    roc = roc_plot(vector_data,iscog, 'threshold_type', threshold_type, 'nooutput');
    set(roc.line_handle, 'Color', 'g');
    hh(2) = roc.line_handle(2);
    
    roc = roc_plot(vector_data,isemo, 'threshold_type', threshold_type, 'nooutput');
    set(roc.line_handle, 'Color', 'b');
    hh(3) = roc.line_handle(2);
    
    set(gca,'fontsize', fontsize)
    legend(hh, {'Pain vs Other' 'Cog vs. Other' 'Emo vs. Other'});
    
    printhdr('Stats for individual studies')
    
    subplot(1,3,1:2)
    
    set(gca, 'FontSize', fontsize);
    barplot_columns(pattern_response,'names',labels,'colors',colors,'nofig');
    ylabel 'NPS Response (Cosine Similarity)'
    hold on; 
    
    if roc_pain.accuracy_p < .05
        plot([0 19],[roc_pain.class_threshold roc_pain.class_threshold],'r--')
    end
    
    if roc_cog.accuracy_p < .05
        plot([0 19],[roc_cog.class_threshold roc_cog.class_threshold],'g--')
    end
    
    if roc_emo.accuracy_p < .05
        plot([0 19],[roc_emo.class_threshold roc_emo.class_threshold],'b--')
    end
    
else
    % no plot, but print stats
    
    printhdr('Stats for individual studies')
    barplot_columns(pattern_response,'names',labels, 'colors',colors, 'skipallplots');
    
    
end % plot


% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% Nested functions
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

    function canlab_parse_inputs_subfcn()
        % INPUT PARSER: This is a nested function.
        % It modifies variables in the calling workspace directly.
        % It validates inputs.
        % It uses inputParser, but also allows for special keywords that change
        % variables, which are not handled by inputParser.
        
        % ----------------------------------------------------------------------
        % Parse inputs
        % ----------------------------------------------------------------------
        p = inputParser;
        
        % Validation functions - customized for each type of input
        % ----------------------------------------------------------------------
        
        valfcn_scalar = @(x) validateattributes(x, {'numeric'}, {'nonempty', 'scalar'});
        valfcn_number = @(x) validateattributes(x, {'numeric'}, {'nonempty'}); % scalar or vector
        valfcn_cell = @(x) validateattributes(x, {'cell'}, {'nonempty'}, {'numel'}, num_cols);
        
        valfcn_image_vector = @(x) isa(x, 'image_vector');
        
        % Required inputs
        % ----------------------------------------------------------------------
        p.addRequired('obj', valfcn_image_vector);
        
        % Optional inputs - key-value pairs
        % ----------------------------------------------------------------------
        % Pattern: keyword, value, validation function handle
        
        p.addParameter('fontsize', fontsize);
        p.addParameter('colors', colors);
        p.addParameter('sim_string', sim_string);
        p.addParameter('threshold_type', threshold_type);
        
        % Optional inputs - Logical flags and keywords
        % ----------------------------------------------------------------------
        % Identify and remove these before parsing, because they do not follow name-value pair pattern
        wh = strcmp(varargin, 'noverbose');
        if any(wh), doverbose = false; varargin(wh) = []; end
        
        wh = strcmp(varargin, 'noplots');
        if any(wh), doplot = false; varargin(wh) = []; end
        
        wh = strcmp(varargin, 'noplot');
        if any(wh), doplot = false; varargin(wh) = []; end
        
        % Parse inputs and distribute out to variable names in workspace
        % ----------------------------------------------------------------------
        p.parse(varargin{:});
        
        IN = p.Results;
        fn = fieldnames(IN);
        
        for i = 1:length(fn)
            str = sprintf('%s = IN.(''%s'');', fn{i}, fn{i});
            eval(str)
        end
        
    end  % Parse subfunction


end % main function
