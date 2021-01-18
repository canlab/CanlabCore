%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

function run_examples(whichOnes, verbose)
N_EXAMPLES = 21;

if nargin<1
    verbose='none';
    whichOnes=1:N_EXAMPLES;
else
    if nargin<2
        if ischar(whichOnes) ...
                && (strcmp(whichOnes, 'none') ...
                    || strcmp(whichOnes, 'text') ...
                    || strcmp(whichOnes, 'graphic'))
                verbose=whichOnes;
                whichOnes=2:N_EXAMPLES;
        else
            verbose='graphic';
        end
    end
    if ischar(whichOnes)
        whichOnes=str2double(whichOnes);
    end
    if ~isnumeric(whichOnes) || any(isnan(whichOnes)) || any(whichOnes<0) || any(whichOnes>N_EXAMPLES)
        error(['run_examples argument must be nothing or numbers from 1 to '...
            num2str(N_EXAMPLES) '!']);
    end
end
if ~strcmpi(verbose, 'none')
    if any(whichOnes < 19)
        UmapExamples.DisplayPub(3);
    end
    if any(whichOnes == 19)
        UmapExamples.DisplayPub(1);
    end
end
if all(whichOnes==0)
    whichOnes = 1:N_EXAMPLES;
end
    
if ismember(1, whichOnes)
    disp('run_umap Example 1 starting...');
    run_umap;
    disp('run_umap Example 1 completed with no MATLAB exceptions!');
end
if ismember(2, whichOnes)
    disp('run_umap Example 2 starting...');
    run_umap('sample30k.csv', 'save_template_file', 'utBalbc2D.mat', 'verbose', verbose);
    disp('run_umap Example 2 completed with no MATLAB exceptions!');
end
if ismember(3, whichOnes)
    disp('run_umap Example 3 starting...');
    run_umap('sample130k.csv', 'template_file', 'utBalbc2D.mat', 'verbose', verbose);
    disp('run_umap Example 3 completed with no MATLAB exceptions!');
end
if ismember(4, whichOnes)
    disp('run_umap Example 4 starting...');
    run_umap('sampleBalbcLabeled55k.csv', 'label_column', 11, 'label_file', 'balbcLabels.properties', 'save_template_file', 'ustBalbc2D.mat', 'verbose', verbose);
    disp('run_umap Example 4 completed with no MATLAB exceptions!');
end
if ismember(5, whichOnes)
    disp('run_umap Example 5 starting...');
    run_umap('sampleRag148k.csv', 'template_file', 'ustBalbc2D.mat', 'match_supervisors', 1, 'verbose', verbose);    
    disp('run_umap Example 5 completed with no MATLAB exceptions!');
end
if whichOnes==5.1
    disp('run_umap Example 5 starting with small example...');
    run_umap('sample10k.csv', 'template_file', 'ustBalbc2D.mat', 'verbose', verbose);
    disp('run_umap Example 5 completed with no MATLAB exceptions!');
end
if whichOnes==5.11
    disp('run_umap Example 5 starting with small example and ''verbose''==''none''...');
    run_umap('sample10k.csv', 'template_file', 'ustBalbc2D.mat', 'verbose', 'none');
    disp('run_umap Example 5 completed with no MATLAB exceptions!');
end
if ismember(6, whichOnes)
    disp('run_umap Example 6 starting...');
    [~,~, clusterIds]=run_umap('sample30k.csv', 'cluster_output', verbose, 'verbose', verbose);
    disp('run_umap Example 6 completed with no MATLAB exceptions!');
end
if ismember(7, whichOnes)
    disp('run_umap Example 7 starting...');
    [~, ~, clusterIds]=run_umap('sample30k.csv', 'n_components', 3, 'save_template_file', 'utBalbc3D.mat', 'verbose', verbose);
    disp('run_umap Example 7 completed with no MATLAB exceptions!');
end
if ismember(8, whichOnes)
    disp('run_umap Example 8 starting...');
    run_umap('sample130k.csv', 'template_file', 'utBalbc3D.mat', 'verbose', verbose);
    disp('run_umap Example 8 completed with no MATLAB exceptions!');
end
if ismember(9, whichOnes)
    disp('run_umap Example 9 starting...');
    run_umap('sampleRagLabeled60k.csv', 'label_column', 11, 'label_file', 'ragLabels.properties', 'save_template_file', 'ustRag2D.mat', 'verbose', verbose);
    disp('run_umap Example 9 completed with no MATLAB exceptions!');
    
end
if ismember(10, whichOnes)
    disp('run_umap Example 10 starting...');
    run_umap('sample30k.csv', 'template_file', 'ustRag2D.mat', 'verbose', verbose);
    disp('run_umap Example 10 completed with no MATLAB exceptions!');
    
end
if ismember(11, whichOnes)
    disp('run_umap Example 11 starting...');
    run_umap('sample30k.csv', 'template_file', 'ustRag2D.mat', 'method', 'Java', 'joined_transform', true, 'verbose', verbose);
    disp('run_umap Example 11 completed with no MATLAB exceptions!');
end
if ismember(12, whichOnes)
    disp('run_umap Example 12 starting...');
    run_umap('sampleRag148k.csv', 'template_file', 'ustBalbc2D.mat', 'qf_tree', true, 'qf_dissimilarity', true, 'see_training', true, 'verbose', verbose);
    disp('run_umap Example 12 completed with no MATLAB exceptions!');
end
if ismember(13, whichOnes)
    disp('run_umap Example 13 starting...');
    run_umap('sample30k.csv', 'verbose', verbose);
    run_umap('sample30k.csv', 'python', true, 'verbose', verbose);
    disp('run_umap Example 13 completed with no MATLAB exceptions!');
    
end
if ismember(14, whichOnes)
    disp('run_umap Example 14 starting (just MEX first)...');
    run_umap('sampleBalbcLabeled55k.csv', 'label_column', 11, 'label_file', 'balbcLabels.properties', 'save_template_file', 'ustBalbc2D.mat', 'verbose', verbose);
    disp('run_umap Example 14 (just MEX) completed with no MATLAB exceptions!');
end
if ismember(15, whichOnes)
    disp('run_umap Example 15 starting (just MEX first)...');
    run_umap('sampleRag55k.csv', 'template_file', 'ustBalbc2D.mat', 'verbose', verbose);
    disp('run_umap Example 15 (just MEX) completed with no MATLAB exceptions!');
end
if ismember(14, whichOnes)
    disp('run_umap Example 14 starting (with Python)...');
    run_umap('sampleBalbcLabeled55k.csv', 'label_column', 11, 'label_file', 'balbcLabels.properties', 'python', true, 'save_template_file', 'pyUstBalbc2D.mat', 'verbose', verbose);
    disp('run_umap Example 14 (with Python) completed with no MATLAB exceptions!');
end
if ismember(15, whichOnes)
    disp('run_umap Example 15 starting (with Python)...');
    run_umap('sampleRag55k.csv', 'template_file', 'pyUstBalbc2D.mat', 'verbose', verbose);
    disp('run_umap Example 15 (with Python) completed with no MATLAB exceptions!');
    
end
if ismember(16, whichOnes)
    disp('run_umap Example 16 starting...');
    run_umap('sampleBalbcLabeled55k.csv', 'label_column', 11, 'label_file', 'balbcLabels.properties', 'qf_tree', true, 'n_components', 3, 'save_template_file', 'ustBalbc3D.mat', 'verbose', verbose);
    [reduction, umap, clusterIdentifiers,extras]=run_umap('sample10k.csv', 'template_file', 'ustBalbc3D.mat', 'qf_tree', true, 'qf_dissimilarity', true, 'see_training', true, 'cluster_output', verbose, 'verbose', verbose);
    disp('run_umap Example 16 completed with no MATLAB exceptions!');
end
if ismember(17, whichOnes)
    disp('run_umap Example 17 starting...');
    run_umap('sampleBalbcLabeled12k.csv', 'template_file', 'ustBalbc2D.mat', 'label_column', 'end', 'label_file', 'balbcLabels.properties', 'match_scenarios', 4, 'verbose', verbose, 'see_training', true, 'color_file', 'colorsByName.properties');
    disp('run_umap Example 17 completed with no MATLAB exceptions!');
end
if ismember(18, whichOnes)
    disp('run_umap Example 18 starting...');    
    run_umap('sampleBalbcLabeled12k.csv', 'template_file', 'ustBalbc2D.mat', 'label_column', 'end', 'label_file', 'balbcLabels.properties', 'match_scenarios', 4, 'match_histogram_fig', false, 'see_training', true, 'verbose', verbose, 'false_positive_negative_plot', true);
    disp('run_umap Example 18 completed with no MATLAB exceptions!');
end
if ismember(19, whichOnes)
    disp('run_umap Example 19 starting...');    
    run_umap('s1_samusikImported_29D.csv', 'label_column', 'end', 'label_file', 's1_29D.properties', 'qf_tree', true, 'n_components', 3, 'save_template_file', 'ust_s1_samusikImported_29D_15nn_3D.mat', 'verbose', verbose);
    run_umap('s2_samusikImported_29D.csv', 'template_file', 'ust_s1_samusikImported_29D_15nn_3D.mat', 'label_column', 'end', 'label_file', 's2_samusikImported_29D.properties', 'match_scenarios', [1 2 4],  'match_histogram_fig', false, 'see_training', true, 'false_positive_negative_plot', true, 'match_supervisors', [3 1 4], 'verbose', verbose);
    disp('run_umap Example 19 completed with no MATLAB exceptions!');
end
if whichOnes==19.1
    disp('run_umap Example 19.1 starting...');    
    run_umap('s2_samusikImported_29D.csv', 'template_file', 'ust_s1_samusikImported_29D_15nn_3D.mat', 'label_column', 'end', 'label_file', 's2_samusikImported_29D.properties', 'match_scenarios', [1 2 4],  'match_histogram_fig', false, 'see_training', true, 'false_positive_negative_plot', true, 'match_supervisors', [3 1 4], 'verbose', verbose);
    disp('run_umap Example 19.1 completed with no MATLAB exceptions!');
end
if ismember(20, whichOnes)
    disp('run_umap Example 20 starting...first with NO nndescent accleration');
    tic;
    run_umap('cytofExample.csv', 'nn_descent_min_rows', 0, 'verbose', verbose);
    toc;
    disp('slow half of Example 20 completed with no MATLAB exceptions!');
    tic;
    disp('run_umap Example 20 now WITH nndescent accleration');
    run_umap('cytofExample.csv', 'verbose', verbose);
    disp('run_umap Example 20 completed with no MATLAB exceptions!');
    toc;
end
if ismember(21, whichOnes)
    run_umap('s1_samusikImported_29D.csv', 'label_column', 'end', 'label_file', 's1_29D.properties', 'n_components', 3, 'save_template_file', 'ust_s1_samusikImported_minkowski_1.80_29D_15nn_3D.mat', 'metric', 'minkowski', 'P', 1.8, 'verbose', verbose);
    run_umap('s2_samusikImported_29D.csv', 'template_file', 'ust_s1_samusikImported_minkowski_1.80_29D_15nn_3D.mat', 'label_column', 'end', 'label_file', 's2_samusikImported_29D.properties', 'match_scenarios', 4,  'see_training', true, 'match_table_fig', false, 'match_histogram_fig', false, 'false_positive_negative_plot', true, 'match_supervisors', 3, 'verbose', verbose);
    disp('run_umap Example 21 completed with no MATLAB exceptions!');
end