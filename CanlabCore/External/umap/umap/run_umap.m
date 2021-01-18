function [reduction, umap, clusterIdentifiers, extras]=run_umap(varargin)
%%RUN_UMAP reduces data matrices with 3+ parameters down to fewer
%   parameters using the algorithm UMAP (Uniform Manifold Approximation and
%   Projection).
%
%   [reduction,umap,clusterIdentifiers,extras]=RUN_UMAP(csv_file_or_data,...
%   'NAME1',VALUE1, 'NAMEN',VALUEN) 
%   
%
%   OUTPUT ARGUMENTS
%   Invoking run_umap returns these values:
%   1)  reduction, the actual data that UMAP reduces from the data 
%       specified by the input argument csv_file_or_data; 
%   2)  umap, an instance of the UMAP class made ready for the invoker 
%       to save in a MATLAB file for further use as a template.
%   3)  clusterIdentifiers, identifiers of clusters found by dbscan 
%       or DBM methods when run on the reduction of umap.
%   4)  extras, an instance of the class UMAP_extra_results.
%       See properties comments in UMAP_extra_results.m.
%
%
%   Input Argument
%   The argument csv_file_or_data is either 
%   A) a char array identifying a CSV text file containing the data 
%      to be reduced. 
%   B) the actual data to be reduced; a numeric matrix.
%
%   If A) then the CSV file must have data column names in the first line.
%   These annotate the parameters which UMAP reduces.  If B) then
%   parameter names are needed by the name-value pair argument 
%   'parameter_names'
%   when creating or running a template.
%
%   Invoke run_umap with no arguments to download CSV files that 
%   our examples below rely upon.
%
%
%   Name-Value Pair Arguments
%   Some of these are identical to those in the original Python
%   implementation documented by the inventors in their document "Basic
%   UMAP parameters" which can be retrieved at
%   https://umap-learn.readthedocs.io/en/latest/parameters.html.
%   The optional argument name/value pairs are:
%
%   Name                    Value
%
%   'min_dist'              Controls how tightly UMAP is allowed to pack
%                           points together as does the same input argument
%                           for the original implementation. Modifying this
%                           value requires the Curve Fitting Toolbox.
%                           Default is 0.3.
%
%   'n_neighbors'           Controls local and global structure as does the
%                           same input argument for the original
%                           implementation.
%                           Default is 15. 
%   
%   'metric' or 'Distance'  Controls how distance is computed in the
%                           ambient space as does the same input argument
%                           for the original Python implementation. Accepted
%                           values for metric include 
%                  'euclidean'   - Euclidean distance (default).
%                  'seuclidean'  - Standardized Euclidean distance. Each
%                                  coordinate difference between X and a
%                                  query point is scaled by dividing by a
%                                  scale value S. The default value of S
%                                  is the standard deviation computed from
%                                  X, S=NANSTD(X). To specify another
%                                  value for S, use the 'Scale' argument.
%                  'cityblock'   - City Block distance.
%                  'chebychev'   - Chebychev distance (maximum coordinate
%                                    difference).
%                  'minkowski'   - Minkowski distance. The default
%                                    exponent is 2. To specify a different
%                                    exponent, use the 'P' argument.
%                  'mahalanobis' - Mahalanobis distance, computed using a
%                                  positive definite covariance matrix C.
%                                  The default value of C is the sample
%                                  covariance matrix of X, as computed by
%                                  NANCOV(X). To specify another value for
%                                  C, use the 'Cov' argument.
%                  'cosine'      - One minus the cosine of the included
%                                  angle between observations (treated as
%                                  vectors).
%                  'correlation' - One minus the sample linear
%                                  correlation between observations
%                                  (treated as sequences of values).
%                  'spearman'    - One minus the sample Spearman's rank
%                                  correlation between observations
%                                  (treated as sequences of values).
%                  'hamming'     - Hamming distance, percentage of
%                                  coordinates that differ.
%                  'jaccard'     - One minus the Jaccard coefficient, the
%                                  percentage of nonzero coordinates that
%                                  differ.
%                  function      - A distance function specified using @
%                                  (for example @KnnFind.ExampleDistFunc). 
%                                  A distance function must be of the form
%   
%                                    function D2 = DISTFUN(ZI, ZJ),
%   
%                                  taking as arguments a 1-by-N vector ZI
%                                  containing a single row of X or Y, an
%                                  M2-by-N matrix ZJ containing multiple
%                                  rows of X or Y, and returning an
%                                  M2-by-1 vector of distances D2, whose
%                                  Jth element is the distance between the
%                                  observations ZI and ZJ(J,:).
%
%
%   'NSMethod'              Nearest neighbors search method. Values:
%                           'kdtree'
%                               Instructs run_umap to use knnsearch with
%                               a kd-tree to find nearest neighbors. 
%                               This is only valid when 'metric' is 
%                               'euclidean', 'cityblock', 'minkowski' or
%                               'chebychev'.
%                           'exhaustive'
%                               Instructs run_umap to use knnsearch with 
%                               the exhaustive search algorithm.
%                               The distance values from all the points
%                               in X to each point in Y are computed to
%                               find nearest neighbors.
%                           'nn_descent'
%                               Instructs run_umap to use KnnFind.Approximate 
%                               which uses the nn_descent C++ mex function.
%                               This tends to deliver the fastest search
%                               given certain data conditions and 
%                               name-value pair arguments.  Any speedup
%                               benefit however comes at the cost of a
%                               a slight loss of accuracy usually < 1%.
%                               This is only valid if 'metric' is NOT 
%                               'spearman', 'hamming', 'seuclidean', ...
%                               'jaccard', or a function.
%                           Default is 'nn_descent' when n_neighbors<=45
%                           and the unreduced data is not a sparse matrix
%                           and has rows>=40,000 & cols>10.
%                           If 'metric'=='mahalanobis' then this nn_descent 
%                           lower limit for rows is 5,000 and for cols is 3.
%                           Otherwise 'kdtree' is the default if cols<=10, 
%                           the unreduced data is not a sparse matrix, and 
%                           the distance metric is 'euclidean', 'cityblock', 
%                           'minkowski' or 'chebychev'.
%                           Otherwise 'exhaustive' is the default.
%
%   'P'                     A positive scalar indicating the exponent of 
%                           Minkowski distance. This argument is only 
%                           valid when 'metric' (or 'Distance') is
%                           'minkowski'. Default is 2.
%   
%   'Cov'                   A positive definite matrix indicating the 
%                           covariance matrix when computing the 
%                           Mahalanobis distance. This argument is only 
%                           valid when 'metric' (or 'Distance') is
%                           'mahalanobis'. Default is NANCOV(X).
%   
%   'Scale'                 A vector S containing non-negative values, 
%                           with length equal to the number of columns 
%                           in X. Each coordinate difference between 
%                           X and a query point is scaled by the 
%                           corresponding element of Scale. 
%                           This argument is only valid when 'Distance' 
%                           is 'seuclidean'. Default is NANSTD(X).  
%
%   'IncludeTies'           A logical value indicating whether knnsearch 
%                           will include all the neighbors whose distance 
%                           values are equal to the Kth smallest distance. 
%                           Default is false.
%
%   'BucketSize'            The maximum number of data points in the leaf 
%                           node of the kd-tree (default is 50). This 
%                           argument is only meaningful when kd-tree is 
%                           used for finding nearest neighbors.
%
%   'randomize'             true/false.  If false run_umap invokes
%                           MATLAB's "rng default" command to ensure the
%                           same random sequence of numbers between
%                           invocations.
%                           Default is false.
%   'template_file'         This identifies a .mat file with a saved
%                           instance of the UMAP class that run_umap
%                           previously produced. The instance must be be a
%                           suitable "training set" for the current "test
%                           set" of data supplied by the argument
%                           csv_file_or_data. Template processing
%                           accelerates the UMAP reduction and augments
%                           reproducibility. run_umap prechecks the
%                           suitability of the template's training set for
%                           the test set by checking the name and standard
%                           deviation distance from the mean for each
%                           parameter (AKA data column).
%                           Default is empty ([]...no template).
%
%   'see_training'          true/false to see/hide plots of both the
%                           supervising data and the supervised data with
%                           label coloring and legend. This takes effect
%                           when applying a UMAP template of a supervised
%                           reduction and when the input argument
%                           verbose='graphic'. Examples 5, 10, 11, 12 and
%                           16 apply a supervised template.  Example 16
%                           illustrates this.
%                           Default is false.
%
%   'parameter_names'       Cell of char arrays to annotate each column
%                           of the data matrix specified by
%                           csv_file_or_data. This is only needed if a
%                           template is being used or saved.
%                           Default is {}.
%                           
%   'verbose'               Accepted values are 'graphic', 'text', or
%                           'none'. If verbose='graphic' then the data
%                           displays with probability coloring and contours
%                           as is conventional in flow cytometry analysis.
%                           If method='Java' or method='MEX', then the
%                           display refreshes as optimize_layout progresses
%                           and a progress bar is shown along with a handy
%                           cancel button. If verbose='text', the progress
%                           is displayed in the MATLAB console as textual
%                           statements.
%                           Default is 'graphic'.
%                           
%   'method'                Selects 1 of 7 implementations for UMAP's
%                           optimize_layout processing phase which does
%                           stochastic gradient descent.  Accepted values
%                           are 'MEX', 'C++', 'Java', 'C', 'C vectorized',
%                           'MATLAB' or 'MATLAB Vectorized'. 'MEX' is our
%                           fastest & most recent implementation. The
%                           source
%                           umap/sgdCpp_files/mexStochasticGradientDescent.cpp
%							provides an illustration of the simplicity and
%							power of MATLAB's C++ MEX programming
%							framework.
%                           The other methods are provided for educational
%                           value.  They represent our iterative history of
%                           speeding up the slowest area of our translation
%                           from Python. We found stochastic gradient
%                           descent to be the least vectorizable. 'C' and
%                           'C vectorized', produced by MATLAB's "C coder"
%                           app, were our first attempts to accelerate our
%                           MATLAB programming. We were surprised to find
%                           our next attempt with 'Java' was faster than
%                           the code produced by C Coder.  Thus we
%                           proceeded to speed up with the 'C++' and then
%                           'MEX' implementations. 'C++', our 2nd fastest,
%                           is a separate spawned executable.  The build
%                           script and cpp source file are found in
%                           umap/sgdCpp_files.
%							Note that MathWorks open source license
%                           prevented the 'C' and 'C vectorized' modules to
%							be distributed.  You can download them too from
%                           http://cgworkspace.cytogenie.org/GetDown2/demo/umapDistribution.zip
%                           MEX, Java and C++ support the progress plots
%                           and cancellation options given by argument
%                           verbose='graphic'.
%                           Default is 'MEX'.
%
%  'progress_callback'      A MATLAB function handle that run_umap
%                           invokes when method is 'Java', 'C++', or 'MEX'
%                           and verbose='graphic'. The input/output
%                           expected of this function is
%                           keepComputing=progress_report(objectOrString).
%                           The function returns true/false to tell the
%                           reduction to keep computing or stop computing.
%                           The objectOrString argument is either a 
%                           status description before stochastic
%                           gradient descent starts or an object
%                           with properties (getEmbedding, getEpochsDone 
%                           and getEpochsToDo) which convey the state of
%                           progress. The function function 
%                           progress_report here in run_umap.m exemplifies
%                           how to write a callback.
%                           Default is the function progress_report
%                           in run_umap.m.
%
%   'ask_to_save_template'  true/false instructs run_umap to ask/not ask
%                           to save a template PROVIDING method='Java',
%                           verbose='graphic', and template_file is empty.
%                           Default is false.
%
%   'label_column'          number identifying the column in the input data
%                           matrix which contains numeric identifiers to
%                           label the data for UMAP supervision mode.
%                           If the value is 'end' then the last column in
%                           the matrix is the label_column;
%   `                       Default is 0, which indicates no label column.
%
%   'label_file'            the name of a properties file that contains
%                           the label names.  The property name/value
%                           format is identifier=false.
%                           Default is [].
%
%   'n_components'          The dimension of the space into which to embed
%                           the data.
%                           Default is 2.
%
%   'epsilon'               The epsilon input argument used by MATLAB's
%                           dbscan algorithm. 
%                           Default is 0.6.
%
%   'minpts'                The minpts input argument used by MATLAB's 
%                           dbscan.
%                           Default is 5.
%
%   'dbscan_distance'       The distance input argument used by MATLAB's
%                           dbscan.
%                           Default is 'euclidean'.
%
%
%   'cluster_output'        Allowed values: 'none', 'numeric', 'graphic'.  
%                           When the value~='none' && nargout>2 
%                           cluster results are returned in the 3rd output
%                           argument clusterIdentifiers.
%                           Default is 'numeric'.
%
%   'cluster_method_2D'     Clustering method when n_components==2.
%                           Allowed values are 'dbscan' or 'dbm'.
%                           Default is our own method 'dbm'.
%
%   'cluster_detail'        Used when (nargout>2 and the input 
%                           argument 'cluster_output'~='none') OR 
%                           'cluster_output'=='graphic'.  
%                           Allowed values are 'very low', 'low', 'medium',
%                           'high', 'very high', 'most high', 'adaptive' or
%                           'nearest neighbor' or 'dbscan arguments'
%                           if 'dbscan arguments' then run_umap uses the 
%                           input arguments 'epsilon' and 'minpts' to
%                           determine cluster detail IF the dbscan method 
%                           is needed.  If needed and 'cluster_detail'
%                           value is 'adaptive' or 'nearest neighbor' 
%                           then run_umap replaces with 'dbscan arguments'.
%                           Default is 'very high'.
%
%   'save_template_file'    Fully qualified path of the file to save the
%                           resulting UMAP object as a template.  One
%                           can also save run_umap's 2nd output argument.
%                           Default is [].
%
%   'match_supervisors'     A number indicating how to relabel data points 
%                           in the embedding data if the UMAP reduction 
%                           is guided by a template that in turn is guided 
%                           by supervisory labels.
%                           0 matches supervised and supervising data
%                             groupings by distance of medians. Supervising
%                             groupings are data points in the template's
%                             embedding that have the same supervisory
%                             label. Supervised groupings are DBM clusters
%                             in the final template-guided embedding. The
%                             publication that introduces DBM clustering is
%                             http://cgworkspace.cytogenie.org/GetDown2/demo/dbm.pdf.
%                           1 (default) matches groupings by quadratic 
%                             form dissimilarity.  The publication 
%                             that introduces QF dissimilarity is
%                             https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5818510/.
%                           2 matches supervised DBM clusters by assigning
%                             the label of the closest supervising data 
%                             point to each supervised data point and then 
%                             choosing the most frequent label in the 
%                             cluster.  Closeness is based on Euclidean 
%                             distance in the supervised and supervising
%                             embedding data spaces.
%                           3 is similar to 2 except it only uses closeness
%                             to the supervising data point to relabel the
%                             supervised data points without the aid of DBM
%                             clustering.  Thus supervised groupings in the
%                             embedding space may have small fragments
%                             of data points that occur in distant 
%                             islands/clusters of the embedding.
%
%   'match_3D_limit'        The lower limit for the # of data rows before
%                           3D progress plotting avoids supervisor label
%                           matching.  This applies only when reducing with
%                           a supervised template and the n_components>2
%                           and verbose=graphic. If > limit then supervisor
%                           matching ONLY occurs in the final plot
%                           ...otherwise supervisor label matching occurs
%                           during progress plotting before epochs finish.
%                           Default is 20000. 
%                           
%   'qf_dissimilarity'      Show QF dissimilarity scores between data
%                           groupings in the supervised and supervising 
%                           embeddings. The showing uses a sortable data 
%                           table as well as a histogram.
%                           Default is false.
%                           run_umap only consults this argument when it
%                           guides a reduction with a supervised
%                           template.
%                           
%   'qf_tree'               Show a dendrogram plot that represents the
%                           relatedness of data groupings in the
%                           supervising and supervised embeddings. The
%                           above documentation for the match_supervisors
%                           argument defines "data groupings". The
%                           publication that introduces the QF tree is
%                           https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6586874/.
%                           This uses phytree from MATLAB's Bioinformatics
%                           Toolbox, hence changing this value to true
%                           requires the Bioinformatics Toolbox.
%                           Default is false.
%                           run_umap only consults this argument when it
%                           guides a reduction with a supervised
%                           template.
%                           
%   'joined_transform'      true/false for a new transform method to avoid
%                           false positives when applying a template whose
%                           training data differs too much from test set
%                           data. This feature is not part of UMAP's
%                           original Python implementation. 
%                           Currently this not supported when
%                           method is not Java. 
%                           Default is false.
%
%   'python'                true/false to use UMAP's original
%                           implementation written in Python instead of
%                           this one written in MATLAB, C++ and Java.  The
%                           Python implementation is from Leland McInnes,
%                           John Healy, and James Melville.
%                           If true then certain arguments are ignored:
%                           joined_transform, method, verbose, 
%                           and progress_callback.
%                           Default is false.
%
%   'nn_descent_min_rows'   the # of input data rows needed before 
%                           UMAP version 2.0 engages its  NEW fast fuzzy 
%                           simplicial set processing. 
%                           Default is 40,000.
%
%   'nn_descent_min_cols'   the # of input data columns needed before 
%                           UMAP version 2.0 engages its  NEW fast fuzzy 
%                           simplicial set processing. 
%                           Default is 11
%
%   'nn_descent_transform_queue_size' 
%                           a factor of "slack" used for fuzzy simplicial
%                           set optimization with umap supervised templates.
%                           1 means no slack and 4 means 400% slack.
%                           The more slack the greater accuracy but the 
%                           less the acceleration.
%                           Default is 1.35.
%
%   'nn_descent_max_neighbors'
%                           the maximum # of n_neighbors after which the
%                           NEW acceleration of fuzzy simplicial set 
%                           processing UMAP version 2.0 becomes too slow.
%							The default is 45.
%
%                           The above 4 nn_descent* arguments guide accelerants 
%                           of fuzzy simplicial set processing released in
%                           version 2.0 of our UMAP for MatLab.  The 
%                           mex accelerants engage when metric is 
%                           anything other than mahalanobis or spearman, 
%                           when n_neighbors <= 30 and the input 
%                           data matrix has rows >= 65,000 & columns>=11.  
%                           NOTE: there could be a slight loss of
%                           accuracy (usually < 1%), so you may want
%                           to set this option off.
%                           
%   'match_scenarios'       Used for "ust/UST" scenarios where UMAP 
%                           uses a previously created supervised template.
%                           This parameter produce a table of comparision
%                           statistics for each class in the supervision.
%                           2 is the typical comparison scenario between
%                               the classification of the training set
%                               and the classification ust produces on 
%								the test set.
%                           1 compares the classification 
%                               of the training set with a prior 
%                               classification of the test set if (and
%                               only if) the test set input data has a label
%                               column denoting this prior classification.
%                           3 compares a prior classification  of the 
%                               test set with the classification which ust 
%                               produces.  The test set input data
%                               must include a label_column denoting the 
%                               prior classification.  The comparison
%                               metric used is QF dissimilarity.  See
%                               https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6586874/
%                               https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5818510/
%                           4 same as 3 except the comparison metric used
%                               as F-measure.  See 
%                               https://en.wikipedia.org/wiki/F-score
%                           Default is 0
%
% 'match_histogram_fig'     This parameter applies to match scenarios
%                           described above.  If true then run_umap shows
%                           histograms for the 2 comparison metrics
%                           of QF dissimilary and F-measure
%
% 'false_positive_negative_plot'
%                           Used for parameter match_scenarios 4 where 
%                           run_umap uses F-measure metric to compare 
%                           the classification done by
%                           ust on a test set with a prior classification 
%                           for the test set.  The test set input data
%                           must include a label_column denoting the prior
%							classification. The false +/- displays includes 
%                           several graphs to illustrate how the 
%                           ust classification compares to the prior 
%                           classification which is asssumed to be more 
%                           correct.
%                           Default is false.
%
%  'override_template_args' If false use the umap settings found in the
%                           template rather than those found in the 
%                           provided as arguments ro run_umap.  This 
%                           affects the arguments metric, P, Cov, Scale, 
%                           n_neighbors and min_dist.
%                           The arguments for n_components and parameter_names
%                           are ALWAYS ignored when using a template 
%                           (either supervised or unsupervised).
%                           Default is false.
%
%   EXAMPLES 
%   Note these examples assume your current MATLAB folder is where
%   run_umap.m is stored.
%
%   1.  Download the example CSV files and run sample10k.csv.
%
%       run_umap
%
%   2.  Reduce parameters for sample30k.csv and save as UMAP template (ut).
%
%       run_umap('sample30k.csv', 'save_template_file', 'utBalbc2D.mat');
%
%   3.  Reduce parameters for sample130k.csv using prior template.
%
%       run_umap('sample130k.csv', 'template_file', 'utBalbc2D.mat');
%
%   4.  Reduce parameters for sampleBalbcLabeled55k.csv supervised by
%       labels produced by EPP and save as a UMAP supervised template 
%       (ust), EPP is a conservative clustering technique described at
%       https://www.nature.com/articles/s42003-019-0467-6. EPP stands for
%       "Exhaustive Projection Pursuit".  By clustering exhaustively in 2
%       dimension pairs, this technique steers more carefully away from the
%       curse of dimensionality than does UMAP or t-SNE.
%
%       To use EPP you can download AutoGate from cytogenie.org which
%       contains tutorials on using EPP.
%
%       run_umap('sampleBalbcLabeled55k.csv', 'label_column', 11, 'label_file', 'balbcLabels.properties', 'save_template_file', 'ustBalbc2D.mat');
%
%   5.  Reduce parameters for sampleRag148k.csv using template that is
%       supervised by EPP.  This takes the clusters created by EPP on the
%       lymphocytes of a normal mouse strain (BALB/c) and applies them via
%       a template to a mouse strain (RAG) that has neither T cells nor B
%       cells.
%
%       run_umap('sampleRag148k.csv', 'template_file', 'ustBalbc2D.mat', 'match_supervisors', 1);
%
%   6.  Reduce parameters for sample30k.csv and return & plot cluster 
%       identifiers using density-based merging described at 
%       http://cgworkspace.cytogenie.org/GetDown2/demo/dbm.pdf.
%
%       [~,~, clusterIds]=run_umap('sample30k.csv', 'cluster_output', 'graphic');
%
%   7.  Repeat sample 2 but for 3D output and return cluster identifiers
%       and save the result as 3D template.
%
%       [~, ~, clusterIds]=run_umap('sample30k.csv', 'n_components', 3, 'save_template_file', 'utBalbc3D.mat');
%
%   8.  Repeat example 3 in 3D.
%
%       run_umap('sample130k.csv', 'template_file', 'utBalbc3D.mat');
%
%   9.  Reduce parameters and save template for sampleRagLabeled60k.csv
%       using labels produced by an expert biologist drawing manual gate 
%       sequences on lymphocyte data taken from a RAG mouse strain which 
%       has no T cells or B cells.
%
%       run_umap('sampleRagLabeled60k.csv', 'label_column', 11, 'label_file', 'ragLabels.properties', 'save_template_file', 'ustRag2D.mat');
%
%   10. Reduce parameters for lymphocyte data taken from a BALB/c mouse
%       strain using template created in example 9.  This takes the
%       clusters created on the lymphocyte data of a knockout mouse strain
%       (RAG) with no B cells or T cells and applies them to a normal mouse
%       strain (BALB/c) which has both cell types.  This illustrates logic
%       to prevent false positives for data not seen when training/creating
%       supervised templates.  Choose to re-supervise to see effect.
%
%       run_umap('sample30k.csv', 'template_file', 'ustRag2D.mat');
%
%   11. Repeat example 10 but use joined_transform.  Currently 'method'==
%       'Java' is the only support for this.
%
%        run_umap('sample30k.csv', 'template_file', 'ustRag2D.mat', 'method', 'Java', 'joined_transform', true);
%
%   12. Run example 5 again showing training/test set plot pair, QF tree
%       and QF dissimilarity plots.
%
%       run_umap('sampleRag148k.csv', 'template_file', 'ustBalbc2D.mat', 'qf_tree', true, 'qf_dissimilarity', true, 'see_training', true);
%
%   13. Compare our implementation to the original Python implementation by
%       repeating example 2 as follows.
%       
%       run_umap('sample30k.csv');
%       run_umap('sample30k.csv', 'python', true);
%
%   14. Compare our implementation with MEX method to the original Python 
%       implementation by repeating example 4 as follows.
%
%       run_umap('sampleBalbcLabeled55k.csv', 'label_column', 11, 'label_file', 'balbcLabels.properties', 'save_template_file', 'ustBalbc2D.mat');
%       run_umap('sampleBalbcLabeled55k.csv', 'label_column', 11, 'label_file', 'balbcLabels.properties', 'python', true, 'save_template_file', 'pyUstBalbc2D.mat');
%
%   15. Compare our implementation to the original Python 
%       implementation by repeating example 5 as follows.
%
%       run_umap('sampleRag148k.csv', 'template_file', 'ustBalbc2D.mat');
%       run_umap('sampleRag148k.csv', 'template_file', 'pyUstBalbc2D.mat');
%
%   16. Combining aspects of previous examples, this one creates
%       a UMAP supervised template for 3D output, then applies this
%       template to a different example. The final run_umap returns all
%       possible outputs, including the extras argument that contains
%       supervisor matching labels (1 per row of input data matrix) and
%       qf_tree and qf_dissimilarity arguments. The main plot shows
%       training/test set plot pair.
%
%       run_umap('sampleBalbcLabeled55k.csv', 'label_column', 11, 'label_file', 'balbcLabels.properties', 'qf_tree', true, 'n_components', 3, 'save_template_file', 'ustBalbc3D.mat');
%       [reduction, umap, clusterIdentifiers,extras]=run_umap('sample10k.csv', 'template_file', 'ustBalbc3D.mat', 'qf_tree', true, 'qf_dissimilarity', true, 'see_training', true, 'cluster_output', 'graphic');
%
%   17. Reduce parameters for sampleBalbCLabeled12k.csv using template that is
%       supervised by EPP and invoke match_scenarios 4 to analyze 
%       dissimilarity between subsets defined by umap supervised templates 
%       and the previously  classified subsets in the test set sample.
%
%       run_umap('sampleBalbcLabeled12k.csv', 'template_file', 'ustBalbc2D.mat', 'label_column', 'end', 'label_file', 'balbcLabels.properties',  'match_scenarios', 4,  'see_training', true);
%
%   18. Same as example 17 but add a false positive/negative plot
%
%       run_umap('sampleBalbcLabeled12k.csv', 'template_file', 'ustBalbc2D.mat', 'label_column', 'end', 'label_file', 'balbcLabels.properties',  'match_scenarios', 4, 'match_histogram_fig', false, 'see_training', true, 'false_positive_negative_plot', true);
%
%   19. Do false positive/negative test on 3 differnt types of match_supervisors: 
%           1) clusters with 'most high' cluster_detail
%           3) nearest neighbors in reduced space (2D)
%           4) nearest neighbors in non-reduced space (29D)
%       The samples are from the panoramic data set used by the Nolan Lab 
%       at Stanford University.  This is often referenced in gate
%       automation publications for flow cytometry including FlowCAP.
%       Leland McInnes references this in his UMAP publication:
%`          https://www.nature.com/articles/nbt.4314?utm_source=feedburner&utm_medium=feed&utm_campaign=Feed%3A+nbt%2Frss%2Fcurrent+%28Nature+Biotechnology+-+Issue%29
%       Nikolay Samusik references it in:
%           https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4896314/
%       The first time you run this example you must be online to
%       get the data from http://cgworkspace.cytogenie.org
%
%       run_umap('s1_samusikImported_29D.csv', 'label_column', 'end', 'label_file', 's1_29D.properties', 'n_components', 3, 'save_template_file', 'ust_s1_samusikImported_29D_15nn_3D.mat');
%       run_umap('s2_samusikImported_29D.csv', 'template_file', 'ust_s1_samusikImported_29D_15nn_3D.mat', 'label_column', 'end', 'label_file', 's2_samusikImported_29D.properties', 'match_scenarios', [1 2 4],  'match_histogram_fig', false, 'see_training', true, 'false_positive_negative_plot', true, 'match_supervisors', [3 1 4]);
%
%   20.  Check the speed difference processing 87,772 rows by 29 columns
%        with UMAP.m version 2.0. No acceleration as with previous versions
%
%        run_umap('cytofExample.csv', 'nn_descent_min_rows', 0);
%               
%        WITH acceleration
%   
%       run_umap('cytofExample.csv');
%
%   21. Determine if metric=minkowski and P=1.8 produces better
%       false positive/negative results than those reported in example 19.
%       P=1.8 causes minkowski space to combine aspects of
%       cityblock and euclidean space.
%
%       run_umap('s1_samusikImported_29D.csv', 'metric', 'minkowski', 'P', 1.8, 'label_column', 'end', 'label_file', 's1_29D.properties', 'n_components', 3, 'save_template_file', 'ust_s1_samusikImported_minkowski_1.80_29D_15nn_3D.mat');
%       run_umap('s2_samusikImported_29D.csv', 'template_file', 'ust_s1_samusikImported_minkowski_1.80_29D_15nn_3D.mat', 'label_column', 'end', 'label_file', 's2_samusikImported_29D.properties', 'match_scenarios', 4,  'see_training', true, 'match_table_fig', false, 'match_histogram_fig', false, 'false_positive_negative_plot', true, 'match_supervisors', 3);
%
%
%   NOTE that you can do supervised UMAP and templates with n_components
%       ...but we have not had time to update the 3D GUI to show where the 
%       supervised regions fall.
%
%
%   REQUIRED PATHS
%   This distribution has 2 folders:  umap and util.  
%   You must set paths to these folders plus the java inside of umap.jar.
%   Assume you have put these 2 folders under /Users/Stephen.
%   The commands that MATLAB requires would be:
%
%   addpath /Users/Stephen/umap
%   addpath /Users/Stephen/util
%   javaaddpath('/Users/Stephen/umap/umap.jar');
%
%
%   ALGORITHMS
%   UMAP is the invention of Leland McInnes, John Healy and James Melville
%   at Canada's Tutte Institute for Mathematics and Computing.  See
%   https://umap-learn.readthedocs.io/en/latest/.
%
%   AUTHORSHIP
%   Primary Developer+math lead: Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer:  Stephen Meehan <swmeehan@stanford.edu> 
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University. 
%   License: BSD 3 clause
%
%   IMPLEMENTATION NOTES
%   This is a total rewrite of the original Python implementation from
%   Leland McInnes, John Healy and James Melville. This implementation is
%   written in MATLAB, C++, and Java. The source is distributed openly on
%   MathWorks File Exchange. This implementation follows a very similar
%   structure to the Python implementation, and many of the function
%   descriptions are nearly identical. Leland McInnes has looked over it
%   and considered it "a fairly faithful direct translation of the original
%   Python code (except for the nearest neighbor search)". If you have
%   UMAP's Python implementation you can check how faithful and fast this
%   re-implementation is by using the argument python=true. When python is
%   false and method is MEX we observe superior performance on our Mac and
%   Windows laptops in most cases. For the cases of template-guided and
%   supervised parameter reduction the performance is significantly faster
%   than the Python implementation regardless of data size.
%
%   If you wish to have a simple user GUI to run these UMAP features, 
%   download AutoGate at CytoGenie.org.
%   
%
curPath=fileparts(mfilename('fullpath'));
if isequal(pwd, curPath)
    prev=which('CytoGate.m');
    if ~isempty(prev) && false
        msg(Html.Wrap(['<b><center><font size="6" color="red">Can not use '...
            'run_umap!</font></center></b><br>You must remove AutoGate '...
            'folders from MatLab path!<hr>']));
        reduction=[];
        umap=[];
        clusterIdentifiers=[];
        extras=[];
        return;
    end
end
if ~CheckUmapFolder(curPath, 'FileBasics.m', true)
    return;
end
if ~CheckUmapFolder(curPath, UmapUtil.LocateMex)
    return;
end
warning('OFF', 'MATLAB:ui:javaframe:PropertyToBeRemoved');
[args, argued]=UmapUtil.Initialize(varargin{:});
args=UmapUtil.CheckArgs(args, argued);
clusterIdentifiers=[];
globals=BasicMap.Global;
homeFolder=globals.appFolder;
try
    props=fullfile(homeFolder, 'globals.mat');
    globals.load(props);
catch
end
reductions=Map;
reductions.load(fullfile(homeFolder, 'reductions.mat'));
args.reduction=reductions.newId;
reduction=[];
umap=[];
beGraphic=strcmpi(args.verbose, 'graphic');
csv_file_or_data=args.csv_file_or_data;
save_template_file = args.save_template_file;
curAxes=[];
testSetLabels=[];
labels=[];
beQuiet=strcmpi(args.verbose, 'none');
extras=UMAP_extra_results;
firstQf=[];
if beGraphic
    xLabel=[];
    yLabel=[];
    zLabel=[];
    fig=figure('name', 'Running UMAP ...');
    extras.fig=fig;
    if args.qf_tree 
        movegui(fig, 'south')
    elseif args.matchingUmap || args.matchingUst
        movegui(fig, 'center')
    else
        movegui(fig, 'onscreen');
    end
    curAxes=axes('Parent', fig);
    if isempty(csv_file_or_data)
        [answer, cancelled]=Gui.Ask(Html.Wrap(...
            ['Options for getting examples & accelerants<br>'...
            'from the Herzenberg Lab@Stanford University.<hr>']),...
            {'<html>Download from Google Drive <i>and exit</i></html>', ...
            '<html>Download examples+accelerants then <i>keep running</i></html>!!'...
            '<html>Download accelerants <b>only</b> then <i>keep running</i></html>'}, ...
            'runUmapDownload', 'Getting our examples...', 2 );
        if ~cancelled
            if answer==3
                if UmapUtil.DownloadAdditions
                    csv_file_or_data='sample10k.csv';
                end
            elseif answer==2
                if UmapUtil.DownloadAdditions(false)
                    csv_file_or_data=downloadCsv;
                end
            elseif answer==1
                UmapUtil.GoogleDrive([], false);
            end
        end
        if isempty(csv_file_or_data)
            if beGraphic
                delete(fig);
            end
            globals.save;
            return;
        end
        if answer==2 && ~askYesOrNo(Html.Wrap([...
                'Test csv files have been downloaded:<ol>'...
                ' <li>sample10k<li>sample30k<li>sampleBalbcLabeled55k'...
                '<li>sample130k<li>sampleRag148k.csv<li>sampleRag55k.csv'...
                '</ol><br><center>Run UMAP on <b>sample10k</b> now?<hr></center>']))
            if beGraphic
                delete(fig);
            end
            globals.save;
            return;
        end
        args.csv_file_or_data=csv_file_or_data;
    end
end
args=UmapUtil.RelocateExamples(args);
csv_file_or_data=args.csv_file_or_data;
if nargout>=3 && args.n_components>2
    % check for presence of DBSCAN
    if ~Density.HasDbScan(false)
        if beGraphic
            if ~askYesOrNo(Html.WrapHr(['DBSCAN for clustering in 3+D is '...
                    '<br>not downloaded ...Continue?']))
                delete(fig);
                globals.save;
                return;
            end
        end
        dispNoDbScan;
    end
end
if ischar(csv_file_or_data)
    if ~exist(csv_file_or_data, 'file')
        if beGraphic
            delete(fig);
        end
        if startsWith(csv_file_or_data, UmapUtil.LocalSamplesFolder)
            [~,fn, ext]=fileparts(csv_file_or_data);
            if askYesOrNo(struct('icon', 'error.png', ...
                    'msg', ['<html>Can not access the example "' ...
                    fn ext '"<br><br>' globals.h2Start '<center>Open '...
                    'our shared Google Drive?</center>' globals.h2End '<hr></html>']),...
                    'Example not found')
                UmapUtil.GoogleDrive([], false, true);
                msg('All sample data is in examples sub folder');
            end
        else
            msg(['<html>The csv file <br>"<b>' ...
                globals.smallStart csv_file_or_data globals.smallEnd ...
                '</b>"<br><center><font color="red"><i>can not be found !!' ...
                '</i></font><hr></center></html>'], 25, 'center', ...
                'Error...', 'error.png');
        end
        globals.save;
        return;
    end
    [inData, parameter_names]=File.ReadCsv(csv_file_or_data);
else
    inData=csv_file_or_data;
    parameter_names=args.parameter_names;
end
newSubsetIdxs=[];
template_file=args.template_file;
[nRows, nCols]=size(inData);
if ~KnnFind.CheckDistArgs(nCols, args)
    error('Incorrect supplementary metric args (P, Cov, Scale) ');
end
sCols=num2str(nCols);
if nRows*nCols<15
    if beGraphic
        delete(fig);
    end
    if length(inData)==1
        if inData>1&&inData<=19
            if askYesOrNo(['Run example #' num2str(inData) '??'])
                run_examples(inData)
                return
            end
        end
    end
    if isempty(inData)
        if ischar(csv_file_or_data) && startsWith(csv_file_or_data, ...
                UmapUtil.LocalSamplesFolder) && ...
                askYesOrNo(Html.WrapHr('Remove corrupt example file?'))
            tempfile=[tempname '.html'];
            movefile(csv_file_or_data, tempfile);
            if askYesOrNo(['<html>' Html.H2('File removed ...')...
                    '<br>Inspect the contents in your browser?</html>'])
                Html.BrowseFile(tempfile);
            end
        else
            msgError('Table of numbers NOT provided');
        end
    else
        msgError(sprintf(...
            'Too little data: %d row(s) X  %d col(s)??!', ...
            nRows, nCols));
    end
    reduction=[];
    umap=[];
    clusterIdentifiers=[];
    extras=[];
    return
end
if strcmpi(args.label_column, 'end')
    args.label_column=nCols;
end
if args.label_column>0 
    if  args.label_column>nCols
        msg(Html.WrapHr(['The input data has ' sCols ' columns ...<br>'...
            'THUS the label_column must be >=1 and <= ' sCols]));
        assert(args.label_column>0 && args.label_column<=nCols, [...
            'label_column must be >=1 and <= ' sCols]);
    end
    labelCols=1;
    if args.matchingTestLabels
        testSetLabels=inData(:, args.label_column);
    end
    if ~isempty(template_file) 
        if args.label_column<=length(parameter_names)
            parameter_names(args.label_column)=[];
        end
    else
        labels=inData(:, args.label_column);        
    end
    inData(:, args.label_column)=[];
else
    labelCols=0;
end
firstPlot=true;
if ~isempty(template_file)
    if ischar(template_file)
        if ~exist(template_file, 'file')
            if beGraphic
                delete(fig);
            end
            msg(['<html>The template file <br>"<b>' ...
                globals.smallStart template_file globals.smallEnd ...
            '</b>"<br><center><font color="red"><i>can not be found !!' ...
            '</i></font><hr></center></html>'], 25, 'center', ...
            'Error...', 'error.png');
            globals.save;
            return;
        end
        if isempty(parameter_names)
            warning('No parameter names supplied; template may become mismatched. Good luck!');
        elseif length(parameter_names)~=size(inData, 2)
            showMsg(Html.WrapHr(sprintf(['<b>Can not create '...
                'or use template</b> ...<br>'...
                '%d parameter_names... but data has %d parameters?'], ...
                length(parameter_names), size(inData,2))));
            if beGraphic
                delete(fig)
            end
            globals.save;
            return;
        end
        [umap, ~, canLoad, reOrgData, paramIdxs]=...
            Template.Get(inData, parameter_names, ...
            template_file, 3);
        if ~isempty(reOrgData)
            % column label order differed
            inData=reOrgData;
            if ~isempty(parameter_names) && ~isempty(paramIdxs)
                parameter_names=parameter_names(paramIdxs);
            end
        end
    elseif isa(template_file, 'UMAP')
        umap=template_file;
        canLoad = true;
    end
    if isempty(umap)
        if ~canLoad
            if beGraphic
                showMsg(Html.WrapHr(['No template data found in <br>"<b>', ...
                    template_file '</b>"']));
            else
                disp(['No template data found in ' template_file]);
            end
        end
        if beGraphic
            delete(fig);
        end
        globals.save;
        return;
    else
        args.n_components=umap.n_components;
        if ~isempty(umap.supervisors)
            sprv=umap.supervisors;
            sprv.verbose=args.verbose;
            sprv.overrideColors(args.color_file, beQuiet);
            if beQuiet
                LabelBasics.Frequency(sprv.labels, sprv.labelMap, true)
            else
                LabelBasics.Frequency(sprv.labels, sprv.labelMap, true, [])
            end
            
            sprv.description=args.description;
            sprv.context=args.context;
            
            %Connor's NEW joined_transform immunizes reduction from
            %false positives if items in the test set are TOO different
            %from the training set
            
            if ~args.joined_transform
                [percNewSubsets, unknownIdxs]=...
                    Template.CheckForUntrainedFalsePositives(umap, inData);
                if percNewSubsets>13 && beGraphic
                    [choice, cancelled]=Template.Ask(percNewSubsets);
                    if cancelled
                        if beGraphic
                            delete(fig);
                        end
                        globals.save;
                        return;
                    end
                    if choice==2
                        umap.clearLimits;
                        newSubsetIdxs=unknownIdxs;
                        template_file=[];
                    end
                end
            end
            sprv.initClustering(args.cluster_detail{1}, ...
                args.cluster_method_2D, args.minpts, ...
                args.epsilon, args.dbscan_distance);
            sprv.initPlots(args.contour_percent);
        end
    end
else
    umap = UMAP;
end

isSupervising=isprop(umap, 'supervisors') && ~isempty(umap.supervisors);
if isSupervising
    if args.n_components==2
        progressMatchType=0;
    else
        limit=args.match_3D_limit;
        if limit<nRows
            progressMatchType=-1;
        else
            progressMatchType=3;
        end
    end
else
    if any(args.match_supervisors>0)
        if argued.contains('match_supervisors')
            if ~isequal(1, args.match_supervisors)
                warning('match_supervisors only affects supervised template reduction');
                args.match_supervisors=1;
            end
        end
    end
end
UmapUtil.SetArgsTemplateCanOverride(umap, args, argued, parameter_names);
umap.n_epochs=args.n_epochs;
umap.nn_descent_min_rows=args.nn_descent_min_rows;
umap.nn_descent_min_cols=args.nn_descent_min_cols;
umap.nn_descent_max_neighbors=args.nn_descent_max_neighbors;
umap.nn_descent_transform_queue_size=args.nn_descent_transform_queue_size;
if strcmpi('Java', args.method)
    if ~initJava
        args.method='MEX';
        showMsg(Html.WrapHr('Could not load umap.jar for Java method'), ...
            'Problem with JAVA...', 'south west', false, false);
    end
elseif strcmpi('Mex', args.method)
    exeSGD=fullfile(curPath, UmapUtil.LocateMex('sgd'));
    exeNN=fullfile(curPath, UmapUtil.LocateMex);
    if ~exist(exeSGD, 'file') || ~exist(exeNN, 'file')
        UmapUtil.OfferFullDistribution(true)
        globals.save;
        if ~exist(exeSGD, 'file') || ~exist(exeNN, 'file')
            if ~askYesOrNo(['<html>Continue more slowly <br>'...
                    'without accelerants?<hr></html>'])
                if beGraphic
                    close(fig);
                end
                return;
            end
        end
    end
end
        
method=umap.setMethod(args.method);
umap.verbose=~beQuiet;
umap.randomize=args.randomize;
tick=tic;

labelMap=[];
nParams=length(parameter_names);

good=nParams==0||nParams==nCols || (args.label_column>0 &&...
    (nParams==nCols-1 || nParams==nCols));
if ~good
    if args.label_column>0
        preAmble=sprintf(['# data columns=%d, # parameter_names=%d '...
            'since label_column=%d <br># parameter_names must be '...
            '%d or %d '],  nCols, nParams, args.label_column, ...
            nCols, nCols-1);
    else
        preAmble=sprintf(['# of data columns(%s) must equal '...
            '# of parameter_names(%d)'], sCols, nParams);
    end
    msg(Html.WrapHr(preAmble));
    assert(nParams==0||nParams==nCols || (args.label_column>0 &&...
        (nParams==nCols-1 || nParams==nCols)), preAmble);    
end
if ~isempty(newSubsetIdxs)
    hasLabels=true;
    labelCols=0;
    [labels, labelMap]=resupervise(umap, inData, newSubsetIdxs);
    nLabels=length(unique(labels));
elseif args.label_column>0 && isempty(template_file) && ~args.matchingUmap
    hasLabels=true;
    labelCols=1;
    good=args.label_column>0 && args.label_column<=nCols;
    if ~good
        msg(Html.WrapHr(['The input data has ' sCols ' columns ...<br>'...
            'THUS the label_column must be >=1 and <= ' sCols]));
        assert(args.label_column>0 && args.label_column<=nCols, [...
            'label_column must be >=1 and <= ' sCols]);
        
    end
    nLabels=length(unique(labels));
    if nLabels > .5*nRows
        preAmble='%d is a LOT of distinct labels for a %dx%d matrix!';
        msg(['WARNING:  ' sprintf(preAmble, nLabels, nRows, nCols)]);
        warning(preAmble, nLabels, nRows, nCols);
    end
    if args.label_column<=nParams
        parameter_names(args.label_column)=[];
    end
    umap.dimNames=parameter_names;
    nLabels=length(unique(labels));
else
    hasLabels=false;
    nLabels=0;
end
if args.label_column>0 
    if isSupervising && ~isempty(args.label_file) && isempty(testSetLabels)
        if args.label_column==0
            warning(['label_file has no effect when reducing '...
                'with supervised template without test set labels']);
        elseif isempty(find(args.match_scenarios==1, 1)) .... 
                && isempty(find(args.match_scenarios>2,1))
            warning(['test set labels only needed for umap'...
                ' supervised templates IF specifying '...
                'match_scenarios 1 3 or 4 ']);
        end
    else
        if ~isempty(testSetLabels)
            [labelMap, halt]=getLabelMap(testSetLabels);
        else
            [labelMap, halt]=getLabelMap(labels);
            if hasLabels
                ColorsByName.Override(labelMap, args.color_file, beQuiet);
            end
        end
        if halt
            delete(fig);
            globals.save;
            return;
        end
        if ~isempty(labelMap) && ~args.buildLabelMap ...
                && isSupervising && ~isempty(testSetLabels)
            testSetLabels=LabelBasics.RelabelIfNeedBe(testSetLabels, ...
                sprv.labelMap, labelMap);
        end

    end
end
nanRows=any(isnan(inData'));
badRows=sum(nanRows);
if badRows>0
    if beGraphic
        if askYesOrNo(Html.WrapHr(['Data matrix has ' ...
                String.Pluralize2('row', badRows) ...
                'with NAN values <br>which cause odd effects on '...
                'UMAP!<br>Try to remove nan values?']))
            inData=inData(~nanRows,:);
        end        
    else
        warning(['Data matrix has ' ...
                String.Pluralize2('row', badRows) ...
                'with NAN values!']);
        inData=inData(~nanRows,:);
    end
    if any(isnan(inData(:)))
        showMsg(Html.WrapHr(['Sorry...<br>can not proceed<br>'...
            '<br>NAN values exist... SIGH!']));
        globals.save;
        return;
    end
end
args.hiD=nCols-labelCols;
info=[String.encodeInteger(nRows) 'x' String.encodeInteger(nCols-labelCols)];
if ischar(csv_file_or_data)
    [~, fileName]=fileparts(csv_file_or_data);
    info=['UMAP on ' String.ToLaTex(fileName) ', ' info];
else
    info=['[UMAP on ' info];
end
if args.python
    info=[info ', Python'];
else
    info=[info ', ' method];
end
if beGraphic
    set(fig, 'NumberTitle', 'off', 'name', info );
    drawnow;
end
pause(.01);
info2=['(optimize\_layout method=' method ')'];
if strcmpi(method, 'C++')
    if ~StochasticGradientDescent.IsAvailable
        if ~askYesOrNo(Html.Wrap(...
                ['This C++ executable is missing or corrupt:'...
            '<br>"<b>' StochasticGradientDescent.GetCmd '</b>"'...
            '<br><br>Maybe try rebuilding by changing clang++ '...
            'to g++ in the build scripts in the same folder...<br>'...
            '<br><center>Try <b>method=Java</b> instead?</center><hr>']))
            return;
        end
        method='MEX';
    end
end
if strcmpi(method, 'Java') || strcmpi(method, 'C++') || strcmpi(method, 'MEX')
    if beGraphic
        if isempty(args.progress_callback)
            umap.progress_callback=@(javaObject)progress_report(javaObject);
        else
            umap.progress_callback=args.progress_callback;
        end
        set(fig, 'NumberTitle', 'off', 'name', info);
        try
            nTh=edu.stanford.facs.swing.StochasticGradientDescent.EPOCH_REPORTS+3;
            figure(fig);
            if args.qf_tree
                puLocation='north++';
            else
                if args.see_training
                    puLocation='south west++';
                else
                    puLocation='south++';
                end
            end
            path=BasicMap.Path;
            pu=PopUp(Html.WrapHr(sprintf(['Using UMAP to reduce '...
                ' <b>%d</b> parameters down to ' ...
                num2str(args.n_components) '...'], nCols-labelCols)), ...
                puLocation, 'Reducing parameters...', false, true, ...
                fullfile(path, 'genieSearch.png'));
            pu.initProgress(nTh);
            pu.pb.setStringPainted(true);
            pu.setTimeSpentTic;
            drawnow;
        catch
            args.method='MEX';
            method=umap.setMethod(args.method);
            showMsg(Html.WrapHr(['Could not load umap.jar for Java method'...
                '<br><br>Switching optimize_layout method to "MEX" ']), ...
                'Problem with JAVA...', 'south west', false, false);
        end
    end
end
tic;
if beGraphic
    if ispc
        left=.21;
        width=.61;
        height=.145;
        lbl=annotation(fig, 'textbox','String', {['\color{blue}Running '...
            info], ['\fontsize{9}' info2]}, 'units', 'normalized', ...
            'position', [left .4 width height], 'fontSize', 11, ...
            'HorizontalAlignment', 'center');
    else
        left=.21;
        width=.61;
        height=.131;
        lbl=annotation(fig, 'textbox','String', {['\color{blue}Running '...
            info], ['\fontsize{11}' info2]}, 'units', 'normalized', ...
            'position', [left .4 width height], 'fontSize', 13, ...
            'HorizontalAlignment', 'center');
    end
    updatePlot;
end
paramAnnotation=[];
if ~beQuiet
    strMetric=umap.metric;
    if ~ischar(strMetric)
        strMetric='custom';
    end
    if ~isempty(umap.dist_args)
        if strcmpi(strMetric,'minkowski')
            strMetric=[strMetric ' P=' String.encodeBank(umap.dist_args)];
        elseif strcmpi(strMetric,'mahalanobis')
            strMetric=[strMetric ' Cov']; 
        elseif strcmpi(strMetric, 'seuclidian')
            strMetric=[strMetric ' Scale'];
        end
    end
    txt=sprintf(['n\\_neighbors=\\color{blue}%d\\color{black}, '...
        'min\\_dist=\\color{blue}%s\\color{black}, '...
        'metric=\\color{blue}%s\\color{black},'...
        'randomize=\\color{blue}%d\\color{black}, '...
        'labels=\\color{blue}%d'], ...
        umap.n_neighbors, num2str(umap.min_dist), strMetric,...
        umap.randomize, nLabels); 
    if beGraphic
        paramAnnotation=annotation(fig, 'textbox','String', txt,...
            'units', 'normalized', 'position', [.03 .94 .92 .05],...
            'fontSize', 9, 'HorizontalAlignment', 'center');
        drawnow;
    end
end

%READY TO START REDUCING HiD to LoD!!
if isSupervising
    args.reductionType=UMAP.REDUCTION_SUPERVISED_TEMPLATE;
else    
    if isempty(template_file)
        if hasLabels
            args.reductionType=UMAP.REDUCTION_SUPERVISED;
        else
            args.reductionType=UMAP.REDUCTION_BASIC;
        end
    else
        args.reductionType=UMAP.REDUCTION_TEMPLATE;
    end
    if args.matchingUst
        warning('qf_dissimilarity=true and match_scenarios=1 or 2 ONLY affects supervised template reduction')
        args.match_scenarios(args.ustMatches)=[];
        args.matchingUst=false;
    end
end
extras.args=args;
strReduction=[UmapUtil.GetReductionLongText(args.reductionType) ' reduction'];
reportProgress(['Running ' strReduction ', v' UMAP.VERSION], true);
%Map.SetStruct(reductions, args.reduction, args);
reductions.save;
if ~beQuiet
    disp(UMAP.REDUCTION_TYPES);
end
if ~isempty(template_file)
    if ~isempty(umap.pythonTemplate)
        args.python=true;
    end
    if ~args.python 
        if ~args.joined_transform
            reduction=umap.transform(inData);
        else
            reduction=umap.transform2(inData);
        end
    else
        if isempty(umap.pythonTemplate) || ~exist(umap.pythonTemplate, 'file')
            reduction=[];
            msg(Html.WrapHr('Python template file not found'),  8,...
                'south', 'Error...', 'error.png');
        else
            inFile=[tempname '.csv'];
            reduction=UmapPython.Go(inData,inFile, [], ...
                lower(umap.metric), umap.n_neighbors, ...
                umap.min_dist, umap.n_components, [], ...
                umap.pythonTemplate, args.verbose);
        end
    end
else
    if ~args.python
        if ~hasLabels
            reduction = umap.fit_transform(inData);
        else
            reduction = umap.fit_transform(inData, labels);
            if ~isempty(reduction)
                if ~isempty(labelMap)
                    umap.setSupervisors(labels, labelMap, curAxes);
                end
            end
        end
    else
        inFile=[tempname '.csv'];
        if ~hasLabels
            labels=[];
        end
        reduction=UmapPython.Go(inData,inFile, [], ...
            lower(umap.metric), umap.n_neighbors, ...
            umap.min_dist, umap.n_components, labels, ...
            [], args.verbose);
        pythonTemplate=fullfile(fileparts(inFile), ...
            [UmapPython.PYTHON_TEMPLATE '.umap']);
        if ~isempty(reduction)
            umap.embedding=reduction;
            umap.raw_data=inData;
            if ~isempty(labelMap)
                umap.setSupervisors(labels, labelMap, curAxes);
            end
        end
    end
    if ~isempty(umap.supervisors)
        sprv=umap.supervisors;
        sprv.description=args.description;
    end
end
if ~isempty(paramAnnotation)
    try
        set(paramAnnotation, 'visible', 'on');
    catch
    end
end
if ~isempty(reduction)
    if beGraphic
        figure(fig);
        if ~exist('pu', 'var')
            pu=PopUp('Updating plot', 'north east', [], false);
        end
        delete(lbl);
        if strcmpi(method, 'Java') || strcmpi(method, 'C++') || strcmpi(method, 'MEX')
            pu.pb.setString('All done');
        end
        updatePlot(reduction, true)
        annotation(fig, 'textbox', 'String', ['Compute time=\color{blue}' ...
            String.MinutesSeconds(toc(tick))],'units', 'normalized', ...
            'position', [.65 .01 .33 .05], 'fontSize', 9)
        if isempty(template_file) && args.ask_to_save_template
            if isequal('Yes', questdlg({'Save this UMAP reduction', ...
                    'as template to accelerate reduction', ...
                    'for compatible other data sets?'}))
                if length(parameter_names)~=size(inData, 2)
                    showMsg(Html.WrapHr(sprintf(['<b>Can not create '...
                        'template</b> ...<br>'...
                        '%d parameter_names ...but data has %d parameters?'], ...
                        length(parameter_names), size(inData,2))));
                else
                    umap.prepareForTemplate(curAxes);
                    if ischar(csv_file_or_data)
                        Template.Save(umap, csv_file_or_data);
                    else
                        Template.Save(umap, fullfile(pwd, 'template.csv'));
                    end
                end
            end
        end
    end
    testBasicReduction;
    if ~beQuiet
        fprintf('UMAP reduction finished (cost %s)\n', ...
            String.MinutesSeconds(toc(tick)));
    end
    reportProgress(['Finished ' strReduction]);
    if beGraphic
        if exist('pu', 'var') && isa(pu, 'PopUp')
            pu.stop;
            pu.dlg.dispose;
        end
    else
        if isSupervising
            if ~beQuiet
                pu=PopUp('Matching results', 'west+', [], false, [], [],...
                    false, args.parent_popUp);
            else
                pu=[];
            end
            if nargout>3
                if ~beQuiet
                    disp('Setting supervisor labels');
                end
                nSupervisors=sprv.computeAndMatchClusters( ...
                    reduction, args.match_supervisors(1), pu);
                if nSupervisors>0
                    extras.supervisorMatchedLabels=sprv.supervise(...
                        reduction, false, args.match_supervisors(1));
                end
            end
            doQfs(reduction);
            if ~beQuiet
                pu.close;
            end
        end
    end
else
    msg('Parameter reduction was cancelled or not done');
    if exist('pu', 'var') && isa(pu, 'PopUp')
        pu.stop;
        pu.dlg.dispose;
    end
end

if (nargout>1 || ~isempty(save_template_file)) && ~isempty(reduction)
    if beGraphic
        umap.prepareForTemplate(curAxes);
    else
        umap.prepareForTemplate;
    end
    if args.python
        if  ~isempty(save_template_file)
            [f1, f2]=fileparts(save_template_file);
            if isempty(f1)
                f1=pwd;
                if ischar(args.csv_file_or_data) 
                    f3=fileparts(args.csv_file_or_data);
                    if ~isempty(f3)
                        f1=f3;
                    end
                end
            end
        elseif ischar(csv_file_or_data)
            [f1, f2]=fileparts(csv_file_or_data);
            if isempty(f1)
                f1=pwd;
            end
            f2=[f2 '.umap'];
        else
            f1=[];
        end
        if ~isempty(f1)
            umap.pythonTemplate=fullfile(f1, [f2 '.python']);
            movefile(pythonTemplate, umap.pythonTemplate, 'f');
        end
    end
    if  ~isempty(save_template_file)
        f1=fileparts(save_template_file);
        if isempty(f1)
            f1=pwd;
            if ischar(args.csv_file_or_data)
                f2=fileparts(args.csv_file_or_data);
                if ~isempty(f2)
                    f1=f2;
                end
            end
            save_template_file=fullfile(f1, save_template_file);
        end
        save(save_template_file, 'umap');
    end
end    

if nargout>2 || ~strcmpi(args.cluster_output, 'none')
    if nargout<3 && ~strcmpi(args.cluster_output, 'graphic')
        warning('No clusterIdentifiers output argument');
    elseif ~strcmpi(args.cluster_output, 'ignore')
        clusterIdentifiers=doClusters(reduction);
        if isempty(clusterIdentifiers)
            dispNoDbScan;
        end
    end
end
if ~beGraphic
    if nargout<4 %don't delete figures if expecting them in extras
        extras.closeMatchFigs;
        extras.closeTreeFigs;
    end
end
globals.save;

    function testBasicReduction
        if isempty(testSetLabels)
            return
        end
        if strcmp(args.reductionType, UMAP.REDUCTION_SUPERVISED_TEMPLATE)
            return;
        end
        createdPu=false;
        if ~exist('pu', 'var')
            if beQuiet
                pu=[];
            else
                createdPu=true;
                pu=PopUp('Matching results', 'west+', [], false, ...
                    [],[],false, args.parent_popUp);
            end
        end
        nCluDtls=length(args.cluster_detail);
        for c=1:nCluDtls
            for ms=1:length(args.match_scenarios)
                scenario=args.match_scenarios(ms);
                reportProgress(sprintf(...
                    'Match clusters=%s, scenario=%d:"%s"', ...
                    args.cluster_detail{c},scenario, ...
                    UmapUtil.GetMatchScenarioText(scenario, ...
                    args.reductionType)));        
                if scenario==3
                    matchStrategy=1; % match by qf dissimilarity
                elseif scenario==4
                    matchStrategy=2; % match by F measure overlap
                else
                    continue;
                end
                [clusterIds, numClusters, density]=UmapUtil.Cluster(...
                    reduction, args.cluster_detail{c}, pu, ...
                    args.cluster_method_2D, args.minpts, ...
                    args.epsilon, args.dbscan_distance);
                [qft, tNames]=UmapUtil.Match(args, inData, ...
                    testSetLabels, labelMap, clusterIds,  ...
                    args.cluster_detail{c}, matchStrategy, false, pu);
                if ~isempty(qft)
                    extras.qfd{end+1}=qft;
                    if beGraphic
                        if args.match_table_fig
                            if ~isempty(qft.fig)
                                set(qft.fig, 'visible','on');
                            end
                        end
                        if args.match_histogram_fig
                            if ~isempty(qft.qHistFig)
                                set(qft.qHistFig, 'visible','on');
                            end
                            if ~isempty(qft.fHistFig)
                                set(qft.fHistFig, 'visible','on');
                            end
                        end
                        matchedLbls=UmapUtil.GetMatches(reduction, qft.qf, ...
                            tNames, labelMap, density, clusterIds, numClusters);
                        Supervisors.Plot(reduction, matchedLbls,...
                            labelMap, nCols-labelCols, umap, curAxes, ...
                            true, false, args.contour_percent);
                        if args.n_components==2
                            UmapUtil.DrawClusterBorders(curAxes, density, ...
                                [.5 0 .65]);
                        end
                    end
                end
            end
        end
        if ~isempty(extras.qfd)
            extras.doMatchOutput(nCols-labelCols);
            if beGraphic
                extras.seeMatches(args.match_html);
            end
        end
        if createdPu
            pu.close;
        end
    end


    function [map, halt]=getLabelMap(lbls)
        halt=false;
        map=[];
        if isempty(args.label_file)
            warning(['label_column without label_file '...
                'to match/supervise, will use default names/colors']);
            args.buildLabelMap=true;
        end
        if args.buildLabelMap
        else
            if exist(args.label_file, 'file')
                map=File.ReadProperties( args.label_file);
                if isempty(map)
                    problem='load';
                end
            elseif ~isempty(args.label_file)
                problem='find';
            end
            if isempty(map)
                if askYesOrNo(['<html>Can not ' problem ' the '...
                        ' label file <br><br>"<b>' globals.smallStart ...
                        args.label_file globals.smallEnd '</b>"<br><br>'...
                        '<center>Use default names & colors?</center>'...
                        '<hr></html>'], 'Error', 'north west', true)
                    args.buildLabelMap=true;
                else
                    halt=true;
                end
            end
        end
        if args.buildLabelMap
            map=java.util.Properties;
            u=unique(lbls)';
            nU=length(u);
            for i=1:nU
                key=num2str(u(i));
                map.put(java.lang.String(key), ['Subset #' key]);
                map.put([key '.color'], num2str(Gui.HslColor(i, nU)));
            end
        end
        if args.color_defaults
            ColorsByName.Override(map, args.color_file, beQuiet);
        end
    end

    function clues=doClusters(data)
        if isempty(data)
            clues=[];
        else
            [mins, maxs]=Supervisors.GetMinsMaxs(data);
            [~, clues]=Density.FindClusters(reduction, ...
                args.cluster_detail{1}, args.cluster_method_2D, ~beQuiet,...
                args.epsilon, args.minpts, args.dbscan_distance, mins, maxs);
            if strcmpi(args.cluster_output, 'graphic')
                if ~isempty(clues)
                    if ~exist('xLabel', 'var')
                        dimInfo=sprintf('  %dD\\rightarrow%dD', nCols-labelCols, ...
                            args.n_components);
                        xLabel=['UMAP-X' dimInfo];
                        yLabel=['UMAP-Y' dimInfo];
                        zLabel=['UMAP-Z' dimInfo];
                    end
                    cp=ClusterPlots.Go([], data, clues, [], xLabel, ...
                        yLabel, zLabel, true, [], false, true, false);
                    if nCols-labelCols==2
                        if isequal('dbscan', args.cluster_method_2D)
                            clue='dbscan';
                        else
                            clue='dbm';
                        end
                    else
                        clue='dbscan';
                    end
                    annotateClues(get(cp.ax, 'Parent'), args.cluster_detail{1}, clue, args.epsilon, args.minpts, args.dbscan_distance);
                end
            end
        end
    end

    function lbl=annotateClues(fig, detail, clue, epsilon, minpts, dist)
        X=.005;
        Y=.872;
        W=.58;
        H=.115;
        [epsilon, minpts]=Density.GetDbscanParameters(detail, ...
            epsilon, minpts);
        info=['clue method="', clue '", detail="' detail '"'];
        info2=['epsilon=' num2str(epsilon) ', minpts=' num2str(minpts) ...
            ', dbscan distance="' dist '"'];
        lbl=annotation(fig, 'textbox','String', ...
            {['\color{blue} ' info], ['\fontsize{10} ' info2]}, ...
            'units', 'normalized', 'position', [X Y W H],...
            'fontSize', 11, 'HorizontalAlignment', 'center');
    end

    function updatePlot(data, lastCall)
        doingDensity3D=false;
        labelsDone=true;
        if nargin<2
            lastCall=false;
        end
        if nargin>0
            if isempty(xLabel)
                dimInfo=sprintf('  %dD\\rightarrow%dD', nCols-labelCols, ...
                    args.n_components);
                xLabel=['UMAP-X' dimInfo];
                yLabel=['UMAP-Y' dimInfo];
                if args.n_components>2
                    zLabel=['UMAP-Z' dimInfo];
                end
            end
            if args.n_components>2
                nD=size(data, 2);
                assert(nD==args.n_components);
                if ~plotLabels(data, lastCall)
                    labelsDone=false;
                    if args.frequencyDensity3D
                        doingDensity3D=true;
                        Gui.PlotDensity3D(curAxes, data, 64, 'iso',...
                            xLabel, yLabel, zLabel);
                    else
                        Gui.PlotNeighDist3D(curAxes, data, ...
                            args.n_neighbors);
                    end
                end
                if args.n_components>3
                    title(curAxes, ['NOTE:  Only 3 of \color{red}' ...
                        num2str(args.n_components) ...
                        ' dimensions being shown...']);
                end
            else
                if plotLabels(data, lastCall)
                    if ~isempty(paramAnnotation)
                        set(paramAnnotation, 'visible', 'off');
                    end
                else
                    labelsDone=false;
                    if lastCall
                        ProbabilityDensity2.Draw(curAxes, data, ...
                            true, true, true, .05);
                    else
                        ProbabilityDensity2.Draw(curAxes, data);
                    end
                end
            end
        end
        if ~labelsDone
            if nargin>0
                umap.adjustLims(curAxes, data );
            else
                umap.adjustLims(curAxes);
            end
            xlabel(curAxes, xLabel);
            ylabel(curAxes, yLabel);
            if args.n_components>2
                zlabel(curAxes, zLabel);
            end
            grid(curAxes, 'on')
            set(curAxes, 'plotboxaspectratio', [1 1 1])
            if lastCall
                if ~doingDensity3D
                    Gui.StretchLims(curAxes, data, .04);
                end
            end
        end
        if lastCall
            if args.n_components==2
                if isSupervising
                    if ~hasLabels
                        sprv.drawClusterBorders(curAxes);
                    end
                end
            end
        end
        drawnow;
    end

    function ok=plotLabels(data, lastCall)
        if hasLabels
            ok=true;
            if lastCall
                Supervisors.Plot(data, labels, labelMap, nCols-labelCols, ...
                    umap, curAxes, true, false, args.contour_percent);
                Gui.StretchLims(curAxes, data, .04);
                if ~isempty(sprv)
                    sprv.prepareForTemplate;
                    if args.qf_tree
                        sprv.inputData=umap.raw_data;
                        [~,qft]=sprv.qfTreeSupervisors(true, ...
                            [], 'UMAP training set');
                        if ~isempty(qft) && ~isempty(qft.fig)
                            extras.qft=qft;
                        end
                    end
                end
            else
                Supervisors.Plot(data, labels, labelMap, nCols-labelCols, ...
                    umap, curAxes, false, false, args.contour_percent);
            end
        elseif isSupervising
            ok=true;
            if lastCall
                if args.match_supervisors(1)==4
                    sprv.inputData=umap.raw_data;
                    sprv.computeNearestNeighborsUnreduced(inData, pu)
                end
                if isempty(sprv.embedding)
                    sprv.embedding=umap.embedding;
                end
                if firstPlot && args.see_training
                    [curAxes, ~,~,extras.supervisorMatchedLabels]...
                        =sprv.plotTrainingAndTestSets(...
                        data, curAxes, umap, pu, args.match_supervisors(1), ...
                        true);
                else
                    [~,extras.supervisorMatchedLabels, firstQf]=...
                        sprv.plotTestSet(umap, curAxes, data, pu, ...
                        args.match_supervisors(1), true, true);
                end
                doQfs(data);
                drawnow;
            else
                if firstPlot && args.see_training
                    curAxes=sprv.plotTrainingAndTestSets(...
                        data, curAxes, umap, pu, progressMatchType, false);
                else
                    sprv.plotTestSet(umap, curAxes, ...
                        data, pu, progressMatchType, false);
                end
                firstPlot=false;
            end
        else
            ok=false;
        end
    end

    function doQfs(data)
        if args.qf_tree || all(args.match_scenarios>0)
            creatingPu=~exist('pu', 'var') || isempty(pu);
            if creatingPu
                if beQuiet
                    pu=[];
                else
                    pu=PopUp('Matching results', 'north west+', [], false);
                end
            end
            sprv.inputData=umap.raw_data;
            cascading={};
            hasFig=exist('fig', 'var');
            if hasFig
                scrFig=fig;
            else
                scrFig=[];
            end
            if args.qf_tree
                if ~beQuiet
                    disp('Computing QF tree(s)');
                end
                [~,qft]=sprv.qfTreeSupervisors(false, pu);
                if ~isempty(qft) && ~isempty(qft.fig)
                    extras.qftSupervisors=qft;
                    if beGraphic
                        cascading{end+1}=qft.fig;
                        Gui.CascadeFigs(cascading, false, true, 15, 2, ...
                            true, false, scrFig);
                        if hasFig
                            figure(scrFig);
                        end
                    end
                    [~,qft]=sprv.qfTreeSupervisees(data, ...
                        inData, false, pu);
                    if ~isempty(qft)
                        extras.qft=qft;
                        if beGraphic
                            cascading{end+1}=qft.fig;
                            Gui.CascadeFigs(cascading, false, true, 40, 2, ...
                                true, false, scrFig);
                        end
                    end
                end
            end
            if all(args.match_scenarios>0)
                if ~beQuiet
                    matchProgress('Matching UST results');
                end
                if hasFig
                    figure(scrFig);
                end
                matchTypes=unique(args.match_supervisors);
                matchTypes=[args.match_supervisors(1) ...
                    matchTypes(matchTypes ~= args.match_supervisors(1))];
                nMatchTypes=length(matchTypes);
                nCluDtls=length(args.cluster_detail);
                for c=1:nCluDtls
                    for mi=1:nMatchTypes
                        if c>1 || (mi>1 || ~beGraphic)
                            matchType=matchTypes(mi);
                            if c>1
                                if matchType>=3 % nearest neighbor no clustering
                                    continue;
                                elseif mi==1
                                    sprv.initClustering(...
                                        args.cluster_detail{c}, ...
                                        args.cluster_method_2D, ...
                                        args.minpts, ...
                                        args.epsilon, ...
                                        args.dbscan_distance);
                                    sprv.computeAndMatchClusters(data,...
                                        matchType, pu);
                                end
                            end
                            if ~beQuiet
                                if matchType<3
                                    word=[' @ "' args.cluster_detail{c} '"'];
                                else
                                    word='';
                                end
                                matchProgress(...
                                    sprintf('New match type %s %s', ...
                                    UmapUtil.GetMatchTypeLongText(matchType, ...
                                    args.reductionType, args.n_components, ...
                                    nCols-labelCols), word), ...
                                    sprintf('Match type=%d, clu=%s', ...
                                    matchType, args.cluster_detail{c}));
                            end
                            if matchType==4
                                sprv.changeMatchType(inData, ...
                                    matchType, pu);
                            else
                                sprv.changeMatchType(data, ...
                                    matchType, pu);
                            end
                        end
                        for msi=1:length(args.match_scenarios)
                            scenario=args.match_scenarios(msi);
                            matchStrategy=1;
                            report=sprintf(...
                                'Match scenario=%d:"%s"', ...
                                scenario, UmapUtil.GetMatchScenarioText(...
                                scenario, args.reductionType));                            
                            if scenario==2
                                %match training set prior classification 
                                %to ust trained re-classification of test set
                                [~,qfd]=sprv.qfDissimilarity(data, ...
                                    inData, false, pu, firstQf);
                                firstQf=[];
                                reportProgress(report);
                            else
                                if scenario==1
                                    %match training set classification to 
                                    %prior classification of test set
                                    %ONLY needed once
                                    if mi>1 || c>1
                                        continue;
                                    end
                                    withTraining=true;
                                else
                                    %match ust re-classification of test 
                                    %set to prior classification of test set
                                    
                                    withTraining=false;
                                    if scenario==4
                                        matchStrategy=2;
                                    end
                                end
                                if isempty(testSetLabels)
                                    warning(...
                                    ['Can not do qf dissimilarity if'...
                                        ' label_column is not provided']);
                                    continue;
                                end
                                reportProgress(report);
                                qfd=sprv.qfDissimilarityTestSetPrior(...
                                    data, inData, testSetLabels, ...
                                    withTraining, false, pu,[],matchStrategy);
                            end
                            if ~isempty(qfd)
                                extras.qfd{end+1}=qfd;
                                if beGraphic
                                    if args.match_table_fig
                                        cascading{end+1}=qfd.fig;
                                    end
                                    if args.match_histogram_fig
                                        if ~isempty(qfd.qHistFig) ...
                                                && matchStrategy==1
                                            cascading{end+1}=qfd.qHistFig;
                                        end
                                        if ~isempty(qfd.fHistFig) && matchStrategy==2
                                            cascading{end+1}=qfd.fHistFig;
                                        end
                                    end
                                end
                            else
                                if creatingPu
                                    pu.close;
                                end
                                return;
                            end
                            
                        end
                    end
                end
                extras.doMatchOutput(nCols-labelCols);
                if beGraphic
                    if args.false_positive_negative_plot
                        ff=extras.seeFalsePosNeg;
                        if ~isempty(ff)
                            cascading{end+1}=ff;
                        end
                    end
                    if ~isempty(cascading)
                        Gui.CascadeFigs(cascading, false, true, 70, 2, ...
                            true, false, scrFig, args.cascade_x);
                    end
                    if ischar(args.csv_file_or_data) && ischar(args.template_file)
                        h1=['<h3>' args.csv_file_or_data '<br>'...
                            args.template_file '</h3>'];
                    elseif ischar(args.template_file)
                        h1=['<h3>' args.template_file '</h3>'];
                    else
                        h1=[];
                    end
                    if ~isempty(args.match_file)
                        extras.saveMatchFiles(h1);
                        if args.match_html==1
                            extras.seeMatches(2, h1);
                        else
                            extras.seeMatches(-1, h1);
                        end
                    else
                        extras.seeMatches(args.match_html, h1);
                    end
                end
                if hasFig
                    figure(scrFig);
                end
            end    
            if creatingPu
                if ~isempty(pu)
                    pu.close;
                end
            end
        end
    end

    function matchProgress(s, ttl)
        if ~beQuiet
            if isempty(args.parent_popUp) && exist('pu', 'var') ...
                    && ~isempty(pu)
                if nargin>1
                    pu.dlg.setTitle(ttl);
                else
                    pu.dlg.setTitle(s);
                end
            end
            fprintf('%s %s\n',args.parent_context,  s);
        end
    end

    
    function keepComputing=progress_report(objectOrString)
        if ischar(objectOrString)
            if ~isequal(objectOrString, ...
                    StochasticGradientDescent.FINDING_ISLANDS) 
                drawnow;
                if ~String.StartsWith(objectOrString, KnnFind.PROGRESS_PREFIX)
                    pu.pb.setValue(pu.pb.getValue+1);
                end
            end
            keepComputing=~pu.cancelled;
            pu.pb.setString(objectOrString);
            pu.pack;
            pu.showTimeSpent;
            return;
        end 
        keepComputing=~pu.cancelled;
        done=objectOrString.getEpochsDone-1;
        toDo=objectOrString.getEpochsToDo;
        pu.pb.setValue(3+(pu.pb.getMaximum*(done/toDo)));
        pu.pb.setMaximum(toDo);
        pu.pb.setString(sprintf('%d/%d epochs done', done, toDo));
        if isvalid(lbl)
            delete(lbl);
        end
        updatePlot(objectOrString.getEmbedding);
        pu.showTimeSpent;
    end

    
    function csvFile=downloadCsv
        csvFile=[];

        zip=fullfile(UmapUtil.LocalSamplesFolder, 'samples.zip');
        if ~isempty(WebDownload.GetZipIfMissing(zip, ...
                WebDownload.ResolveUrl))
            csvFile='sample10k.csv';
            msg(Html.WrapHr(['Samples files are stored in<br>'...
                UmapUtil.LocalSamplesFolder]), 8, 'south east+');
        end
    end


    function reportProgress(report, starting)
        if ~isempty(args.parent_popUp)
            if nargin>1 && starting
                args.parent_popUp.setText(['Test ' ...
                    args.parent_context ' ' args.description]);
                args.parent_popUp.dlg.pack;
            else
                disp([num2str([args.parent_popUp.pb.getValue+1 args.parent_popUp.pb.getMaximum]) ' ' report]);
                args.parent_popUp.incrementProgress;
            end
            args.parent_popUp.setText2(report);
        end
        if ~beQuiet
            fprintf('%s %s\n',args.parent_context, report);
        end
        
    end

    function dispNoDbScan
        warning(['dbscan for clustering in 3+D is not available ... '...
            '\nDownload from MathWorks File Exchange: '...
            'https://www.mathworks.com/matlabcentral/fileexchange/52905-dbscan-clustering-algorithm']);
    end
end
