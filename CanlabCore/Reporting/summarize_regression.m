function S = summarize_regression(OUT, varargin)
% summarize_regression  Generate comprehensive regression analysis report
%
% Creates a complete report for regression analysis results including statistical
% maps, brain region interpretations, diagnostic information, and publication-quality
% figures. The function generates thresholded statistical maps, coverage tables,
% region tables, full-brain montages, and detailed interpretations using hierarchical
% atlas labels.
%
% KEY FEATURES:
% =============
% • Single-file implementation - All helper functions included, no dependencies
% • Professional PDF reports - Uses exportgraphics (no MATLAB Report Generator needed)
% • Multi-image support - Processes both t-statistics and beta coefficients
% • Parallel processing - Optional parallelization across terms for speed
% • Hierarchical interpretation - Generates human-readable summaries using atlas labels
% • HPC compatible - Robust error handling and software OpenGL for cluster nodes
% • Publication-ready figures - High-resolution images with consistent formatting
% • Comprehensive diagnostics - Includes VIF tables and model summaries
%
% INPUTS:
% =======
% OUT - Regression results structure containing:
%   • .t      - statistic_image object with t-statistics (required)
%   • .b      - statistic_image object with beta coefficients (optional, recommended)
%   • .analysis_name - String describing the analysis
%   • .variable_names - Cell array of term names
%   • .diagnostics - Structure containing diagnostic information
%   • .input_parameters - Structure with model parameters
%
% OPTIONAL PARAMETERS (name-value pairs):
% =======================================
% Basic Analysis:
%   'atlas'          - Atlas object for region labeling (default: canlab2024 if available)
%   'p'              - P-value threshold (default: 0.05)
%   'corr'           - Multiple comparison correction: 'fdr', 'unc', 'bfr', 'bonf'
%                     (default: 'fdr')
%   'coverage'       - Minimum coverage percentage for region tables (default: 0)
%   'subdivide'      - Whether to subdivide regions in tables (default: true)
%
% Output Control:
%   'outdir'         - Output directory (default: 'summarize_regression_YYYYMMDD_HHMMSS')
%   'save_figs'      - Save individual figure images (default: true)
%   'make_pdf'       - Generate PDF report (default: true)
%   'pdf_name'       - PDF filename (default: 'summarize_regression_report')
%
% Visualization:
%   'make_montage'   - Create region-centered montages (default: true)
%   'montage_topk'   - Number of top regions to show in montage (default: 40)
%   'montage_labels' - Label regions in montages using atlas (default: true)
%   'make_fullbrain' - Create full-brain HCP-style montages (default: true)
%
% Parallel Processing:
%   'parallel'       - Enable parallel processing across terms (default: false)
%   'n_workers'      - Number of parallel workers (0 = use default pool) (default: 0)
%
% Interpretation Settings:
%   'max_labels1'    - Maximum level 1 labels per level 2 region (default: 8)
%   'max_labels2'    - Maximum level 2 labels per level 3 region (default: 6)
%   'max_labels3'    - Maximum level 3 regions to report (default: 4)
%
% Diagnostics:
%   'vif_flag'       - VIF threshold for flagging collinearity (default: 5)
%
% OUTPUTS:
% ========
% S - Structure containing all analysis results with fields:
%   • .outdir          - Output directory path
%   • .logfile         - Path to text log file
%   • .analysis_name   - Analysis name from input
%   • .variable_names  - Variable names from input
%   • .n_terms         - Number of terms processed
%   • .options         - Options used for analysis
%   • .vif_table       - Table of Variance Inflation Factors
%   • .vif_csv         - Path to CSV file with VIF table
%   • .pdf             - Path to generated PDF report (if make_pdf=true)
%   • .terms{i}        - Cell array of term structures, each containing:
%       - .name              - Term name
%       - .dir               - Term-specific output directory
%       - .n_surviving_voxels- Number of voxels surviving threshold
%       - .threshold_obj     - Thresholded statistic_image object
%       - .coverage_table    - Atlas coverage table
%       - .coverage_csv      - Path to coverage CSV file
%       - .region_table      - Region table with hierarchical labels
%       - .regions_csv       - Path to region CSV file
%       - .rpos/.rneg        - Positive/negative region objects
%       - .flavor_text       - Human-readable interpretation
%       - .fullbrain_t_png   - Full-brain t-stat montage path
%       - .fullbrain_beta_png- Full-brain beta coefficient montage path
%       - .montage_t_png     - Top regions t-stat montage path
%       - .montage_beta_png  - Top regions beta montage (same slices as t-stat)
%       - Error fields (.threshold_error, .coverage_error, etc.) if any step fails
%
% REPORT FORMAT:
% ==============
% The function generates a multi-page PDF report with:
%   PAGE 1: Title Page
%     • Analysis name and generation timestamp
%     • Model summary in plain language
%     • VIF table with collinearity assessment
%
%   PAGES 2-3, 4-5, etc.: Per-term results (2 pages per term)
%     Page A: Four-panel image display (2×2 grid):
%       • Top-left: Full-brain t-statistics montage
%       • Top-right: Full-brain beta coefficients montage
%       • Bottom-left: Top regions t-statistics montage
%       • Bottom-right: Top regions beta coefficients montage
%                     (shown on same slices as t-stat for direct comparison)
%
%     Page B: Interpretation text with:
%       • Color-coded ACTIVATIONS (green) and DEACTIVATIONS (red) sections
%       • Hierarchical region labels (Level 3 → Level 2 → Level 1)
%       • Automatic multi-page continuation for long interpretations
%       • Proper text wrapping and bullet-point formatting
%
% USAGE EXAMPLES:
% ===============
% % Basic usage with default settings
% S = summarize_regression(regression_results);
%
% % Specify output directory and custom threshold
% S = summarize_regression(regression_results, ...
%     'outdir', '/path/to/my_analysis', ...
%     'p', 0.001, 'corr', 'fdr');
%
% % Use custom atlas and enable parallel processing
% S = summarize_regression(regression_results, ...
%     'atlas', my_custom_atlas, ...
%     'parallel', true, 'n_workers', 4);
%
% % Generate minimal output (no PDF, no figures)
% S = summarize_regression(regression_results, ...
%     'make_pdf', false, 'save_figs', false, ...
%     'make_montage', false, 'make_fullbrain', false);
%
% % Detailed interpretation with more labels
% S = summarize_regression(regression_results, ...
%     'max_labels3', 6, 'max_labels2', 8, 'max_labels1', 10);
%
% % HPC cluster usage with software OpenGL
% S = summarize_regression(regression_results, ...
%     'parallel', true, 'n_workers', 12, ...
%     'outdir', '/scratch/user/analysis_001');
%
% NOTES FOR HPC/CLUSTER USE:
% ===========================
% • Parallel figure creation can be finicky on some compute nodes. If you see
%   worker crashes, set 'parallel' to false.
% • The function automatically sets OpenGL to 'software' mode for compatibility
%   with headless nodes and virtual environments.
% • Memory usage scales with number of workers and image size. Monitor memory
%   when processing many terms or large images.
% • Output directories are created automatically with timestamp to prevent
%   overwriting previous results.
%
% TECHNICAL DETAILS:
% ==================
% • Image Processing: Uses CANlab tools (threshold, region, montage methods)
% • PDF Generation: Uses MATLAB's exportgraphics with vector text and raster images
% • Text Wrapping: Custom algorithm handles hierarchical labels and long region lists
% • Error Handling: Each term processed independently; errors in one term don't
%   stop processing of others
% • Beta Montages: Beta coefficient montages use the same anatomical slices as
%   t-stat montages for direct comparison (when t-image regions are available)
%
% Author: Michael Sun, PhD
% Date: 12/28/2025 (Merry Christmas)

% -----------------------
% Parse inputs
% -----------------------
p = inputParser;
p.addRequired('OUT', @(x) isstruct(x) || isa(x,'cell') || isa(x,'statistic_image') || isa(x,'fmri_data') || isa(x,'fmri_data_st'));

p.addParameter('atlas', [], @(x) isempty(x) || isstruct(x) || isobject(x));
p.addParameter('p', .05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('corr', 'fdr', @(x) ischar(x) || isstring(x));
p.addParameter('coverage', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('subdivide', true, @(x) islogical(x) && isscalar(x));
p.addParameter('outdir', '', @(x) ischar(x) || isstring(x));
p.addParameter('save_figs', true, @(x) islogical(x) && isscalar(x));
p.addParameter('make_montage', true, @(x) islogical(x) && isscalar(x));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.addParameter('montage_topk', 42, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('montage_labels', true, @(x) islogical(x) && isscalar(x));
p.addParameter('make_pdf', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('pdf_name', 'summarize_regression_report', @(x)ischar(x)||isstring(x));
p.addParameter('vif_flag', 5, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('write_images', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('write_format', 'nii', @(x)ischar(x)||isstring(x)); % (nii / nii.gz if supported)

% NEW: parallel
p.addParameter('parallel', false, @(x)islogical(x)&&isscalar(x));
p.addParameter('n_workers', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);

% Interpretation controls
p.addParameter('max_labels3', 4, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('max_labels2', 6, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('max_labels1', 8, @(x)isnumeric(x)&&isscalar(x));

% Fullbrain
p.addParameter('make_fullbrain', true, @(x) islogical(x) && isscalar(x));

p.parse(OUT, varargin{:});
opt = p.Results;

% Default atlas if not provided
if isempty(opt.atlas)
    if evalin('base', 'exist(''canlab2024'', ''var'')')
        opt.atlas = evalin('base', 'canlab2024');
    else
        try
            opt.atlas = canlab2024; %#ok<NASGU>
            opt.atlas = canlab2024;
        catch
            warning('summarize_regression:NoAtlas', ...
                'No atlas provided and canlab2024 not found. Atlas-based tables will be skipped.');
            opt.atlas = [];
        end
    end
end

% Output directory
if strlength(string(opt.outdir)) == 0
    stamp = datestr(now, 'yyyymmdd_HHMMSS');
    opt.outdir = fullfile(pwd, ['summarize_regression_' stamp]);
end
opt.outdir = char(string(opt.outdir));
if ~exist(opt.outdir, 'dir'); mkdir(opt.outdir); end

% Headless safety (ok if redundant)
try, opengl('save','software'); catch, end
try, opengl software; catch, end

% Report log (serial write)
logfile = fullfile(opt.outdir, 'report.txt');
fid = fopen(logfile, 'w');
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

% -----------------------
% Extract t images + names
% -----------------------
[timgs, tnames, meta] = local_extract_t_images(OUT);
local_log(fid, opt.verbose, 'Extracted %d t-images\n', numel(timgs));

% number of terms (ALWAYS define this once)
nTerms = numel(timgs);

% Try to extract beta images if available
[betaimgs, betanames] = local_extract_beta_images(OUT, meta);
local_log(fid, opt.verbose, 'Extracted %d beta images\n', numel(betaimgs));

% Check if we have beta images
for i = 1:numel(betaimgs)
    if ~isempty(betaimgs{i})
        local_log(fid, opt.verbose, 'Beta image %d: %s - valid: %d\n', ...
            i, betanames{i}, isa(betaimgs{i}, 'fmri_data'));
    else
        local_log(fid, opt.verbose, 'Beta image %d: %s - EMPTY\n', i, betanames{i});
    end
end

% --- ensure betaimgs aligns with timgs (length + order)
if isempty(betaimgs), betaimgs = cell(nTerms,1); end

if numel(betaimgs) < nTerms
    betaimgs(end+1:nTerms) = {[]};
elseif numel(betaimgs) > nTerms
    betaimgs = betaimgs(1:nTerms);
end


% Summary struct
S = struct();
S.outdir = opt.outdir;
S.logfile = logfile;
S.analysis_name = meta.analysis_name;
S.variable_names = meta.variable_names;
S.n_terms = nTerms;
S.options = opt;

% -----------------------
% Header info (serial)
% -----------------------
local_log(fid, opt.verbose, '=== Summarize Regression ===\n');
local_log(fid, opt.verbose, 'Outdir: %s\n', opt.outdir);
local_log(fid, opt.verbose, 'Analysis: %s\n', meta.analysis_name);

if ~isempty(meta.input_parameters)
    local_log(fid, opt.verbose, '\n--- Input Parameters ---\n');
    local_log(fid, opt.verbose, '%s\n', local_strip_html(evalc('disp(meta.input_parameters)')));

    local_log(fid, opt.verbose, '\n--- Model summary (plain language) ---\n');
    model_summary = local_summarize_input_parameters(meta.input_parameters, opt);
    local_log(fid, opt.verbose, '%s\n', model_summary);
end

Tvif = [];
if ~isempty(meta.diagnostics)
    local_log(fid, opt.verbose, '\n--- Diagnostics ---\n');
    Tvif = local_make_vif_table(meta, opt.vif_flag);
    S.vif_table = Tvif;

    if ~isempty(Tvif) && istable(Tvif) && height(Tvif) > 0
        local_log(fid, opt.verbose, '\n--- VIFs (mapped to variable_names) ---\n');
        local_log(fid, opt.verbose, '%s\n', local_strip_html(local_table_to_text(Tvif)));

        vifPath = fullfile(opt.outdir, 'VIF_table.csv');
        writetable(Tvif, vifPath);
        S.vif_csv = vifPath;
        local_log(fid, opt.verbose, 'Saved VIF table: %s\n', vifPath);

        local_log(fid, opt.verbose, '\n--- VIF interpretation ---\n');
        local_log(fid, opt.verbose, '%s\n', local_vif_flavor(Tvif, opt));
    end
end

if ~isempty(meta.warnings)
    local_log(fid, opt.verbose, '\n--- warnings ---\n');
    local_log_warnings(fid, meta.warnings);
end

local_log(fid, opt.verbose, '\n--- terms detected: %d ---\n', numel(timgs));
for i = 1:numel(timgs)
    local_log(fid, opt.verbose, '  %d) %s\n', i, tnames{i});
end

% -----------------------
% Per-term processing (serial or parallel)
% -----------------------
termStructs = cell(nTerms, 1);
termMsgs    = cell(nTerms, 1);

% Prepare per-term dirs up front (serial)
termDirs = cell(numel(timgs), 1);
termSafeNames = cell(numel(timgs), 1);
for i = 1:numel(timgs)
    termSafeNames{i} = local_safename(tnames{i});
    termDirs{i} = fullfile(opt.outdir, sprintf('%02d_%s', i, termSafeNames{i}));
    if ~exist(termDirs{i}, 'dir'); mkdir(termDirs{i}); end
end

doParallel = opt.parallel;
if doParallel
    % Try to start a pool if none
    try
        pool = gcp('nocreate');
        if isempty(pool)
            if opt.n_workers > 0
                parpool(opt.n_workers);
            else
                parpool;
            end
        elseif opt.n_workers > 0 && pool.NumWorkers ~= opt.n_workers
            % user asked for specific size; best-effort restart
            delete(pool);
            parpool(opt.n_workers);
        end
    catch ME
        warning('summarize_regression:ParallelUnavailable', ...
            'Parallel requested but could not start pool (%s). Falling back to serial.', ME.message);
        doParallel = false;
    end
end

% Update the per-term processing call:
if doParallel
    parfor i = 1:numel(timgs)
        [termStructs{i}, termMsgs{i}] = local_process_one_term( ...
            timgs{i}, tnames{i}, i, opt, opt.atlas, [], false, opt.outdir, betaimgs{i});
    end
else
    for i = 1:numel(timgs)
        [termStructs{i}, termMsgs{i}] = local_process_one_term( ...
            timgs{i}, tnames{i}, i, opt, opt.atlas, fid, opt.verbose, opt.outdir, betaimgs{i});
    end
end

% Write per-term logs (serial)
for i = 1:numel(timgs)
    local_log(fid, opt.verbose, '%s', termMsgs{i});
end

% after loop
S.terms = termStructs;  % or vertcat if you standardized fields


% -----------------------
% PDF
% -----------------------
if opt.make_pdf
    try
        pdfPath = local_write_pdf_report_exportgraphics(S, meta, opt);
        S.pdf = pdfPath;
        local_log(fid, opt.verbose, '\nSaved PDF report: %s\n', pdfPath);
    catch ME
        local_log(fid, opt.verbose, '\nPDF FAILED: %s\n', ME.message);
        S.pdf_error = ME;
    end
end

end % main


% =====================================================================
% Per-term worker (safe for parfor)
% =====================================================================
function [Tstruct, termMsg] = local_process_one_term(timg, tname, i, opt, atlas, fid, verbose, outdir, betaimg)
% Updated to handle beta image

% --- defaults / safety
if nargin < 9 || isempty(betaimg)
    betaimg = [];
end
if nargin < 8 || isempty(outdir)
    outdir = opt.outdir;
end
if nargin < 7 || isempty(verbose)
    verbose = opt.verbose;
end

termName = local_safename(tname);
termDir  = fullfile(outdir, sprintf('%02d_%s', i, termName));
if ~exist(termDir, 'dir'); mkdir(termDir); end

% --- init outputs so they ALWAYS exist
Tstruct = struct();
Tstruct.name = tname;
Tstruct.dir  = termDir;

% Initialize msg using sprintf (more reliable than + operator)
msg = '';
msg = sprintf('%s\n============================\n', msg);
msg = sprintf('%sTERM %d: %s\n', msg, i, string(tname));
msg = sprintf('%s============================\n', msg);

thr = [];
thr_pos = [];
thr_neg = [];

% ---- threshold
try
    % Defensive: ensure single-image
    if isprop(timg,'dat') && ~isempty(timg.dat) && size(timg.dat,2) > 1
        timg = get_wh_image(timg, 1);
    end

    thr = threshold(timg, opt.p, char(opt.corr));

    thr_pos = thr; thr_pos.dat(thr_pos.dat < 0) = 0;
    thr_neg = thr; thr_neg.dat(thr_neg.dat > 0) = 0;

    % Count surviving voxels
    if isa(thr,'statistic_image') && isprop(thr,'sig') && ~isempty(thr.sig)
        nVox = nnz(thr.sig);
    else
        nVox = nnz(thr.dat ~= 0);
    end

    Tstruct.n_surviving_voxels = nVox;
    Tstruct.threshold_obj      = thr;

    if nVox == 0
        msg = sprintf('%sNo voxels survive threshold (p=%.4f, %s) – atlas tables & montage skipped.\n', ...
            msg, opt.p, char(opt.corr));
        Tstruct.coverage_skipped  = true;
        Tstruct.regions_skipped   = true;
        Tstruct.montage_skipped   = true;
        Tstruct.fullbrain_skipped = true;

        termMsg = msg;
        return;
    end

    msg = sprintf('%sThreshold OK: p=%.4f, corr=%s\n', msg, opt.p, char(opt.corr));

catch ME
    msg = sprintf('%sThreshold FAILED: %s\n', msg, ME.message);
    Tstruct.threshold_error = ME;

    termMsg = msg;
    return;
end

% ---- fullbrain (t-statistics)
if opt.make_fullbrain && opt.save_figs
    try
        fbPath_t = local_save_fullbrain_hcp_png(thr, tname, termDir, 't');
        Tstruct.fullbrain_t_png = fbPath_t;
        msg = sprintf('%sSaved fullbrain t-stat montage: %s\n', msg, fbPath_t);
        
        % Also save beta fullbrain if available
        if ~isempty(betaimg)
            try
                fbPath_beta = local_save_fullbrain_hcp_png(betaimg, [tname '_beta'], termDir, 'beta');
                Tstruct.fullbrain_beta_png = fbPath_beta;
                msg = sprintf('%sSaved fullbrain beta montage: %s\n', msg, fbPath_beta);
            catch ME2
                msg = sprintf('%sBeta fullbrain montage FAILED: %s\n', msg, ME2.message);
            end
        end
    catch ME
        msg = sprintf('%sFullbrain montage FAILED: %s\n', msg, ME.message);
        Tstruct.fullbrain_error = ME;
    end
end

% ---- coverage
if ~isempty(atlas)
    try
        Tcov = thr.table_of_atlas_regions_covered('atlas', atlas, 'coverage', 0);

        if istable(Tcov) && any(strcmpi(Tcov.Properties.VariableNames, 'Coverage'))
            Tcov2 = Tcov(Tcov.Coverage > opt.coverage, :);
        else
            Tcov2 = Tcov;
        end

        if istable(Tcov2) && all(ismember({'max_abs_region','Coverage'}, Tcov2.Properties.VariableNames))
            Tcov2 = sortrows(Tcov2, {'max_abs_region','Coverage'}, {'descend','descend'});
        elseif istable(Tcov2) && ismember('Coverage', Tcov2.Properties.VariableNames)
            Tcov2 = sortrows(Tcov2, 'Coverage', 'descend');
        end

        covPath = fullfile(termDir, sprintf('%s_coverage.csv', termName));
        writetable(Tcov2, covPath);

        Tstruct.coverage_table = Tcov2;
        Tstruct.coverage_csv   = covPath;
        msg = sprintf('%sSaved coverage table: %s\n', msg, covPath);

    catch ME
        msg = sprintf('%sCoverage table FAILED: %s\n', msg, ME.message);
        Tstruct.coverage_error = ME;
    end
end

% ---- regions
if ~isempty(atlas)
    try
        rObj = region(thr);

        if opt.subdivide
            [rpos, rneg, rtbl] = rObj.table('subdivide', 'atlas_obj', atlas);
        else
            [rpos, rneg, rtbl] = rObj.table('atlas_obj', atlas);
        end

        regPath = fullfile(termDir, sprintf('%s_regions.csv', termName));
        if istable(rtbl)
            writetable(rtbl, regPath);
        else
            save(fullfile(termDir, sprintf('%s_regions.mat', termName)), 'rtbl', 'rpos', 'rneg');
        end

        Tstruct.rpos = rpos;
        Tstruct.rneg = rneg;
        Tstruct.region_table = rtbl;
        Tstruct.regions_csv  = regPath;

        msg = sprintf('%sSaved region table: %s\n', msg, regPath);

    catch ME
        msg = sprintf('%sRegion table FAILED: %s\n', msg, ME.message);
        Tstruct.regions_error = ME;
    end
end

% ---- montage (t-statistics)
do_montage = opt.make_montage && opt.save_figs && ~isempty(atlas);
if do_montage
    try
        figPath_t = local_save_montage_topk(thr, atlas, tname, termDir, opt.montage_topk, opt.montage_labels, 't');
        Tstruct.montage_t_png = figPath_t;
        msg = sprintf('%sSaved montage (top-%d, t-stat): %s\n', msg, opt.montage_topk, figPath_t);
        
        % Also save beta montage if available - PASS THE THRESHOLDED T-IMAGE
        if ~isempty(betaimg)
            try
                figPath_beta = local_save_montage_topk(betaimg, atlas, [tname '_beta'], termDir, opt.montage_topk, opt.montage_labels, 'beta', thr);
                Tstruct.montage_beta_png = figPath_beta;
                msg = sprintf('%sSaved montage (top-%d, beta): %s\n', msg, opt.montage_topk, figPath_beta);
            catch ME2
                msg = sprintf('%sBeta montage FAILED: %s\n', msg, ME2.message);
            end
        end
    catch ME
        msg = sprintf('%sMontage FAILED: %s\n', msg, ME.message);
        Tstruct.montage_error = ME;
    end
end

% ---- interpretation (single place)
try
    if isempty(atlas)
        flav = "Interpretation skipped (no atlas provided).";
    else
        flav = local_term_flavor_hier(thr_pos, thr_neg, atlas, opt);
    end

    Tstruct.flavor_text = char(flav);
    msg = sprintf('%s\n--- interpretation ---\n%s\n', msg, string(flav));

catch ME
    msg = sprintf('%s\n--- interpretation ---\nFAILED: %s\n', msg, ME.message);
    Tstruct.flavor_error = ME;
end

termMsg = msg;

close all; % close all figure windows.
end

% =====================================================================
% Helpers
% =====================================================================
function figPath = local_save_montage_topk(img, atlas, termTitle, termDir, topk, do_labels, imgType, varargin)
% imgType: 't' for t-statistics, 'beta' for beta coefficients

if nargin < 7 || isempty(imgType)
    imgType = 't';
end

termSafe = local_safename(termTitle);
figPath  = fullfile(termDir, sprintf('%s_montage_%s.png', termSafe, imgType));

persistent did_init
if isempty(did_init)
    try, opengl('save','software'); catch, end
    did_init = true;
end

% For beta images with t-image reference
t_img_for_reference = [];
if nargin > 7 && ~isempty(varargin{1})
    t_img_for_reference = varargin{1};
end

% Create regions differently based on image type
if strcmp(imgType, 't') || isempty(t_img_for_reference)
    % For t-statistics or beta without reference
    rObj = region(img);
    
    if isempty(rObj) || numel(rObj) == 0
        fh = figure('Visible','off','Color','w');
        ax = axes('Parent',fh); axis(ax,'off');
        text(0.5,0.5,'No clusters found', ...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14);
        drawnow;
        rgb = print(fh, '-RGBImage', '-r150');
        imwrite(rgb, figPath);
        close(fh);
        return
    end
    
    % subset to top-k by peak |value|
    rPlot = rObj;
    if numel(rObj) > topk
        peaks = nan(numel(rObj),1);
        for rr = 1:numel(rObj)
            try
                d = rObj(rr).dat;
                if isnumeric(d)
                    peaks(rr) = max(abs(d(:)));
                elseif iscell(d)
                    vals = [];
                    for k = 1:numel(d)
                        if isnumeric(d{k}) && ~isempty(d{k}), vals = [vals; d{k}(:)]; end %#ok<AGROW>
                    end
                    if ~isempty(vals), peaks(rr) = max(abs(vals)); end
                end
            catch
                peaks(rr) = NaN;
            end
        end
        if all(isnan(peaks))
            rPlot = rObj(1:topk);
        else
            [~,ix] = sort(peaks,'descend','MissingPlacement','last');
            rPlot = rObj(ix(1:topk));
        end
    end
    
    % Create montage
    fh = figure('Visible','off','Color','w','Units','pixels','Position',[50 50 1400 900]);
    set(fh,'Renderer','opengl');
    fh_before = fh;
    
    try
        if do_labels && exist('autolabel_regions_using_atlas', 'file')
            autolabel_regions_using_atlas(rPlot, atlas).montage('regioncenters');
        else
            rPlot.montage('regioncenters');
        end
        
        title_str = strrep(termTitle,'_','\_');
        if strcmp(imgType, 'beta')
            title([title_str ' (beta coefficients)'], 'Interpreter','tex', 'FontSize', 12);
        else
            title([title_str ' (t-statistics)'], 'Interpreter','tex', 'FontSize', 12);
        end
        
        drawnow; pause(0.05);
        
        fh_after = gcf;
        rgb = print(fh_after, '-RGBImage', '-r200');
        imwrite(rgb, figPath);
        
    catch ME
        try, close(fh); catch, end
        try, close(gcf); catch, end
        warning('Montage failed: %s', ME.message);
        
        % Create error figure
        fh = figure('Visible','off','Color','w');
        ax = axes('Parent',fh); axis(ax,'off');
        text(0.5,0.5,'Montage generation failed', ...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);
        drawnow;
        rgb = print(fh, '-RGBImage', '-r150');
        imwrite(rgb, figPath);
        close(fh);
        return
    end
    
    try, if isvalid(fh_before), close(fh_before); end, catch, end
    try, if isvalid(fh_after) && fh_after ~= fh_before, close(fh_after); end, catch, end
    
elseif strcmp(imgType, 'beta') && ~isempty(t_img_for_reference)
    % For beta images using t-image slices
    % First, get regions from t-image
    rObj_t = region(t_img_for_reference);
    
    if isempty(rObj_t) || numel(rObj_t) == 0
        % Fallback to regular beta montage
        figPath = local_save_montage_topk(img, atlas, termTitle, termDir, topk, do_labels, 'beta');
        return
    end
    
    % Get top regions from t-image
    rPlot_t = rObj_t;
    if numel(rObj_t) > topk
        peaks = nan(numel(rObj_t),1);
        for rr = 1:numel(rObj_t)
            try
                d = rObj_t(rr).dat;
                if isnumeric(d)
                    peaks(rr) = max(abs(d(:)));
                end
            catch
                peaks(rr) = NaN;
            end
        end
        if all(isnan(peaks))
            rPlot_t = rObj_t(1:topk);
        else
            [~,ix] = sort(peaks,'descend','MissingPlacement','last');
            rPlot_t = rObj_t(ix(1:topk));
        end
    end
    
    % Create montage using t-region centers, but show beta values
    fh = figure('Visible','off','Color','w','Units','pixels','Position',[50 50 1400 900]);
    set(fh,'Renderer','opengl');
    
    try
        % Create display object with region centers from t-image
        if do_labels && exist('autolabel_regions_using_atlas', 'file')
            o2 = autolabel_regions_using_atlas(rPlot_t, atlas).montage('regioncenters');
        else
            o2 = rPlot_t.montage('regioncenters');
        end
        
        % Remove t-blobs and add beta blobs
        o2 = removeblobs(o2);
        o2 = addblobs(o2, img);

        % Redisplay
        % montage(o2);
        
        title_str = strrep(termTitle,'_','\_');
        title([title_str ' (beta coefficients - same slices as t-stat)'], 'Interpreter','tex', 'FontSize', 12);
        
        drawnow; pause(0.05);
        
        % Save figure
        fh = gcf;
        rgb = print(fh, '-RGBImage', '-r200');
        imwrite(rgb, figPath);
        
        % Close figure
        close(fh);
        
    catch ME
        try, close(fh); catch, end
        warning('Beta montage with t-slices failed: %s', ME.message);
        
        % Fallback to regular beta montage
        figPath = local_save_montage_topk(img, atlas, termTitle, termDir, topk, do_labels, 'beta');
    end
end
end

function Tvif = local_make_vif_table(meta, vif_flag)
Tvif = table();
if nargin < 2 || isempty(vif_flag), vif_flag = 5; end
if isempty(meta) || ~isfield(meta,'diagnostics') || isempty(meta.diagnostics), return; end
if ~isfield(meta.diagnostics,'Variance_inflation_factors'), return; end

vif = meta.diagnostics.Variance_inflation_factors(:);
if isempty(vif), return; end

vn = {};
if isfield(meta,'variable_names') && iscell(meta.variable_names)
    vn = meta.variable_names(:);
end

% remove intercept if vn is longer by 1
if ~isempty(vn) && (numel(vn) == numel(vif) + 1)
    intercept_idx = find(contains(lower(string(vn)), "intercept") | contains(lower(string(vn)), "const"), 1, 'last');
    if isempty(intercept_idx), intercept_idx = numel(vn); end
    vn(intercept_idx) = [];
end

% truncate safely
if ~isempty(vn)
    n = min(numel(vn), numel(vif));
    vn = vn(1:n);
    vif = vif(1:n);
else
    vn = arrayfun(@(i) sprintf('X%02d', i), (1:numel(vif))', 'UniformOutput', false);
end

Tvif = table(string(vn), vif, 'VariableNames', {'TERM','VIF'});
Tvif.Flag = strings(height(Tvif), 1);
Tvif.Flag(Tvif.VIF >= vif_flag) = "*";
Tvif = sortrows(Tvif, 'VIF', 'descend');
end


function [timgs, tnames, meta] = local_extract_t_images(OUT)
meta = struct('analysis_name','','variable_names',[],'diagnostics',[],'warnings',[],'input_parameters',[]);
if isstruct(OUT)
    if isfield(OUT,'analysis_name'); meta.analysis_name = OUT.analysis_name; end
    if isfield(OUT,'variable_names'); meta.variable_names = OUT.variable_names; end
    if isfield(OUT,'diagnostics'); meta.diagnostics = OUT.diagnostics; end
    if isfield(OUT,'warnings'); meta.warnings = OUT.warnings; end
    if isfield(OUT,'input_parameters'); meta.input_parameters = OUT.input_parameters; end
end

timgs = {};
tnames = {};

if isstruct(OUT) && isfield(OUT,'t') && ~isempty(OUT.t)
    nimg = [];
    try
        if isprop(OUT.t, 'dat') && ~isempty(OUT.t.dat)
            nimg = size(OUT.t.dat, 2);
        end
    catch
        nimg = [];
    end
    if isempty(nimg)
        try, nimg = numel(OUT.t); catch, nimg = 1; end
    end

    if nimg > 1
        timgs = cell(nimg, 1);
        tnames = cell(nimg, 1);
        for i = 1:nimg
            timgs{i} = get_wh_image(OUT.t, i);
            tnames{i} = sprintf('t_%02d', i);
        end
    else
        timgs = {OUT.t};
        tnames = {'t'};
    end

    if iscell(meta.variable_names) && numel(meta.variable_names) == numel(timgs)
        tnames = meta.variable_names(:);
    end
    return;
end

if isa(OUT,'statistic_image') || isa(OUT,'fmri_data') || isa(OUT,'fmri_data_st')
    nimg = 1;
    try
        if isprop(OUT,'dat') && ~isempty(OUT.dat)
            nimg = size(OUT.dat, 2);
        else
            nimg = numel(OUT);
        end
    catch
        nimg = 1;
    end

    timgs = cell(nimg, 1);
    tnames = cell(nimg, 1);
    for i = 1:nimg
        if nimg > 1
            timgs{i} = get_wh_image(OUT, i);
        else
            timgs{i} = OUT;
        end
        tnames{i} = sprintf('img_%02d', i);
    end

    if strlength(string(meta.analysis_name)) == 0
        meta.analysis_name = 'summarize_regression_input_image';
    end
    return;
end

error('summarize_regression:UnknownInput', ...
    'Could not detect t-images in input. Expected struct with .t or an image_vector/statistic_image input.');
end


function local_log(fid, verbose, varargin)
msg = sprintf(varargin{:});

% Write to file if fid is valid
if exist('fid','var') && ~isempty(fid) && isnumeric(fid) && isscalar(fid) && fid > 1
    fprintf(fid, '%s', msg);
end

% Write to console if verbose
if exist('verbose','var') && ~isempty(verbose) && verbose
    fprintf('%s', msg);
end
end

function local_log_warnings(fid, w)
if isempty(w)
    fprintf(fid, '(none)\n');
    return;
end
if iscell(w)
    for i = 1:numel(w)
        fprintf(fid, ' - %s\n', string(w{i}));
    end
else
    fprintf(fid, ' - %s\n', string(w));
end
end


function s = local_safename(s)
s = char(string(s));
s = regexprep(s, '\s+', '_');
s = regexprep(s, '[^\w\-]', '');
if isempty(s); s = 'term'; end
end


function txt = local_summarize_input_parameters(ip, opt)
if isempty(ip) || ~isstruct(ip)
    txt = 'Model settings were not available.';
    return
end

getf = @(f, d) local_getfield(ip, f, d);

brain_is_predictor = getf('brain_is_predictor', NaN);
do_robust          = getf('do_robust', NaN);
grandmeanscale     = getf('grandmeanscale', NaN);
do_intercept       = getf('do_intercept', NaN);
do_resid           = getf('do_resid', NaN);

thr_str = "";
try
    if isfield(ip,'initial_statistical_threshold') && ~isempty(ip.initial_statistical_threshold)
        t = ip.initial_statistical_threshold;
        if iscell(t) && numel(t) >= 2
            thr_str = sprintf('%s q/p < %.3g', upper(string(t{2})), t{1});
        end
    end
catch
end
if strlength(thr_str)==0
    thr_str = sprintf('%s q/p < %.3g', upper(string(opt.corr)), opt.p);
end

parts = strings(0);

if isequal(brain_is_predictor,1)
    parts(end+1) = "The regression models tested here use brain data to predict each regressor.";
elseif isequal(brain_is_predictor,0)
    parts(end+1) = "The regression models tested here use regressors to predict brain activity.";
end

if isequal(do_robust,1)
    parts(end+1) = "Robust regression was used.";
else
    parts(end+1) = "Ordinary least squares regression was used.";
end

if isequal(grandmeanscale,1)
    parts(end+1) = "Images were grand-mean scaled.";
elseif isequal(grandmeanscale,0)
    parts(end+1) = "Grand means were not scaled between images.";
end

if isequal(do_intercept,1)
    parts(end+1) = "An intercept term was included (or added) in the model.";
else
    parts(end+1) = "No intercept term was included in the model.";
end

if isequal(do_resid,1)
    parts(end+1) = "Residual images were requested for output.";
else
    parts(end+1) = "Residual images were not requested for output.";
end

parts(end+1) = "The initial statistical threshold is " + thr_str + ".";

txt = char(strjoin(parts, " "));
end


function v = local_getfield(s, f, d)
if isfield(s,f) && ~isempty(s.(f))
    v = s.(f);
else
    v = d;
end
end


function s = local_vif_flavor(Tvif, opt)
if isempty(Tvif) || ~istable(Tvif) || ~all(ismember({'TERM','VIF'}, Tvif.Properties.VariableNames))
    s = "VIF diagnostics were not available.";
    return
end

hi = Tvif(Tvif.VIF >= opt.vif_flag, :);

if isempty(hi)
    s = "No VIFs exceed " + opt.vif_flag + ". Collinearity does not appear concerning.";
    return
end

names = string(hi.TERM);
if numel(names) > 8
    names = [names(1:8); "..."];
end

s = "VIFs are elevated (≥ " + opt.vif_flag + ") for: " + strjoin(names, ", ") + ...
    ". Consider centering nuisance covariates, removing redundant regressors, combining highly correlated terms, or re-running with a reduced design to test stability. " + note;
end



% Interpretation (FIXED - better handling of negative activations)
% -----------------------
function txt = local_term_flavor_hier(thr_pos, thr_neg, atlas, opt)
    % Get coverage tables directly (more reliable than checking dat)
    try
        T_pos = thr_pos.table_of_atlas_regions_covered('atlas', atlas, 'coverage', 0);
        pos_has_data = ~isempty(T_pos) && istable(T_pos) && height(T_pos) > 0;
    catch
        T_pos = table();
        pos_has_data = false;
    end
    
    try
        [~, T_neg] = thr_neg.table_of_atlas_regions_covered('atlas', atlas, 'coverage', 0);
        neg_has_data = ~isempty(T_neg) && istable(T_neg) && height(T_neg) > 0;
    catch
        T_neg = table();
        neg_has_data = false;
    end
    
    if pos_has_data
        try
            pos_str = local_hier_summary(T_pos, atlas, opt, "activations");
        catch ME
            pos_str = "";
            warning('Failed to generate positive interpretation: %s', ME.message);
        end
    else
        pos_str = "";
    end
    
    if neg_has_data
        try
            neg_str = local_hier_summary(T_neg, atlas, opt, "deactivations");
        catch ME
            neg_str = "";
            warning('Failed to generate negative interpretation: %s', ME.message);
        end
    else
        neg_str = "";
    end

    if strlength(pos_str)==0 && strlength(neg_str)==0
        txt = "No suprathreshold effects to interpret.";
        return
    end

    if strlength(pos_str)==0
        txt = "Thresholded effects show no major activations. " + neg_str;
    elseif strlength(neg_str)==0
        txt = pos_str + "Thresholded effects show no major deactivations.";
    else
        txt = pos_str + char(10) + char(10) + neg_str;
    end

    txt = char(txt);
end

function s = local_hier_summary(T, atlas, opt, label)
% T should be the coverage table, not a statistic_image object

s = "";

if isempty(T) || ~istable(T) || height(T)==0
    return
end

vn = lower(string(T.Properties.VariableNames));

% Find the correct column names (they might have slightly different names)
i1 = find(contains(vn, ["label", "region"]), 1);
i2 = find(contains(vn, ["label_2", "region_2", "level2"]), 1);
i3 = find(contains(vn, ["label_3", "region_3", "level3"]), 1);

if isempty(i1) || isempty(i2) || isempty(i3)
    % Try to find alternative column names
    i1 = find(~cellfun(@isempty, regexp(vn, '^(labels|regions|label|region)$')), 1);
    i2 = find(~cellfun(@isempty, regexp(vn, '^(labels_2|regions_2|level2|label2|region2)$')), 1);
    i3 = find(~cellfun(@isempty, regexp(vn, '^(labels_3|regions_3|level3|label3|region3)$')), 1);
end

if isempty(i1) || isempty(i2) || isempty(i3)
    % If we still can't find the columns, try a different approach
    s = sprintf("Significant %s detected but cannot generate hierarchical summary (missing column names).", label);
    return
end

col1 = T.Properties.VariableNames{i1};
col2 = T.Properties.VariableNames{i2};
col3 = T.Properties.VariableNames{i3};

% Sort by coverage if present
ic = find(contains(vn, ["coverage", "size", "volume"]), 1);
if ~isempty(ic)
    ccol = T.Properties.VariableNames{ic};
    T = sortrows(T, ccol, 'descend');
else
    % Try to sort by peak statistic if coverage not available
    ic = find(contains(vn, ["peak", "stat", "t", "z"]), 1);
    if ~isempty(ic)
        ccol = T.Properties.VariableNames{ic};
        T = sortrows(T, ccol, 'descend');
    end
end

% Extract unique level 3 labels
try
    L3_values = T.(col3);
    if iscell(L3_values)
        L3 = unique(string(L3_values), 'stable');
    else
        L3 = unique(string(L3_values), 'stable');
    end
    L3 = L3(~ismissing(L3) & strlength(L3) > 0);
    L3 = L3(1:min(opt.max_labels3, numel(L3)));
catch
    s = sprintf("Significant %s detected but error extracting hierarchical labels.", label);
    return
end

if isempty(L3)
    return
end

lines = strings(0,1);
lines(end+1) = sprintf("Thresholded effects show major %s in:", label);

for a = 1:numel(L3)
    lines(end+1) = "• " + L3(a);
    
    % Filter for this level 3 region
    if iscell(T.(col3))
        T3 = T(strcmp(string(T.(col3)), L3(a)), :);
    else
        T3 = T(string(T.(col3)) == L3(a), :);
    end
    
    % Get unique level 2 labels
    try
        if iscell(T3.(col2))
            L2 = unique(string(T3.(col2)), 'stable');
        else
            L2 = unique(string(T3.(col2)), 'stable');
        end
        L2 = L2(~ismissing(L2) & strlength(L2) > 0);
        L2 = L2(1:min(opt.max_labels2, numel(L2)));
    catch
        continue
    end
    
    for b = 1:numel(L2)
        % Filter for this level 2 region
        if iscell(T3.(col2))
            T32 = T3(strcmp(string(T3.(col2)), L2(b)), :);
        else
            T32 = T3(string(T3.(col2)) == L2(b), :);
        end
        
        % Get unique level 1 labels
        try
            if iscell(T32.(col1))
                L1 = unique(string(T32.(col1)), 'stable');
            else
                L1 = unique(string(T32.(col1)), 'stable');
            end
            L1 = L1(~ismissing(L1) & strlength(L1) > 0);
            L1 = L1(1:min(opt.max_labels1, numel(L1)));
        catch
            continue
        end
        
        if ~isempty(L1)
            lines(end+1) = "    – " + L2(b) + ": " + strjoin(L1, ", ");
        end
    end
end

s = strjoin(lines, newline);
end


% -----------------------
% Figures (robust headless)
% -----------------------
function figPath = local_save_fullbrain_hcp_png(img, termTitle, termDir, imgType)
% Save full brain montage

if nargin < 4 || isempty(imgType)
    imgType = 't';
end

termSafe = local_safename(termTitle);
figPath  = fullfile(termDir, sprintf('%s_fullbrain_hcp_%s.png', termSafe, imgType));

% This is often required on compute nodes
try, opengl('save','software'); catch, end
try, opengl software; catch, end

% Track figures so we can grab the one montage actually created
figs_before = findall(0,'Type','figure');

% Make montage (may create/activate its own figure)
img.montage('full hcp', 'coordinates', true);

drawnow; pause(0.5); drawnow;

figs_after = findall(0,'Type','figure');
new_figs = setdiff(figs_after, figs_before);

if ~isempty(new_figs)
    f = new_figs(1);
else
    f = gcf;
end

try
    set(f, 'Color','w');
    set(f, 'Renderer','opengl');
catch
end

% Add title if possible
try
    ax = get(f,'CurrentAxes');
    if ~isempty(ax) && isvalid(ax)
        title_str = strrep(termTitle,'_','\_');
        if strcmp(imgType, 'beta')
            title(ax, [title_str ' (beta, unthresholded)'], 'Interpreter','tex', 'FontSize', 12);
        else
            title(ax, title_str, 'Interpreter','tex', 'FontSize', 12);
        end
    end
catch
end

drawnow; pause(0.5); drawnow;

% IMPORTANT: use getframe for OpenGL surfaces
try
    fr = getframe(f);
    imwrite(fr.cdata, figPath);
catch
    % fallback
    rgb = print(f, '-RGBImage', '-r200');
    imwrite(rgb, figPath);
end

% Cleanup: close only the montage figure we used (don't nuke all)
try, close(f); catch, end
end

% -----------------------
% PDF (direct PDF with proper text)
% -----------------------
function pdfPath = local_write_pdf_report_exportgraphics(S, meta, opt)
    outdir = char(S.outdir);
    pdfPath = fullfile(outdir, sprintf('%s.pdf', char(string(opt.pdf_name))));
    if exist(pdfPath,'file'); delete(pdfPath); end
    
    % Create a figure for PDF pages - use tight layout
    fh = figure('Visible','off','Color','w','Units','inches','Position',[0.25 0.25 8 10.5]); % Letter size with small margins
    set(fh, 'Renderer', 'painters'); % Vector graphics for text
    
    % ---------- PAGE 1: TITLE PAGE ----------
    clf(fh);
    
    % Title at top
    ax_title = axes('Position',[0.05 0.88 0.9 0.08], 'Visible','off');
    text(ax_title, 0.5, 0.5, local_strip_html(S.analysis_name), ...
        'FontSize',18,'FontWeight','bold','HorizontalAlignment','center','Interpreter','none');
    
    % Date
    axes('Position',[0.05 0.82 0.9 0.04], 'Visible','off');
    text(0.5, 0.5, ['Generated: ' datestr(now, 'dd-mmm-yyyy HH:MM:SS')], ...
        'FontSize',9,'HorizontalAlignment','center','Interpreter','none');
    
    % Summary box
    ax_summary = axes('Position',[0.05 0.65 0.9 0.15], 'Visible','off');
    text(ax_summary, 0, 1, 'ANALYSIS SUMMARY', 'FontSize',12,'FontWeight','bold','VerticalAlignment','top','Interpreter','none');
    
    model_summary_string = "";
    if ~isempty(meta.input_parameters)
        model_summary_string = local_summarize_input_parameters(meta.input_parameters, opt);
    end
    
    summary_lines = local_wrap_text(char(model_summary_string), 90);
    text(ax_summary, 0, 0.8, summary_lines, 'FontSize',9,'VerticalAlignment','top','Interpreter','none');
    
    % VIF Table if available
    Tvif = table();
    if isfield(S,'vif_table') && ~isempty(S.vif_table), Tvif = S.vif_table; end
    
    if ~isempty(Tvif) && height(Tvif) > 0
        ax_vif = axes('Position',[0.05 0.35 0.9 0.25], 'Visible','off');
        text(ax_vif, 0, 1, 'VARIANCE INFLATION FACTORS (VIF)', 'FontSize',11,'FontWeight','bold','VerticalAlignment','top','Interpreter','none');
        
        % Create formatted table
        y_pos = 0.9;
        text(ax_vif, 0, y_pos, 'Term', 'FontSize',10,'FontWeight','bold','VerticalAlignment','top');
        text(ax_vif, 0.35, y_pos, 'VIF', 'FontSize',10,'FontWeight','bold','VerticalAlignment','top');
        text(ax_vif, 0.5, y_pos, 'Note', 'FontSize',10,'FontWeight','bold','VerticalAlignment','top');
        
        y_pos = y_pos - 0.04;
        line_height = 0.035;
        
        for row = 1:height(Tvif)
            text(ax_vif, 0, y_pos, string(Tvif.TERM(row)), 'FontSize',9,'VerticalAlignment','top','Interpreter','none');
            text(ax_vif, 0.35, y_pos, sprintf('%.2f', Tvif.VIF(row)), 'FontSize',9,'VerticalAlignment','top');
            
            if Tvif.VIF(row) >= opt.vif_flag
                text(ax_vif, 0.5, y_pos, 'High collinearity', 'FontSize',9,'FontWeight','bold','Color',[0.8 0 0],'VerticalAlignment','top');
            else
                text(ax_vif, 0.5, y_pos, 'OK', 'FontSize',9,'Color',[0 0.5 0],'VerticalAlignment','top');
            end
            
            y_pos = y_pos - line_height;
        end
    end
    
    % Footer
    ax_footer = axes('Position',[0.05 0.02 0.9 0.03], 'Visible','off');
    text(ax_footer, 0.5, 0.5, 'Regression Analysis Report • summarize_regression • CANlab Tools', ...
        'FontSize',8,'FontAngle','italic','HorizontalAlignment','center','Interpreter','none');
    
    % Export title page
    try
        exportgraphics(fh, pdfPath, 'ContentType', 'vector', 'Resolution', 300);
    catch ME
        warning('PDF export failed on title page: %s', ME.message);
        close(fh);
        rethrow(ME);
    end
    
    % ---------- TERM PAGES ----------
    for i = 1:numel(S.terms)
        term = S.terms{i};
        term_name = string(term.name);
        
        % PAGE A: IMAGES (2x2 grid) - Maximize space
        clf(fh);
        
        % Term header - minimal space
        ax_header = axes('Position',[0.05 0.96 0.9 0.03], 'Visible','off');
        text(ax_header, 0.5, 0.5, sprintf('TERM %d: %s', i, term_name), ...
            'FontSize',12,'FontWeight','bold','HorizontalAlignment','center','Interpreter','none');
        
        % Check which images are available
        has_t_fullbrain = isfield(term,'fullbrain_t_png') && ~isempty(term.fullbrain_t_png) && exist(term.fullbrain_t_png, 'file');
        has_beta_fullbrain = isfield(term,'fullbrain_beta_png') && ~isempty(term.fullbrain_beta_png) && exist(term.fullbrain_beta_png, 'file');
        has_t_montage = isfield(term,'montage_t_png') && ~isempty(term.montage_t_png) && exist(term.montage_t_png, 'file');
        has_beta_montage = isfield(term,'montage_beta_png') && ~isempty(term.montage_beta_png) && exist(term.montage_beta_png, 'file');
        
        % Create 2x2 grid of images - maximize space
        img_width = 0.44;
        img_height = 0.44;
        
        % Row 1, Col 1: t-stat fullbrain
        ax1 = axes('Position',[0.03 0.48 img_width img_height]);
        if has_t_fullbrain
            try
                img = imread(term.fullbrain_t_png);
                imshow(img);
                title('Full Brain: t-statistics', 'FontSize',8, 'FontWeight', 'normal');
            catch
                text(ax1, 0.5, 0.5, 'Image load error', 'HorizontalAlignment','center', 'FontSize',7, 'Color', [0.5 0.5 0.5]);
            end
        else
            text(ax1, 0.5, 0.5, 'No t-stat fullbrain', 'HorizontalAlignment','center', 'FontSize',7, 'Color', [0.5 0.5 0.5]);
        end
        axis(ax1, 'off');
        
        % Row 1, Col 2: beta fullbrain
        ax2 = axes('Position',[0.53 0.48 img_width img_height]);
        if has_beta_fullbrain
            try
                img = imread(term.fullbrain_beta_png);
                imshow(img);
                title('Full Brain: Beta coefficients', 'FontSize',8, 'FontWeight', 'normal');
            catch
                text(ax2, 0.5, 0.5, 'Image load error', 'HorizontalAlignment','center', 'FontSize',7, 'Color', [0.5 0.5 0.5]);
            end
        else
            text(ax2, 0.5, 0.5, 'No beta fullbrain', 'HorizontalAlignment','center', 'FontSize',7, 'Color', [0.5 0.5 0.5]);
        end
        axis(ax2, 'off');
        
        % Row 2, Col 1: t-stat montage
        ax3 = axes('Position',[0.03 0.02 img_width img_height]);
        if has_t_montage
            try
                img = imread(term.montage_t_png);
                imshow(img);
                title(sprintf('Top %d Regions: t-statistics', opt.montage_topk), 'FontSize',8, 'FontWeight', 'normal');
            catch
                text(ax3, 0.5, 0.5, 'Image load error', 'HorizontalAlignment','center', 'FontSize',7, 'Color', [0.5 0.5 0.5]);
            end
        else
            text(ax3, 0.5, 0.5, 'No t-stat montage', 'HorizontalAlignment','center', 'FontSize',7, 'Color', [0.5 0.5 0.5]);
        end
        axis(ax3, 'off');
        
        % Row 2, Col 2: beta montage
        ax4 = axes('Position',[0.53 0.02 img_width img_height]);
        if has_beta_montage
            try
                img = imread(term.montage_beta_png);
                imshow(img);
                title(sprintf('Top %d Regions: Beta coefficients', opt.montage_topk), 'FontSize',8, 'FontWeight', 'normal');
            catch
                text(ax4, 0.5, 0.5, 'Image load error', 'HorizontalAlignment','center', 'FontSize',7, 'Color', [0.5 0.5 0.5]);
            end
        else
            text(ax4, 0.5, 0.5, 'No beta montage', 'HorizontalAlignment','center', 'FontSize',7, 'Color', [0.5 0.5 0.5]);
        end
        axis(ax4, 'off');
        
        % Page number
        ax_pagenum = axes('Position',[0.05 0.01 0.9 0.01], 'Visible','off');
        text(ax_pagenum, 0.5, 0.5, sprintf('Page %d', 1 + (i-1)*2 + 1), ...
            'FontSize',7,'HorizontalAlignment','center','Interpreter','none', 'Color', [0.3 0.3 0.3]);
        
        % Append images page
        try
            exportgraphics(fh, pdfPath, 'ContentType', 'vector', 'Resolution', 300, 'Append', true);
        catch ME
            warning('PDF export failed on images page for term %d: %s', i, ME.message);
            continue;
        end
        
        % PAGE B: INTERPRETATION - Handle multi-page if needed
        flav = "No interpretation available.";
        if isfield(term,'flavor_text') && ~isempty(term.flavor_text)
            flav = string(term.flavor_text);
        end
        
        % Remove any page numbering like "(1/1)" from the text
        flav = regexprep(flav, '\(\d+/\d+\)', '');
        
        % Split into logical sections for better pagination
        flav_sections = split_interpretation_into_sections(flav);
        
        % Create interpretation pages (as many as needed)
        for section_num = 1:numel(flav_sections)
            clf(fh);
            
            % Term header
            ax_header = axes('Position',[0.05 0.96 0.9 0.03], 'Visible','off');
            if section_num == 1
                header_text = sprintf('TERM %d: %s - Interpretation', i, term_name);
            else
                header_text = sprintf('TERM %d: %s - Interpretation (cont.)', i, term_name);
            end
            text(ax_header, 0.5, 0.5, header_text, ...
                'FontSize',12,'FontWeight','bold','HorizontalAlignment','center','Interpreter','none');
            
            % Text area - maximize space
            ax_text = axes('Position',[0.05 0.05 0.9 0.89], 'Visible','off');
            
            % Display this section
            section_text = flav_sections{section_num};
            display_interpretation_section(ax_text, section_text, section_num == 1);
            
            % Page number
            ax_pagenum = axes('Position',[0.05 0.01 0.9 0.01], 'Visible','off');
            page_num = 1 + (i-1)*2 + 1 + section_num;
            text(ax_pagenum, 0.5, 0.5, sprintf('Page %d', page_num), ...
                'FontSize',7,'HorizontalAlignment','center','Interpreter','none', 'Color', [0.3 0.3 0.3]);
            
            % Append interpretation page
            try
                exportgraphics(fh, pdfPath, 'ContentType', 'vector', 'Resolution', 300, 'Append', true);
            catch ME
                warning('PDF export failed on interpretation page %d for term %d: %s', section_num, i, ME.message);
                break;
            end
        end
    end
    
    close(fh);
end

function sections = split_interpretation_into_sections(interpretation_text)
    % Split interpretation text into manageable sections for pagination
    
    sections = {};
    if isempty(interpretation_text)
        sections{1} = 'No interpretation available.';
        return;
    end
    
    % Convert to char for string operations
    text_str = char(interpretation_text);
    
    % Find activation and deactivation sections
    activation_start = strfind(text_str, 'Thresholded effects show major activations in:');
    deactivation_start = strfind(text_str, 'Thresholded effects show major deactivations in:');
    
    if ~isempty(activation_start) && ~isempty(deactivation_start) && activation_start < deactivation_start
        % Both sections exist and are in order
        activation_section = text_str(activation_start:deactivation_start-1);
        deactivation_section = text_str(deactivation_start:end);
        
        % Split each section if too long
        activation_parts = split_text_into_pages(activation_section, 1500);
        deactivation_parts = split_text_into_pages(deactivation_section, 1500);
        
        sections = [activation_parts, deactivation_parts];
    elseif ~isempty(activation_start)
        % Only activations
        sections = split_text_into_pages(text_str, 1500);
    elseif ~isempty(deactivation_start)
        % Only deactivations
        sections = split_text_into_pages(text_str, 1500);
    else
        % No clear sections
        sections = split_text_into_pages(text_str, 1500);
    end
    
    % Ensure no section is empty
    sections = sections(~cellfun(@isempty, sections));
    if isempty(sections)
        sections{1} = 'Interpretation text unavailable.';
    end
end

function page_parts = split_text_into_pages(text_str, max_chars_per_page)
    % Split text into pages of approximately max_chars_per_page length
    
    page_parts = {};
    
    if isempty(text_str)
        page_parts{1} = '';
        return;
    end
    
    text_length = length(text_str);
    if text_length <= max_chars_per_page
        page_parts{1} = text_str;
        return;
    end
    
    % Split by paragraphs (double newlines)
    paragraphs = split(text_str, newline + newline);
    
    current_page = '';
    for p = 1:numel(paragraphs)
        paragraph = strtrim(paragraphs{p});
        if isempty(paragraph)
            continue;
        end
        
        % Add newline if not first paragraph
        if ~isempty(current_page)
            paragraph = [newline newline paragraph];
        end
        
        % Check if adding this paragraph would exceed page limit
        if length(current_page) + length(paragraph) <= max_chars_per_page
            current_page = [current_page paragraph];
        else
            % Start new page
            if ~isempty(current_page)
                page_parts{end+1} = current_page;
            end
            current_page = paragraph;
        end
    end
    
    % Add the last page
    if ~isempty(current_page)
        page_parts{end+1} = current_page;
    end
end

function display_interpretation_section(ax, section_text, is_first_section)
    % Display interpretation text with proper formatting
    
    axes(ax);
    set(ax, 'Visible', 'off');
    
    % Parse the section text
    lines = splitlines(section_text);
    
    % Starting position (top of available space)
    y_pos = 0.98;
    line_height = 0.025;
    bullet_indent = 0.02;
    subbullet_indent = 0.04;
    
    for line_idx = 1:numel(lines)
        line_text = strtrim(char(lines(line_idx)));
        
        if isempty(line_text)
            y_pos = y_pos - line_height/2; % Small gap for empty lines
            continue;
        end
        
        % Determine line type and format accordingly
        if contains(line_text, 'Thresholded effects show major activations in:')
            % Main header
            text(ax, 0, y_pos, 'ACTIVATIONS', 'FontSize',11, 'FontWeight','bold', ...
                'VerticalAlignment','top', 'Interpreter','none', 'Color', [0 0.5 0]);
            y_pos = y_pos - line_height * 1.5;
            
        elseif contains(line_text, 'Thresholded effects show major deactivations in:')
            % Main header
            text(ax, 0, y_pos, 'DEACTIVATIONS', 'FontSize',11, 'FontWeight','bold', ...
                'VerticalAlignment','top', 'Interpreter','none', 'Color', [0.8 0 0]);
            y_pos = y_pos - line_height * 1.5;
            
        elseif startsWith(line_text, '•')
            % Main bullet point (region)
            bullet_text = extractAfter(line_text, '•');
            text(ax, bullet_indent, y_pos, ['• ' strtrim(bullet_text)], 'FontSize',10, 'FontWeight','bold', ...
                'VerticalAlignment','top', 'Interpreter','none');
            y_pos = y_pos - line_height;
            
        elseif startsWith(line_text, '–') || startsWith(line_text, '-')
            % Sub-bullet point (subregion)
            if startsWith(line_text, '–')
                subbullet_text = extractAfter(line_text, '–');
            else
                subbullet_text = extractAfter(line_text, '-');
            end
            text(ax, subbullet_indent, y_pos, ['  - ' strtrim(subbullet_text)], 'FontSize',9, 'FontWeight','normal', ...
                'VerticalAlignment','top', 'Interpreter','none');
            y_pos = y_pos - line_height;
            
        else
            % Regular text - wrap if too long
            wrapped_lines = splitlines(local_wrap_text(line_text, 100));
            for w = 1:numel(wrapped_lines)
                wrapped_line = strtrim(char(wrapped_lines(w)));
                if ~isempty(wrapped_line)
                    text(ax, 0, y_pos, wrapped_line, 'FontSize',9, ...
                        'VerticalAlignment','top', 'Interpreter','none');
                    y_pos = y_pos - line_height;
                end
            end
        end
        
        % Check if we're running out of space
        if y_pos < 0.02
            % Add continuation note
            text(ax, 0.5, 0.01, '... continues on next page ...', 'FontSize',8, 'FontAngle','italic', ...
                'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'Interpreter','none', 'Color', [0.5 0.5 0.5]);
            break;
        end
    end
end

function wrapped = local_wrap_text(text_str, max_chars)
    % Improved text wrapping that handles long strings of regions
    if isempty(text_str)
        wrapped = '';
        return;
    end
    
    % Replace multiple spaces with single space
    text_str = regexprep(text_str, '\s+', ' ');
    text_str = strtrim(text_str);
    
    % Split by commas for region lists
    if contains(text_str, ',') && length(text_str) > max_chars/2
        % Try to wrap at commas for region lists
        parts = split(text_str, ',');
        lines = {};
        current_line = '';
        
        for i = 1:length(parts)
            part = strtrim(parts{i});
            if i < length(parts)
                part = [part ',']; % Add back comma
            end
            
            if isempty(current_line)
                current_line = part;
            elseif length(current_line) + length(part) + 1 <= max_chars
                current_line = [current_line ' ' part];
            else
                lines{end+1} = current_line;
                current_line = part;
            end
        end
        
        if ~isempty(current_line)
            lines{end+1} = current_line;
        end
        
        wrapped = strjoin(lines, newline);
    else
        % Regular word wrapping
        words = strsplit(text_str, ' ');
        lines = {};
        current_line = '';
        
        for i = 1:length(words)
            word = words{i};
            
            if isempty(current_line)
                current_line = word;
            elseif length(current_line) + length(word) + 1 <= max_chars
                current_line = [current_line ' ' word];
            else
                lines{end+1} = current_line;
                current_line = word;
            end
        end
        
        if ~isempty(current_line)
            lines{end+1} = current_line;
        end
        
        wrapped = strjoin(lines, newline);
    end
end

% -----------------------
% Utilities
% -----------------------
function s = local_strip_html(s)
s = char(string(s));
s = regexprep(s, '<[^>]*>', '');   % remove <strong> etc.
end


function txt = local_table_to_text(T)
% consistent, plain-text table rendering
try
    txt = evalc('disp(T)');
catch
    % fallback minimal
    txt = sprintf('[table: %d x %d]\n', height(T), width(T));
end
txt = local_strip_html(txt);
end

function local_render_term_page(pagePath, termIdx, termName, fullbrainPath, montagePath, flavorText)

fh = figure('Visible','off','Color','w','Units','pixels','Position',[0 0 1654 2339]); % A4 @200dpi
ax = axes('Position',[0 0 1 1],'Visible','off');

text(0.05,0.96, sprintf('TERM %d: %s', termIdx, char(termName)), ...
    'FontSize',22,'FontWeight','bold','Interpreter','none');

% Two image panels
% Left: fullbrain
if strlength(fullbrainPath) > 0 && exist(fullbrainPath,'file')
    I = imread(fullbrainPath);
    axes('Position',[0.06 0.48 0.42 0.42]); imshow(I); axis off
else
    text(0.08,0.68,'(No fullbrain image)','FontSize',12,'Interpreter','none');
end

% Right: montage
if strlength(montagePath) > 0 && exist(montagePath,'file')
    I = imread(montagePath);
    axes('Position',[0.52 0.48 0.42 0.42]); imshow(I); axis off
else
    text(0.74,0.68,'(No montage image)','FontSize',12,'Interpreter','none');
end

% Interpretation block
axes('Position',[0 0 1 1],'Visible','off');
text(0.05,0.42,'Interpretation','FontSize',16,'FontWeight','bold','Interpreter','none');

flav = local_wrap_mono(char(flavorText), 120);
text(0.05,0.40, flav, 'FontName','Courier', 'FontSize',11, ...
    'Interpreter','none', 'VerticalAlignment','top');

print(fh, pagePath, '-dpng','-r200');
close(fh);
end

function out = local_wrap_mono(s, width)
if nargin < 2, width = 100; end
s = char(string(s));
s = regexprep(s, '\r\n|\r', '\n');
lines = splitlines(string(s));

wrapped = strings(0,1);
for i = 1:numel(lines)
    L = char(lines(i));
    while strlength(string(L)) > width
        cut = width;
        % try to break at space
        sp = find(L(1:width)==' ', 1, 'last');
        if ~isempty(sp), cut = sp; end
        wrapped(end+1,1) = string(strtrim(L(1:cut))); %#ok<AGROW>
        L = strtrim(L(cut+1:end));
    end
    wrapped(end+1,1) = string(L); %#ok<AGROW>
end

out = char(strjoin(wrapped, newline));
end

function [betaimgs, betanames] = local_extract_beta_images(OUT, meta)
% Extract beta (unthresholded) images from .b field

betaimgs = {};
betanames = {};

% Default to empty cells
n_terms = numel(meta.variable_names);
betaimgs = cell(n_terms, 1);
betanames = meta.variable_names(:);

if isstruct(OUT) && isfield(OUT,'b') && ~isempty(OUT.b)
    try
        b_img = OUT.b;
        
        % Determine number of images
        nimg = 1;
        if isprop(b_img, 'dat') && ~isempty(b_img.dat)
            nimg = size(b_img.dat, 2);
        else
            nimg = numel(b_img);
        end
        
        if nimg > 1
            betaimgs = cell(nimg, 1);
            betanames = cell(nimg, 1);
            for i = 1:nimg
                try
                    betaimgs{i} = get_wh_image(b_img, i);
                    betanames{i} = sprintf('beta_%02d', i);
                catch
                    betaimgs{i} = [];
                    betanames{i} = sprintf('beta_%02d', i);
                end
            end
        else
            betaimgs = {b_img};
            betanames = {'beta'};
        end
        
        % Use same names as t-images if available
        if iscell(meta.variable_names) && numel(meta.variable_names) == numel(betaimgs)
            betanames = meta.variable_names(:);
        end
        
        local_log([], true, 'Successfully extracted %d beta images\n', numel(betaimgs));
        
    catch ME
        local_log([], true, 'Error extracting beta images: %s\n', ME.message);
        % Return empty cells
        if n_terms > 0
            betaimgs = cell(n_terms, 1);
            betanames = meta.variable_names(:);
        end
    end
else
    local_log([], true, 'No .b field found in input structure\n');
end
end


function local_write_image(img, outpath)
% Robust writer for canlab objects (statistic_image / fmri_data / image_vector)
if isempty(img), error('Image is empty'); end

% Ensure single image
try
    if isprop(img,'dat') && size(img.dat,2) > 1
        img = get_wh_image(img, 1);
    end
catch
end

% Prefer object write methods if available
if ismethod(img,'write')
    img.write(outpath);
    return
end

% Some canlab classes implement write(img, filename)
try
    write(img, outpath);
    return
catch
end

error('No supported write() method found for class %s', class(img));
end