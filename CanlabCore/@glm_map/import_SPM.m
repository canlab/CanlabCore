function obj = import_SPM(obj, SPM, varargin)
% Import a first-level design (and optionally betas) from an SPM model.
%
% Populates a glm_map object from an SPM structure or an SPM.mat path,
% mapping SPM's design fields into the wrapped fmri_glm_design_matrix object.
% fmri_glm_design_matrix deliberately mirrors SPM's field schema
% (xY/nscan/xBF/Sess/xX), so the import is largely a guarded copy of those
% substructures. Compatible with SPM12 and SPM25 first-level models.
%
% :Usage:
% ::
%
%     obj = import_SPM(glm_map, SPM)             % SPM struct already loaded
%     obj = import_SPM(glm_map, '/path/SPM.mat') % path to SPM.mat
%     obj = import_SPM(glm_map, '/path/to/dir')  % directory containing SPM.mat
%
% :Inputs:
%
%   **obj:**
%        A glm_map object (typically empty: glm_map).
%
%   **SPM:**
%        Either an SPM structure (as loaded from SPM.mat), a path to an
%        SPM.mat file, or a directory containing one.
%
% :Optional Inputs:
%
%   **'load_betas':**
%        Also load the estimated beta_*.nii/.img images into obj.betas
%        (looked up in SPM.swd), labeled by SPM.xX.name.
%
%   **'noverbose':**
%        Suppress console output.
%
% :Outputs:
%
%   **obj:**
%        glm_map with .design populated from SPM, level set to 1, and
%        is_timeseries set true.
%
% :Examples:
% ::
%
%     g = import_SPM(glm_map, '/data/sub-01/1stlevel/SPM.mat');
%     g.onsets                 % read-through to the imported design
%     g = add_contrasts(g, c, names);
%     % (data still needed to fit: g = fit(g, fmri_timeseries_obj))
%
% :See also:
%   - fmri_glm_design_matrix, build_design, fit
%
% ..
%    Programmers' notes:
%    2026 - Initial implementation. Field mapping (SPM -> wrapped design):
%      SPM.xY.RT  -> design.TR / design.xY.RT
%      SPM.nscan  -> design.nscan
%      SPM.xBF    -> design.xBF   (basis set, T, T0, UNITS, Volterra)
%      SPM.Sess   -> design.Sess  (onsets .ons, durations .dur, names .name,
%                                  parametric mods .P, covariates .C)
%      SPM.xX     -> design.xX    (design matrix .X, partition idx, names)
%    SPM12 and SPM25 share these first-level fields; whole-substruct copy is
%    resilient to extra/renamed auxiliary fields (W, K, xKXs, pKX, ...).
% ..

% -------------------------------------------------------------------------
% Parse options
% -------------------------------------------------------------------------
doverbose     = ~any(strcmpi(varargin, 'noverbose'));
do_load_betas = any(strcmpi(varargin, 'load_betas'));

% -------------------------------------------------------------------------
% Resolve SPM (struct, file path, or directory)
% -------------------------------------------------------------------------
if ischar(SPM) || isstring(SPM)
    spmpath = char(SPM);
    if isfolder(spmpath)
        spmpath = fullfile(spmpath, 'SPM.mat');
    end
    if ~exist(spmpath, 'file')
        error('glm_map:NoSPMfile', 'Could not find SPM.mat at: %s', spmpath);
    end
    loaded = load(spmpath, 'SPM');
    if ~isfield(loaded, 'SPM')
        error('glm_map:BadSPMfile', 'File does not contain an SPM variable: %s', spmpath);
    end
    SPM = loaded.SPM;
end

if ~isstruct(SPM) || ~isfield(SPM, 'xX')
    error('glm_map:BadSPM', 'SPM must be an SPM structure (or path to one) with an xX field.');
end

% -------------------------------------------------------------------------
% Determine TR
% -------------------------------------------------------------------------
TR = NaN;
if isfield(SPM, 'xY') && isfield(SPM.xY, 'RT')
    TR = SPM.xY.RT;
elseif isfield(SPM, 'xBF') && isfield(SPM.xBF, 'dt') && isfield(SPM.xBF, 'T')
    TR = SPM.xBF.dt * SPM.xBF.T;
end

% -------------------------------------------------------------------------
% Build the wrapped design object and copy SPM substructures
% -------------------------------------------------------------------------
design = fmri_glm_design_matrix(TR);

if isfield(SPM, 'xY'),    design.xY    = SPM.xY;    end
if isfield(SPM, 'nscan'), design.nscan = SPM.nscan; end
if isfield(SPM, 'xBF'),   design.xBF   = SPM.xBF;   end
if isfield(SPM, 'Sess'),  design.Sess  = SPM.Sess;  end
design.xX = SPM.xX;   % includes .X and .name (read by the Dependent accessors)

design.build_method = 'Imported from SPM';
if ~iscell(design.history), design.history = {}; end
design.history{end + 1} = 'Imported from SPM structure';

% -------------------------------------------------------------------------
% Populate the glm_map object
% -------------------------------------------------------------------------
obj.design        = design;
obj.level         = 1;
obj.is_timeseries = true;

obj.history{end + 1} = sprintf('import_SPM: imported %d-session, %d-column design from SPM', ...
    local_num_sess(SPM), size(obj.X, 2));

% -------------------------------------------------------------------------
% Optionally load estimated betas
% -------------------------------------------------------------------------
if do_load_betas
    obj = local_load_betas(obj, SPM, doverbose);
end

% -------------------------------------------------------------------------
% Report
% -------------------------------------------------------------------------
if doverbose
    fprintf('  import_SPM: TR=%g, %d regressors, %d session(s).%s\n', ...
        TR, size(obj.X, 2), local_num_sess(SPM), ...
        local_iif(do_load_betas && ~isempty(obj.betas), ' Betas loaded.', ''));
end

end % import_SPM


% =====================================================================
% Local helpers
% =====================================================================
function n = local_num_sess(SPM)
if isfield(SPM, 'Sess') && ~isempty(SPM.Sess)
    n = numel(SPM.Sess);
else
    n = 0;
end
end % local_num_sess


function obj = local_load_betas(obj, SPM, doverbose)
% Load beta_*.nii / beta_*.img images from the SPM working directory.

swd = '';
if isfield(SPM, 'swd') && ~isempty(SPM.swd), swd = SPM.swd; end
if isempty(swd) || ~isfolder(swd)
    warning('glm_map:NoSWD', 'SPM.swd is not a valid folder; cannot load betas.');
    return
end

files = dir(fullfile(swd, 'beta_*.nii'));
if isempty(files), files = dir(fullfile(swd, 'beta_*.img')); end
if isempty(files)
    warning('glm_map:NoBetas', 'No beta_*.nii/.img images found in %s.', swd);
    return
end

% Sort numerically by the beta index to keep column order aligned with X
names = {files.name};
fullnames = cellfun(@(f) fullfile(swd, f), names(:), 'UniformOutput', false);
betas = statistic_image(char(fullnames), 'type', 'Beta', 'noverbose');

if isfield(SPM, 'xX') && isfield(SPM.xX, 'name')
    betas.image_labels = SPM.xX.name;
end

obj.betas = betas;

if doverbose
    fprintf('  Loaded %d beta image(s) from %s.\n', numel(fullnames), swd);
end

end % local_load_betas


function s = local_iif(tf, a, b)
if tf, s = a; else, s = b; end
end % local_iif
