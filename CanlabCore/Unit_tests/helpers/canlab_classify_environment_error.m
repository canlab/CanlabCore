function category = canlab_classify_environment_error(ME, had_env_skip)
%CANLAB_CLASSIFY_ENVIRONMENT_ERROR Bucket a caught error for CI test handling.
%
%   category = canlab_classify_environment_error(ME)
%   category = canlab_classify_environment_error(ME, had_env_skip)
%
%   Classifies an MException caught while running toolbox example/walkthrough
%   code into one of a small set of buckets, so that test harnesses can decide
%   whether to SKIP (mark Incomplete) or FAIL. The goal is to keep headless CI
%   from going red over missing-display / missing-data / interactive-prompt
%   conditions that are not real regressions, while still surfacing genuine
%   code errors.
%
%   This centralizes the heuristics previously hand-coded inside
%   canlab_test_help_examples.m (skip_on_environment_error) so the walkthrough
%   harness and the per-method help-example tests share one definition.
%
% :Inputs:
%
%   **ME:**
%        an MException (or any struct with .identifier, .message, and
%        optionally .stack fields).
%
%   **had_env_skip:**
%        logical (default false). When true, "undefined variable/function"
%        and "unrecognized field" errors are treated as CASCADE fallout from
%        an earlier skipped section rather than genuine failures. Pass the
%        running "have we already skipped something upstream" flag here.
%
% :Output:
%
%   **category:**
%        one of:
%          'graphics' - needs a display / OpenGL / Java / figure window, or
%                       the error originates inside a graphics-only code path
%                       (orthviews, surface rendering, etc.).
%          'input'    - code tried to prompt for interactive user input,
%                       unavailable in batch/CI.
%          'capability' - an UndefinedFunction error for a known optional
%                       MATLAB toolbox function (e.g. niftiinfo) that is not
%                       provisioned on this runner.
%          'data'     - an optional data file (signature, atlas, feature set)
%                       is not on the CI path.
%          'cascade'  - an undefined-variable/field error following an earlier
%                       environment skip (downstream fallout, not a new bug).
%          'genuine'  - none of the above; treat as a real failure.
%
% :See also: canlab_run_walkthrough_snapshot, canlab_test_help_examples

% ..
%    Author: CANlab. Part of the CanlabCore Unit_tests headless-CI harness.
% ..

if nargin < 2 || isempty(had_env_skip)
    had_env_skip = false;
end

msg = lower(ME.message);
id  = '';
if isprop(ME, 'identifier') || isfield(ME, 'identifier')
    id = ME.identifier;
end

% ---------------------------------------------------------------------
% Graphics: explicit graphics error ids, telltale message tokens, or a
% stack frame inside a known graphics-only function. spm_orthviews and the
% surface renderers fail or behave erratically on a headless runner.
% ---------------------------------------------------------------------
gfx_ids = {'MATLAB:graphics:opengl:Unavailable', ...
           'MATLAB:graphics:initialize', ...
           'MATLAB:class:InvalidHandle'};
is_gfx = any(strcmp(id, gfx_ids)) || ...
         strncmp(id, 'MATLAB:hg:', 10) || ...        % handle-graphics errors
         contains(msg, 'opengl') || contains(msg, 'display') || ...
         contains(msg, 'java')   || contains(msg, 'jvm')     || ...
         contains(msg, 'figure window') || contains(msg, 'no display') || ...
         contains(msg, 'graphics') || ...
         contains(msg, 'invalid or deleted object');  % set/get on a handle
                                                       % whose graphics section
                                                       % was skipped upstream

gfx_fns = {'spm_orthviews', 'orthviews', 'spm_check_registration', ...
           'spm_figure', 'surface', 'render_on_surface', 'addbrain', ...
           'cluster_surf', 'isosurface', 'riverplot', 'tor_3d', ...
           'render_blobs', 'cluster_orthviews'};
if ~is_gfx && isfield_or_prop(ME, 'stack') && ~isempty(ME.stack)
    stack_names = {ME.stack.name};
    if any(ismember(stack_names, gfx_fns))
        is_gfx = true;
    end
end
if is_gfx
    category = 'graphics';
    return
end

% ---------------------------------------------------------------------
% Interactive input: input() / keyboard prompts are unavailable in -batch.
% MATLAB reports MissingRequiredCapability with a "support for user input"
% message in that case.
% ---------------------------------------------------------------------
is_input = strcmp(id, 'MATLAB:services:MissingRequiredCapability') || ...
           contains(msg, 'support for user input') || ...
           contains(msg, 'input is not available');
if is_input
    category = 'input';
    return
end

% ---------------------------------------------------------------------
% Missing optional MATLAB toolbox. An UndefinedFunction error for a known
% MathWorks toolbox function (e.g. niftiinfo/niftiread from the Image
% Processing Toolbox) means that toolbox is not provisioned on this runner,
% not that the code is broken. Keep this list to functions that are *only*
% ever provided by a toolbox, so we never mask a genuine missing-CanlabCore
% function bug.
% ---------------------------------------------------------------------
optional_toolbox_fns = {'niftiinfo', 'niftiread', 'niftiwrite', ...
                        'cfg_getfile'};
if strcmp(id, 'MATLAB:UndefinedFunction')
    undef = regexp(ME.message, "Undefined function '([^']+)'", 'tokens', 'once');
    if ~isempty(undef) && any(strcmpi(undef{1}, optional_toolbox_fns))
        category = 'capability';
        return
    end
end

% ---------------------------------------------------------------------
% Missing optional data files (NPS+ signatures, Neurosynth feature set,
% Bianciardi atlas sources, etc.). load_image_set / annotate_* abort with a
% disp() + bare error('Exiting'), so the message is literally "Exiting".
% ---------------------------------------------------------------------
bare_exiting = isempty(id) && ...
               ismember(strtrim(strrep(ME.message, '.', '')), {'Exiting', 'exiting'});
is_data = contains(msg, 'cannot find images') || ...
          contains(msg, 'find and add the file') || ...
          contains(msg, 'not found in matlab path') || ...
          contains(msg, 'no such file or directory') || ...
          contains(msg, 'cannot find') || ...
          bare_exiting;
if is_data
    category = 'data';
    return
end

% ---------------------------------------------------------------------
% Cascade: an earlier section was skipped, so a variable/field it would
% have defined is now missing. Treat as fallout, not a new regression.
% ---------------------------------------------------------------------
if had_env_skip && ( ...
        strcmp(id, 'MATLAB:UndefinedFunction') || ...
        strcmp(id, 'MATLAB:undefinedVarOrClass') || ...
        strcmp(id, 'MATLAB:nonExistentField') || ...
        contains(msg, 'undefined') || ...
        contains(msg, 'unrecognized field'))
    category = 'cascade';
    return
end

category = 'genuine';

end


function tf = isfield_or_prop(obj, name)
% MException is an object (use isprop); a plain struct uses isfield.
tf = (isobject(obj) && isprop(obj, name)) || (isstruct(obj) && isfield(obj, name));
end
