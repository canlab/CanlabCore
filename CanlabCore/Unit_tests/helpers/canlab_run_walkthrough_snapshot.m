function canlab_run_walkthrough_snapshot(tc, wt__snapshot_path)
%CANLAB_RUN_WALKTHROUGH_SNAPSHOT Run a frozen walkthrough snapshot under CI hardening.
%
%   canlab_run_walkthrough_snapshot(tc, snapshot_path)
%
%   Executes a copied CANlab_help_examples walkthrough script (a "snapshot"
%   bundled under Unit_tests/walkthroughs/private/) one %%-cell at a time, in
%   a single shared workspace, with headless-graphics and fault-tolerant error
%   handling. Reports outcomes to the supplied matlab.unittest TestCase `tc`.
%
%   Why snapshots instead of running the real tutorials. The walkthroughs in
%   the CANlab_help_examples repo are teaching scripts, full of graphics
%   (orthviews, surface), optional data sets, and the occasional interactive
%   prompt. They should stay clean tutorials, not be retrofitted with test
%   scaffolding, and CanlabCore CI should not depend on an external repo's
%   un-hardened scripts. So we keep verbatim copies here and apply the
%   hardening at run time, in one reusable place (this function) rather than
%   editing each copy. Refresh a snapshot by overwriting the file under
%   private/ from example_help_files/ — no re-hardening needed.
%
%   Hardening applied:
%     * Headless graphics: DefaultFigureVisible is forced 'off' for the run
%       (callers also set this; doing it here makes the harness safe to call
%       directly). Figures are closed between sections.
%     * Section isolation: the script is split at its %% cell boundaries and
%       each cell is run in its own try/catch, so a graphics-only section that
%       fails on a headless runner does not abort the compute sections.
%     * Error classification (see canlab_classify_environment_error):
%         - graphics / input / data / cascade errors  -> section SKIPPED,
%           recorded, execution continues.
%         - genuine errors -> the harness stops and the test FAILS with an
%           informative message (which section, error id, message, and the
%           offending top stack frame), because downstream cells usually
%           cascade once a real error occurs.
%     * Outcome mapping to matlab.unittest:
%         - any genuine failure          -> tc.verifyFail (test FAILS)
%         - every section skipped (env)  -> tc.assumeFail (test INCOMPLETE)
%         - at least one section ran ok  -> PASS, with a logged skip summary.
%
%   Variable hygiene: every internal variable in this function is prefixed
%   `wt__` because snapshot code is eval'd in this workspace and would
%   otherwise clobber ordinary loop variables (i, n, t, r, ...).
%
% :Inputs:
%   **tc:**             a matlab.unittest.TestCase (from functiontests).
%   **snapshot_path:**  absolute path to the snapshot .m file to run.
%
% :See also: canlab_classify_environment_error, canlab_run_all_tests

% ..
%    Author: CANlab. Part of the CanlabCore Unit_tests headless-CI harness.
% ..

tc.assertTrue(exist(wt__snapshot_path, 'file') == 2, ...
    sprintf('Walkthrough snapshot not found: %s', wt__snapshot_path));

[~, wt__name] = fileparts(wt__snapshot_path);
wt__src       = fileread(wt__snapshot_path);
wt__sections  = local_split_sections(wt__src);
wt__nsec      = numel(wt__sections);

% Headless: callers set this too, but be self-sufficient.
wt__prev_vis = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
wt__cleanup  = onCleanup(@() set(0, 'DefaultFigureVisible', wt__prev_vis));

wt__skipped      = {};     % cellstr of human-readable skip notes
wt__ran_real     = 0;      % count of sections that executed without error
wt__had_env_skip = false;  % drives cascade classification downstream
wt__failure      = [];     % first genuine failure, if any (struct)

for wt__s = 1:wt__nsec

    wt__code = wt__sections{wt__s};

    % Skip cells that are pure comment/whitespace, or that define a
    % function (snapshots are scripts; a trailing function would not eval).
    if local_is_noop(wt__code)
        continue
    end

    try
        evalc(wt__code);                 % runs in THIS workspace; vars persist
        wt__ran_real = wt__ran_real + 1;
    catch wt__ME
        wt__cat = canlab_classify_environment_error(wt__ME, wt__had_env_skip);
        if strcmp(wt__cat, 'genuine')
            wt__failure = struct( ...
                'section', wt__s, ...
                'id',      wt__ME.identifier, ...
                'message', wt__ME.message, ...
                'where',   local_top_frame(wt__ME), ...
                'snippet', local_first_stmt(wt__code));
            break                        % stop; later cells will cascade
        else
            wt__had_env_skip = true;
            wt__skipped{end+1, 1} = sprintf('  section %d/%d [%s]: %s', ...
                wt__s, wt__nsec, wt__cat, local_oneline(wt__ME.message)); %#ok<AGROW>
        end
    end

    close all force
end

% ---------------------------------------------------------------------
% Map collected outcomes onto the TestCase.
% ---------------------------------------------------------------------
if ~isempty(wt__failure)
    wt__report = sprintf([ ...
        '%s: genuine error at section %d of %d.\n' ...
        '  identifier: %s\n' ...
        '  message:    %s\n' ...
        '  section starts: %s\n' ...
        '  at:         %s\n' ...
        '  (%d earlier section(s) skipped for environment reasons)'], ...
        wt__name, wt__failure.section, wt__nsec, ...
        wt__failure.id, local_oneline(wt__failure.message), ...
        wt__failure.snippet, wt__failure.where, numel(wt__skipped));
    tc.verifyFail(wt__report);

elseif wt__ran_real == 0
    tc.assumeFail(sprintf( ...
        '%s: all %d section(s) skipped for environment reasons:\n%s', ...
        wt__name, wt__nsec, strjoin(wt__skipped, newline)));

else
    if ~isempty(wt__skipped)
        tc.log(1, sprintf('%s: %d of %d section(s) skipped (environment):\n%s', ...
            wt__name, numel(wt__skipped), wt__nsec, strjoin(wt__skipped, newline)));
    end
    tc.verifyTrue(true, sprintf('%s ran %d section(s) without genuine error', ...
        wt__name, wt__ran_real));
end

end


% =====================================================================
% Local helpers
% =====================================================================

function sections = local_split_sections(src)
% Split source at %% cell markers. Each section keeps its leading marker
% line (a comment, harmless to eval). Everything before the first marker is
% section 1.
lines  = regexp(src, '\r\n|\r|\n', 'split');
is_head = ~cellfun('isempty', regexp(lines, '^\s*%%', 'once'));
starts = find(is_head);
if isempty(starts) || starts(1) ~= 1
    starts = [1, starts];
end
starts = unique(starts);
ends   = [starts(2:end) - 1, numel(lines)];
sections = cell(numel(starts), 1);
for k = 1:numel(starts)
    sections{k} = strjoin(lines(starts(k):ends(k)), newline);
end
end


function tf = local_is_noop(code)
% True if the cell has no executable content: blank, comment-only, or a
% function definition (which cannot be eval'd as a statement).
stripped = regexprep(code, '%[^\n]*', '');          % drop line comments
stripped = strtrim(stripped);
tf = isempty(stripped) || ~isempty(regexp(stripped, '^function\b', 'once'));
end


function s = local_oneline(msg)
% Collapse a multi-line error message to a single trimmed line.
s = strtrim(regexprep(msg, '\s+', ' '));
end


function s = local_top_frame(ME)
% "funcname (line N)" for the deepest stack frame that is not this harness
% (errors raised at a section's top level otherwise just point back at the
% evalc call site here, which is unhelpful). Falls back to a marker.
if isempty(ME.stack)
    s = '<top level of section>';
    return
end
this_fn = mfilename;
for k = 1:numel(ME.stack)
    if ~strcmp(ME.stack(k).name, this_fn)
        s = sprintf('%s (line %d)', ME.stack(k).name, ME.stack(k).line);
        return
    end
end
s = '<top level of section>';
end


function s = local_first_stmt(code)
% First non-comment, non-blank line of a section, trimmed, for the report.
lines = regexp(code, '\r\n|\r|\n', 'split');
s = '<none>';
for k = 1:numel(lines)
    ln = strtrim(lines{k});
    if ~isempty(ln) && ~strncmp(ln, '%', 1)
        s = ln;
        if numel(s) > 100, s = [s(1:97) '...']; end
        return
    end
end
end
