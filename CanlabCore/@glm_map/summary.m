function summary(obj)
% Print a one-screen summary of a glm_map object.
%
% Shows the design summary (via disp) plus any accumulated warnings and the
% provenance history.
%
% :Usage:
% ::
%
%     summary(obj)
%
% :Inputs:
%
%   **obj:**
%        A glm_map object.
%
% :See also:
%   - glm_map, diagnostics
%
% ..
%    Programmers' notes:
%    2026 - Initial implementation.
% ..

disp(obj);

if ~isempty(obj.warnings)
    fprintf('  Warnings:\n');
    for i = 1:numel(obj.warnings)
        fprintf('    - %s\n', obj.warnings{i});
    end
    fprintf('\n');
end

if ~isempty(obj.history)
    fprintf('  History:\n');
    for i = 1:numel(obj.history)
        fprintf('    %d. %s\n', i, obj.history{i});
    end
    fprintf('\n');
end

end % summary
