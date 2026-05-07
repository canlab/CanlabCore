function tests = canlab_test_table
%CANLAB_TEST_TABLE Building tables from a thresholded statistic_image.

tests = functiontests(localfunctions);
end


function test_table_runs_on_thresholded_t(tc)
% table() prints to the command window and may return a region array.
% We only check that it runs without erroring on a permissive threshold.
t = canlab_get_sample_thresholded_t();
tc.verifyWarningFree(@() table(t, 'nolegend'));
end


function test_table_returns_something_useful(tc)
% Capture an output if table() returns one.
t = canlab_get_sample_thresholded_t();
try
    r = table(t, 'nolegend');
    if ~isempty(r)
        tc.verifyTrue(isa(r, 'region') || isstruct(r) || istable(r));
    end
catch ME
    tc.verifyFail(['table() should not error: ' ME.message]);
end
end
