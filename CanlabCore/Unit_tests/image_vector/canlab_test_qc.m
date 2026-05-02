function tests = canlab_test_qc
%CANLAB_TEST_QC Quality-control utilities.

tests = functiontests(localfunctions);
end


function test_descriptives_runs(tc)
obj = canlab_get_sample_fmri_data();
out = descriptives(obj);
tc.verifyTrue(isstruct(out));
end


function test_qc_metrics_second_level_runs(tc)
% qc_metrics_second_level computes group-level QC stats; smoke test only.
% It can emit benign warnings (missing optional fields), so we don't
% require warning-free — only that the call returns without an error.
obj = canlab_get_sample_fmri_data();
try
    qc_metrics_second_level(obj, 'noplot');
    tc.verifyTrue(true);
catch ME
    tc.verifyFail(['qc_metrics_second_level errored: ' ME.message]);
end
end
