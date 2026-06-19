function tests = canlab_test_display
%CANLAB_TEST_DISPLAY Display methods on a thresholded t-image (smoke tests).
%
% These exercise the visualization API without verifying pixel content.
% A test passes if the call returns without erroring (and any figures
% created are cleaned up). The GUI-heavy methods (orthviews, surface)
% are guarded so they skip cleanly in headless `matlab -batch` sessions.

tests = functiontests(localfunctions);
end


function setup(tc) %#ok<DEFNU>
% Close any figures left around by other tests.
close all force;
end


function teardown(tc) %#ok<DEFNU>
close all force;
end


function test_montage_runs_on_thresholded_t(tc)
t = canlab_get_sample_thresholded_t();
tc.verifyWarningFree(@() montage(t));
end


function test_slices_runs_on_thresholded_t(tc)
t = canlab_get_sample_thresholded_t();
% slices() prints output and creates a figure; we only assert no error.
fn = @() slices(t);
tc.verifyWarningFree(@() runAndClose(fn));
end


function test_orthviews_runs_on_thresholded_t(tc)
% orthviews uses SPM's spm_orthviews and requires a Graphics window.
% Skip cleanly when MATLAB has no desktop / display.
tc.assumeTrue(usejava('desktop') || usejava('jvm') && feature('ShowFigureWindows'), ...
    'orthviews requires an interactive figure window');
t = canlab_get_sample_thresholded_t();
tc.verifyWarningFree(@() orthviews(t));
end


function test_surface_runs_on_thresholded_t(tc)
% surface() renders meshes; in pure batch this can fail on missing OpenGL.
tc.assumeTrue(usejava('jvm'), 'surface requires Java for rendering');
t = canlab_get_sample_thresholded_t();
try
    surface(t);
catch ME
    if any(strcmp(ME.identifier, ...
            {'MATLAB:graphics:opengl:Unavailable', 'MATLAB:graphics:initialize'})) ...
            || contains(lower(ME.message), 'opengl') ...
            || contains(lower(ME.message), 'display')
        tc.assumeFail(['surface() needs a graphics environment: ' ME.message]);
    else
        rethrow(ME);
    end
end
end


function runAndClose(fn)
fn();
close all force;
end
