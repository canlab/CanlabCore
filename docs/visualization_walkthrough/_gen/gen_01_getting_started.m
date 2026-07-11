function gen_01_getting_started
% Figures for "1. Getting started".
here = fileparts(mfilename('fullpath'));
addpath(here);
figdir = fullfile(fileparts(here), 'figures');
close all force

obj = load_image_set('emotionreg', 'noverbose');
t   = threshold(ttest(obj), .05, 'unc');

% orthviews (SPM) of the thresholded map
orthviews(t);
gwin = spm_figure('FindWin', 'Graphics');
if ~isempty(gwin), wt_save(gwin, fullfile(figdir, '01_orthviews.png')); end
close all force

disp('gen_01 DONE');
end
