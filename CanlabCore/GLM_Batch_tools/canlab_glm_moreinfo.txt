OVERVIEW
canlab_glm_* are a set of functions to handle the running
of glm analyses of fmri data. They use SPM for analysis of
individual subjects and Tor Wager's robust regression
functions (robfit) for group level analysis. Details can
be found in the description of the DSGN structure (below) 
and in the help statements for the functions, but here is a brief
overview:

Before running canlab_glm functions, you will need to have:
- preprocessed functional data
- modeling files (SPM style "conditions" and "regressors" MAT-files)
- a DSGN structure loaded in your workspace or saved as a MAT-file


PREPROCESSING
You can use canlab_preproc_2012 for your preprocessing needs. Whatever 
you use, you must end up with 4D functional image files that have the
same name for each subject and the relative path to each subject's directory.
See DSGN.funcnames


SPM "REGRESSORS"
A simple script can make multiple regressors files, one 
for each run, from the products of canlab_preproc_2012
(several wagerlab studies have used modified version of the same 
script, make_noise_model1.m, for this purpose). If needed, you
can also include regressors of interest in your multiple regressors
files. You can also save multiple files and call a different one
for each analysis you run.
See DSGN.multireg


SPM "CONDITIONS"
Similarly, you must make conditions files, typically using
data in your experiment logs (e.g., .edat files from E-Prime).
This is often the most complicated process in running an
analysis and varies the most from one experiment to the next.
See DSGN.contrasts


LOWER LEVELS IN SPM
With these files in place, you will then need to assemble a
DSGN structure (probably using a script) that will be passed into 
canlab_glm_subject_levels (defines all aspects of the design/model
for the analysis).

canlab_glm_subject_levels will take DSGN as a
variable (from the workspace) or as a MAT-file, so your setup
script can either save DSGN in a MAT-file or simply load the
DSGN variable into the workspace.

Ex:
>> setup_pr_model2
>> canlab_glm_subject_levels(DSGN, 'email', 'ruzic@colorado.edu')
OR (if saved as a MAT-file)
>> canlab_glm_subject_levels('pr_model2.mat')

Make sure you check out the available options in 
help canlab_glm_subject_levels.


GROUP LEVELS WITH ROBFIT
With subject levels completed, you can run group levels with
robfit using canlab_glm_group_levels. This function is built
with defaults that can be overridden. If you did lower levels
using canlab_glm_subject_levels, you can reuse the DSGN structure
(or MAT-file containing it).

>> canlab_glm_group_levels(DSGN, 'o', '../second_level/pr_model2')

(By default it the output will be written to the same
directory as the subject level analyses. The standard
wagerlab organization is to keep group level analyses
together in a separate directory, "second_level", hence
the 'o' option in the above command.)

Make sure you check out the available options in 
help canlab_glm_subject_levels.
