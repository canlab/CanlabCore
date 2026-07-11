%% Load and explore a sample dataset
% This example shows how to load a pre-set example fMRI dataset into
% an fmri_data object.

%% General instructions
%
% Before you start, the CANlab_Core_Tools must be added to your path with
% subfolders. Otherwise, you will get errors.
%
% The script canlab_toolbox_setup can help download and install the
% toolboxes you need.
%
% Sample datasets are in the "Sample_datasets" folder in CANlab_Core_Tools.
% Many tutorials apply pre-trained patterns and masks. 
% These are stored in this Github repository:
% 
% <https://github.com/canlab/Neuroimaging_Pattern_Masks>
% 
% In addition, you can explore these and find more information here:
% <https://sites.google.com/dartmouth.edu/canlab-brainpatterns/home>
%
% This example will use emotion regulation data in the folder: 
% "Wager_et_al_2008_Neuron_EmotionReg"
% The dataset is a series of contrast images from N = 30 participants.
% Each image is a contrast image for [reappraise neg vs. look neg]
% 
% These data were published in:
% Wager, T. D., Davidson, M. L., Hughes, B. L., Lindquist, M. A., 
% Ochsner, K. N.. (2008). Prefrontal-subcortical pathways mediating 
% successful emotion regulation. Neuron, 59, 1037-50.
%
% Here are a couple of helpful functions we will use for display:
% (you can ignore these.)
dashes = '----------------------------------------------';
printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);

%% The fmri_data object
%
% The fmri_data object class is one of the most important, basic types of
% objects in the CANlab object-oriented toolbox. It stores image data in a
% matrix form, which is more space efficient and friendly for analysis with
% various software packages/algorithms.
%
% For philosophy, see:
% <https://canlab.github.io/>
%
% For more info on the object-oriented approach, see:
% <https://canlab.github.io/objectoriented/>
%
% When you call the _class constructor_ fmri_data(), this is what it does:
%
% <<fmri_data_object_diagram.png>>

%% Section 1: The quick and easy way to load a pre-specified dataset
%
% The function load_image_set has a number of pre-defined image sets that
% you can load with one simple command.  This is the easiest way to load
% sample data.  The images must be on your Matlab path for this to work.

[data_obj, subject_names, image_names] = load_image_set('emotionreg');

% data_obj is an fmri_data object containing all 30 images.
% subject_names is a list of short names for each image.
% image_names is a list of the full image names with their path.

%% Section 2: Manual load.  Use filenames to find file names 
%
% First, we need to list the file names in a string matrix.
% Then, we can load them into an fmri_data object
% We will use the filenames function to get the names, and the fmri_data
% object constructor function to load the data.

% First, check whether images are on your path:
% We will search for one image and save the path name.

printhdr('Check that we can find data images:')
myfile = which('con_00810001.img');
mydir = fileparts(myfile);
if isempty(mydir), disp('Uh-oh! I can''t find the data.'), else disp('Data found.'), end

% Now we can list all the file names.

printhdr('Find files and get their names:')
image_names = filenames(fullfile(mydir, 'con_008100*img'), 'absolute');
disp('Done.');

% Now load them into an fmri_data object.

printhdr('Loading the image data into the object:')
data_obj = fmri_data(image_names);

% This is the gateway to doing many other things, which are explained in
% other help files.  But just to get us started, let's run through a few
% basic things we can do. We'll mainly just look at some standard plots of
% the data.

%% Section 3: Plot the data we just loaded
%
% Operations that we can perform on fmri_data objects are called methods.
% You can see a list of methods by typing methods(data_obj).
% Here, we'll call the plot method to visualize the data.

plot(data_obj)
snapnow

%% Section 4: Get help on an object
%
% Objects have associated help files. 
% Type:

help fmri_data

% ...to get help for the fmri_data() class constructor and object.
%
% Some objects also have doc files that contain additional information.
% Try:

doc fmri_data

%% Section 5: Explore methods
%
% Objects have methods associated with them, or things you can do with
% them. To see a list of methods for the fmri_data object, type:

methods(fmri_data)

%%
% You can also use the name of an _instance_ of an object, i.e., a variable
% in the workspace. |methods(data_obj)| produces the same list.
% In addition, you can get help on any of the object's methods by typing:
% |object class name_dot_method name|. e.g.,

help fmri_data.mean

% Note: There may be trouble accessing these help files in Matlab R2019b.
% Works in 2018 and earlier. But maybe this is not a general problem...give
% it a try.

%% Section 6: Explore on your own
%
% 1. Try to run a couple of other methods on your fmri_data object. Some
% require additional inputs, but many do not. Look at the help for more
% info. What do you see? Can you return any meaningful output, and if so,
% what did the method do?
% 

% That's it for this section!!

%% Explore More: CANlab Toolboxes
% Tutorials, overview, and help: <https://canlab.github.io>
%
% Toolboxes and image repositories on github: <https://github.com/canlab>
%
% <html>
% <table border=1><tr>
% <td>CANlab Core Tools</td>
% <td><a href="https://github.com/canlab/CanlabCore">https://github.com/canlab/CanlabCore</a></td></tr>
% <td>CANlab Neuroimaging_Pattern_Masks repository</td>
% <td><a href="https://github.com/canlab/Neuroimaging_Pattern_Masks">https://github.com/canlab/Neuroimaging_Pattern_Masks</a></td></tr>
% <td>CANlab_help_examples</td>
% <td><a href="https://github.com/canlab/CANlab_help_examples">https://github.com/canlab/CANlab_help_examples</a></td></tr>
% <td>M3 Multilevel mediation toolbox</td>
% <td><a href="https://github.com/canlab/MediationToolbox">https://github.com/canlab/MediationToolbox</a></td></tr>
% <td>M3 CANlab robust regression toolbox</td>
% <td><a href="https://github.com/canlab/RobustToolbox">https://github.com/canlab/RobustToolbox</a></td></tr>
% <td>M3 MKDA coordinate-based meta-analysis toolbox</td>
% <td><a href="https://github.com/canlab/Canlab_MKDA_MetaAnalysis">https://github.com/canlab/Canlab_MKDA_MetaAnalysis</a></td></tr>
% </table>
% </html>
% 
% Here are some other useful CANlab-associated resources:
%
% <html>
% <table border=1><tr>
% <td>Paradigms_Public - CANlab experimental paradigms</td>
% <td><a href="https://github.com/canlab/Paradigms_Public">https://github.com/canlab/Paradigms_Public</a></td></tr>
% <td>FMRI_simulations - brain movies, effect size/power</td>
% <td><a href="https://github.com/canlab/FMRI_simulations">https://github.com/canlab/FMRI_simulations</a></td></tr>
% <td>CANlab_data_public - Published datasets</td>
% <td><a href="https://github.com/canlab/CANlab_data_public">https://github.com/canlab/CANlab_data_public</a></td></tr>
% <td>M3 Neurosynth: Tal Yarkoni</td>
% <td><a href="https://github.com/neurosynth/neurosynth">https://github.com/neurosynth/neurosynth</a></td></tr>
% <td>M3 DCC - Martin Lindquist's dynamic correlation tbx</td>
% <td><a href="https://github.com/canlab/Lindquist_Dynamic_Correlation">https://github.com/canlab/Lindquist_Dynamic_Correlation</a></td></tr>
% <td>M3 CanlabScripts - in-lab Matlab/python/bash</td>
% <td><a href="https://github.com/canlab/CanlabScripts">https://github.com/canlab/CanlabScripts</a></td></tr>
% </table>
% </html>
%
% *Object-oriented, interactive approach*
% The core basis for interacting with CANlab tools is through object-oriented framework.
% A simple set of neuroimaging data-specific objects (or _classes_) allows you to perform
% *interactive analysis* using simple commands (called _methods_) that
% operate on objects. 
%
% Map of core object classes:
%
% <<CANlab_object_types_flowchart.png>>

