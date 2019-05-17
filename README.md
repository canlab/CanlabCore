CanlabCore
==========

This repository contains core tools for MRI/fMRI/PET analysis from the Cognitive and Affective Neuorscience Lab (Tor Wager, PI) and our collaborators.  Many of these functions are needed to run other toolboxes, e.g., the CAN lab’s multilevel mediation and Martin Lindquist’s hemodynamic response estimation toolboxes. A brief introduction to the toolbox can be found <a href = "http://canlab.github.io">here</a>. 

The tools include object-oriented tools for doing neuroimaging analysis with simple commands and scripts that provide high-level functionality for neuroimaging analysis.  For example, there is an "fmri_data" object type that contains neuroimaging datasets (both PET and fMRI data are ok, despite the name). If you have created and object called my_fmri_data_obj, then plot(my_fmri_data_obj) will generate a series of plots specific to neuroimaging data, including an interactive brain viewer (courtesy of SPM software).  predict(my_fmri_data_obj) will perform cross-validated multivariate prediction of outcomes based on brain data.  ica(my_fmri_data_obj) will perform independent components analysis on the data, and so forth.

The repository also includes other useful toolboxes, including:
- fMRI design optimization using a genetic algorithm (OptimizeGA)
- fMRI hemodynamic response function estimation (HRF_Est_Toolbox2)
- fMRI analysis with Hierarchical Exponentially Weighted Moving Average change-point analysis (hewma_utility)
- Various fMRI diagnostics (diagnostics)
- Miscellaneous other tools and functions for visualizing brain data

Getting help and additional information:
------------------------------------------------------------
We have several sources of documentation for this toolbox:

1.  For a walk-through of a common basic processing pipeline, see our <a href='https://canlabreposguide.hackpad.com/CANLab-Repository-Guide-aGTiWJr0zbt'>hackpad</a>
2.  For function-by-function help documents on the Core Tools objects and functions, see the <a href = http://canlabcore.readthedocs.org/en/latest/>help pages on Readthedocs</a>.
3.  For brief, documented code examples of some specific functions, and a batch script system that uses the CanlabCore object-oriented tools for second-level neuroimaging analysis, see <a href='https://github.com/canlab/CANlab_help_examples'>CANlab_help_examples github repository</a>

The CANlab website is https://canlabweb.colorado.edu/, and we also maintain a WIKI with more information on some of our toolboxes and fMRI analysis more generally, which is <a href = "https://canlabweb.colorado.edu/wiki/doku.php/help/fmri_tools_documentation">here</a>.  For more information on fMRI analysis generally, see <a href = "https://leanpub.com/principlesoffmri">Martin and Tor's online book</a> and our free Coursera videos and classes <a href = "https://www.coursera.org/learn/functional-mri">Principles of fMRI Part 1</a> and <a href = "https://www.coursera.org/learn/functional-mri-2">Part 2 </a>.

Dependencies: These should be installed to use this toolbox
------------------------------------------------------------
Matlab www.mathworks.com

Matlab statistics toolbox

Matlab signal processing toolbox

Statistical Parametric Mapping (SPM) software https://www.fil.ion.ucl.ac.uk/spm/

<recommended> matlab_bgl (graph theory) and spider (machine learning) toolboxes; these are included in this distribution
  
<recommended> the CANlab Neuroimaging_Pattern_Masks repository https://github.com/canlab/Neuroimaging_Pattern_Masks
  
<recommended> the canlab_help_examples repository  https://github.com/canlab/CANlab_help_examples
  
  
