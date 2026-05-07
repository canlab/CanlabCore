CanlabCore
==========

[![tests](https://github.com/canlab/CanlabCore/actions/workflows/test.yml/badge.svg)](https://github.com/canlab/CanlabCore/actions/workflows/test.yml)
[![tests-walkthroughs](https://github.com/canlab/CanlabCore/actions/workflows/tests-walkthroughs.yml/badge.svg)](https://github.com/canlab/CanlabCore/actions/workflows/tests-walkthroughs.yml)

This repository contains core tools for MRI/fMRI/PET analysis from the Cognitive and Affective Neuorscience Lab (Tor Wager, PI) and our collaborators.  Many of these functions are needed to run other toolboxes, e.g., the CAN lab’s multilevel mediation and Martin Lindquist’s hemodynamic response estimation toolboxes. An introduction to the toolbox and its philosophy, along with walkthroughs and tutorials, can be found at <a href = "http://canlab.github.io">canlab.github.io</a>. 

The tools include object-oriented tools for doing neuroimaging analysis with simple commands and scripts that provide high-level functionality for neuroimaging analysis.  For example, there is an "fmri_data" object type that contains neuroimaging datasets (both PET and fMRI data are ok, despite the name). If you have created and object called my_fmri_data_obj, then plot(my_fmri_data_obj) will generate a series of plots specific to neuroimaging data, including an interactive brain viewer (courtesy of SPM software).  predict(my_fmri_data_obj) will perform cross-validated multivariate prediction of outcomes based on brain data.  ica(my_fmri_data_obj) will perform independent components analysis on the data, and so forth.

📖 **[Object methods reference →](docs/Object_methods.md)** — class-by-class index of every method, with runnable examples and sample figures. The fastest way to learn the API.

The repository also includes other useful toolboxes, including:
- fMRI design optimization using a genetic algorithm (OptimizeGA)
- fMRI hemodynamic response function estimation (HRF_Est_Toolbox2)
- fMRI analysis with Hierarchical Exponentially Weighted Moving Average change-point analysis (hewma_utility)
- Various fMRI diagnostics (diagnostics)
- Miscellaneous other tools and functions for visualizing brain data

Getting help and additional information:
------------------------------------------------------------
Sources of documentation for this toolbox:

1. **[canlab.github.io](https://canlab.github.io/)** — top-level entry point with Setup, Repositories, Philosophy, and Batch system for 2nd-level analysis.
2.  **[Walkthroughs](https://canlab.github.io/walkthroughs/)** — step-by-step analysis tutorials with code. "How to do stuff with Canlab tools".
3. **[Tutorials](https://canlab.github.io/tutorials/)** — longer-form tutorials with more equations and theoretical explanation.
4. **[Object methods reference](docs/Object_methods.md)** class-by-class index (`fmri_data`, `image_vector`, `statistic_image`, `atlas`, `region`, ...) of object methods with per-function code maps, runnable examples and sample figures for selected functions. The fastest way to learn the API.
5. **[CANlab_help_examples](https://github.com/canlab/CANlab_help_examples)** repository — runnable MATLAB scripts (`example_help_files/`) and HTML output with figures (`published_html/`) for specific functions and workflows published in Walkthroughs and Tutorials (plus some extra stuff). This repo also contains the code for the second-level Batch script system described on canlab.github.io.
6. Function-by-function help documents and examples are provided in each function. In Matlab type >>help <function_name>
7. Older function-by-function help docs are in <a href = http://canlabcore.readthedocs.org/en/latest/>help pages on Readthedocs</a>.

For more information on fMRI analysis generally, see <a href = "https://leanpub.com/principlesoffmri">Martin Lindquist and Tor Wager's online book</a> and our free Coursera videos and classes <a href = "https://www.coursera.org/learn/functional-mri">Principles of fMRI Part 1</a> and <a href = "https://www.coursera.org/learn/functional-mri-2">Part 2 </a>. New in 2026: Martin and Tor's expanded fMRI methods book <a href = "https://mitpress.mit.edu/9780262045049/elements-of-functional-magnetic-resonance-imaging/">Elements of fMRI</a>. 

Dependencies: These should be installed to use this toolbox
------------------------------------------------------------
Matlab www.mathworks.com

Matlab statistics toolbox

Matlab signal processing toolbox

Statistical Parametric Mapping (SPM) software https://www.fil.ion.ucl.ac.uk/spm/

<recommended> matlab_bgl (graph theory) and spider (machine learning) toolboxes; these are included in this distribution
  
<recommended> the CANlab Neuroimaging_Pattern_Masks repository https://github.com/canlab/Neuroimaging_Pattern_Masks
  
<recommended> the canlab_help_examples repository  https://github.com/canlab/CANlab_help_examples
  
  
