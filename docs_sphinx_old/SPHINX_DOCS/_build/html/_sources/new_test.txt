
New Test
===============

Installation
------------

Model_building_tools
---------------
.. automodule:: Model_building_tools
   :members:

Dependencies:

* NumPy/SciPy
* pandas
* NiBabel
* ply (optional, for complex structured queries)
* scikit-learn (optional, used in some classification functions)

Usage
-----

Running analyses in Neurosynth is pretty straightforward. We're working on a user manual; in the meantime, you can take a look at the code in the /examples directory for an illustration of some common uses cases (some of the examples are in IPython Notebook format; you can view these online by entering the URL of the raw example on github into the online `IPython Notebook Viewer <http://nbviewer.ipython.org>`_--for example `this tutorial <http://nbviewer.ipython.org/urls/raw.github.com/neurosynth/neurosynth/master/examples/neurosynth_demo.ipynb>`_ provides a nice overview). The rest of this Quickstart guide just covers the bare minimum.

The NeuroSynth dataset resides in a separate submodule. If you installed Neurosynth directly from PyPI (i.e., with pip install), and don't want to muck around with git or any source code, you can manually download the data files from the `neurosynth-data repository <http://github.com/neurosynth/neurosynth-data>`_. The latest dataset is always stored in current_data.tar.gz in the root folder. Older datasets are also available in the archive folder.

