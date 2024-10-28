## Nipype interfaces

This folder should contain interfaces for running canlabCore functions, and
in particular first level GLM modeling related functions, in nipype. Add it
to your python path like so,
```
package_directory = '<path_to_canlabCoreRepo>/CanlabCore/nipype/'
if package_directory not in sys.path:
    sys.path.insert(0, package_directory)
import canlaCore
```

These interfaces will run matlab code from canlabCore on data provided
through nipype and return outputs that are available for downstream use
within a nipype pipeline. In order for the interfaces to work they needs 
access to canlabCore functions.

If canlabCore is not on your matlab path, you need to add it or modify this
code with a call that adds it at runtime. You additionally will likely need
SPM on your path. I (BP) do not normally recommend adding 3rd party packages 
to your permanent path, since this invites code collisions in the long run 
(e.g. two 3rd party packages, unaware of one another, that use the same 
file/function names), but the alternative here would require user-specific 
updates to a community repo or alternatively making a copy of this repo and 
losing version tracking abilities. For details on how to do this neatly see 
here,

https://www.mathworks.com/help/matlab/matlab_env/add-folders-to-matlab-search-path-at-startup.html

The use of the MATLABPATH environment variable in particular looks promising
