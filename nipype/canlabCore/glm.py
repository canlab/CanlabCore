from nipype.interfaces.matlab import MatlabCommand
from nipype.interfaces.base import (
    TraitedSpec,
    BaseInterface,
    BaseInterfaceInputSpec,
    File,
)
from nipype.utils.filemanip import save_json
import nipype.interfaces.fsl as fsl # for FSLCommand and its input/out specs
from traits.api import List, Bool, Float, Int, Any, Either
import os
from string import Template


'''
VIFs requires the following folders on your matlab path:

CanlabCore/CanlabCore/diagnostics/
CanlabCore/CanlabCore/Visualization_functions/
CanlabCore/CanlabCore/OptimizeDesign11/core_functions/
CanlabCore/CanlabCore/Statistics_tools/
CanlabCore/CanlabCore/Data_processing_tools/
CanlabCore/CanlabCore/Misc_utilities/
'''

class VIFsInputSpec(BaseInterfaceInputSpec):
    spm_mat_file = File(exists=True, mandatory=True, 
        desc="SPM.mat file produced by nipype.interfaces.spm.model.EstimateModel interface")
        
    events_only = Bool(False, mandatory=False, usedefault=True,
        desc="Produce plots and diagnostics for only events, not nuisance covariantes or other user-specified regressors")
        
    from_multireg = Either(None, Int(), 
        usedefault=True,
        requires=['events_only'],
        desc="followed by an integer to include n columns for the \'regressors\' matrix")
        
    vif_threshold = Either(None, Float, 
        usedefault=True, 
        desc="Only regressors with VIF > t will be printed in VIF table")
        
    sort_by_vif = Bool(False, usedefault=True,
        desc="Sort regressors in VIF table by VIF.")


class VIFsOutputSpec(TraitedSpec):
    vifs = File('vifs.csv', exists=True)
    
    png = File('Variance_Inflation.png', exists=True)
    
    hpfilt = File('High_pass_filter_analysis.png', exists=True)   
    


class VIFs(BaseInterface):
    """Use canlabCore's diagnostics/scn_spm_design_check to estimate vifs and generate timeseries
        diagnostic plots.
    """

    input_spec = VIFsInputSpec
    output_spec = VIFsOutputSpec

    def _run_interface(self, runtime):
        spm_dir = os.path.dirname(self.inputs.spm_mat_file)
        
        d = dict(spm_dir=spm_dir,
                 events_only=self.inputs.events_only,
                 from_multireg=self.inputs.from_multireg,
                 vif_threshold=self.inputs.vif_threshold,
                 sort_by_vif=self.inputs.sort_by_vif)

        # This is your MATLAB code template
        script = Template(
            """spm_dir = '$spm_dir';
             varargin = {};
             if strcmp('$events_only', 'True')
                varargin{end+1} = 'events_only';
             end
             if ~strcmp('$from_multireg', 'None')
                varargin{end+1} = 'from_multireg';
                varargin{end+1} = int(str2num('$from_multireg'));
             end
             if ~strcmp('$vif_threshold', 'None')
                varargin{end+1} = 'vif_threshold';
                varargin{end+1} = str2num('$vif_threshold');
             end
             if strcmp('$sort_by_vif', 'True')
                varargin{end+1} = 'sort_by_vif';
             end

             vifs = scn_spm_design_check(spm_dir, varargin{:});

             tbl = table(vifs.allvifs(:), 'RowNames', vifs.name);
             writetable(tbl, 'vifs.csv', 'WriteRowNames', true, 'WriteVariableNames', false);
             exit;
          """
        ).substitute(d)

        # mfile = True  will create an .m file with your script and executed.
        # Alternatively
        # mfile can be set to False which will cause the matlab code to be
        # passed
        # as a commandline argument to the matlab executable
        # (without creating any files).
        # This, however, is less reliable and harder to debug
        # (code will be reduced to
        # a single line and stripped of any comments).
        mlab = MatlabCommand(script=script, mfile=True)
        result = mlab.run()
        
        
        self.vifs = os.path.abspath('vifs.csv')
        self.png = os.path.abspath('Variance_Inflation.png')
        self.hpfilt = os.path.abspath('High_pass_filter_analysis.png')
        
        
        return result.runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        for item in outputs.keys():
            outputs[item] = os.path.abspath(getattr(self,item))
        return outputs



class NVolsInputSpec(fsl.base.FSLCommandInputSpec):
    in_file = File(
        exists=True,
        argstr="%s",
        mandatory=True,
        position=1,
        desc="input file to generate stats of",
    )


class NVolsOutputSpec(TraitedSpec):
    nvols = Any(desc="nvols output")
            
class NVols(fsl.base.FSLCommand):
    """Use FSL fslnvol command to read out 4d vol count. Based on nipype FSL.ImageStats
    `FSL info
    <http://www.fmrib.ox.ac.uk/fslcourse/lectures/practicals/intro/index.htm#fslutils>`_


    Examples
    --------

    >>> from nipype.interfaces.fsl import NVol
    >>> from nipype.testing import funcfile
    >>> stats = NVols(in_file=funcfile)
    >>> stats.cmdline == 'fslvols %s'%funcfile
    True


    """

    input_spec = NVolsInputSpec
    output_spec = NVolsOutputSpec

    _cmd = "fslnvols"

    def _format_arg(self, name, trait_spec, value):
        return super()._format_arg(name, trait_spec, value)

    def aggregate_outputs(self, runtime=None, needed_outputs=None):
        outputs = self._outputs()
        # local caching for backward compatibility
        outfile = os.path.join(os.getcwd(), "nvols_result.json")
        if runtime is None:
            try:
                nvols = load_json(outfile)["nvols"]
            except OSError:
                return self.run().outputs
        else:
            nvols = int(runtime.stdout)
            save_json(outfile, dict(nvols=nvols))
        outputs.nvols = nvols
        return outputs


class FMRIConcatInputSpec(BaseInterfaceInputSpec):
    spm_mat_file = File(exists=True, mandatory=True, 
        desc="path to SPM file in main SPM directory (containing the betas)")
        
    nscan = List(
        mandatory=True,
        desc="(list (n,) indicating number of TRs of each constituent scan.")

class FMRIConcatOutputSpec(TraitedSpec):
    spm_mat_file = File(exists=True)


class FMRIConcat(BaseInterface):
    """
    Uses modified spm_fmri_concatenate function to adjust temporal 
    autocorrelation function for concatenation of scans
    """

    input_spec = FMRIConcatInputSpec
    output_spec = FMRIConcatOutputSpec

    def _run_interface(self, runtime):
        import shutil
        import os
        
        basename = os.path.basename(self.inputs.spm_mat_file)
        newSPM = os.path.join(os.getcwd(), basename)
        shutil.copyfile(self.inputs.spm_mat_file, newSPM)
        
        d = dict(spm_mat_file=newSPM,
                 nscan=[int(n) for n in self.inputs.nscan])

        # This is your MATLAB code template
        script = Template(
            """ spm_mat_file = '$spm_mat_file';
                nscan = $nscan;
                try
                    spm_fmri_concatenate(spm_mat_file, nscan);
                catch
                    % eventually this can replace spm_fmri_concatenate
                    % entirely, but it needs more testing for now.
                    warning('support for timeseries models using mixed concatenated/block diagonal designs is experimental. Check outputs of FMRIConcat outputs carefully, in particular SPM.xX properties.')
                    spm_fmri_concatenate_multisess(spm_mat_file, nscan);
                end
            """
        ).substitute(d)

        # mfile = True  will create an .m file with your script and executed.
        # Alternatively
        # mfile can be set to False which will cause the matlab code to be
        # passed
        # as a commandline argument to the matlab executable
        # (without creating any files).
        # This, however, is less reliable and harder to debug
        # (code will be reduced to
        # a single line and stripped of any comments).
        mlab = MatlabCommand(script=script, mfile=True)
        result = mlab.run()

        self.spm_mat_file = os.path.abspath(d['spm_mat_file'])
        
        base, ext = os.path.splitext(basename)
        os.remove(os.path.join(os.getcwd(), base + '_backup' + ext))

        return result.runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['spm_mat_file'] = os.path.abspath(self.spm_mat_file)
        return outputs
