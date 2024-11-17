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
