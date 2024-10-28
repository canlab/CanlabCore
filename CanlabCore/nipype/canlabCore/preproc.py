from nipype.interfaces.matlab import MatlabCommand
from nipype.interfaces.base import (
    TraitedSpec,
    BaseInterface,
    BaseInterfaceInputSpec,
    File,
)
from traits.api import List, Bool
import os
from string import Template


class OutliersInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, 
        desc="path to NIFTI like object")
    movement_file = File(exists=True, mandatory=False,
        desc="Path to movement file that can be read by readtable() in matlab")
    movement_order = List([1,2,3,4,5,6], usedefault=True,
        desc="Indices of movement_file corresponding to rot_x, rot_y, rot_z, trans_x, trans_y, trans_z, in that order")
    movement_radians = Bool(default_value=True, usedefault=True, 
        desc="Are rotations specified in radians?")
    out_file = File('outliers.csv', usedefault=True,
        desc="File indicate indices of volume to censor")


class OutliersOutputSpec(TraitedSpec):
    out_file = File(exists=True)


class Outliers(BaseInterface):
    """Use canlabCore's image_vector/outliers to "find outliers based on a 
       combination of rmssd (DVARS), robust spatial variance (MAD), Mahalanobis 
       distance based on covariance and correlation matrices, and images with 
       >25% missing values

       Serves as a substitute for rapidart. Apply it to preprocessed inputs. In
       particular make sure no further volumes are discarded (e.g. early acquisitions)
       after running this or the indices it returns will not match your data

       Inputs ::

           in_file - path to a nifti(.gz) image to evaluate

           out_file - name of output file to write out csv content to [default: outliers.csv]

       Optional Inputs ::

           movement_file - path to movement file to use for framewise displacement calculation.
                      Must contain rot_x, rot_y, rot_z, trans_x, trans_y, trans_z, in that 
                      order. Rotation in degrees and translation in mm. Expects a csv file or
                      equivalent with no header row.

           movement_order - list specifying a reordering of columns in movement file. E.g.
                      [4,5,6,1,2,3] will swap the location of the first 3 columns with the 
                      location of the next three columns if they're not in the expected order.
                      Any columns >6 are dropped because they're unspecified.

           movement_radians - boolean value indicating whether or not rotations are specified
                      in radians. If True values are multiplied by pi/180
    """

    input_spec = OutliersInputSpec
    output_spec = OutliersOutputSpec

    def _run_interface(self, runtime):
        d = dict(in_file=self.inputs.in_file,
                 movement_file=self.inputs.movement_file,
                 movement_order=self.inputs.movement_order,
                 movement_radians=self.inputs.movement_radians,
                 out_file=self.inputs.out_file)
        # This is your MATLAB code template
        script = Template(
            """in_file = '$in_file';
                             out_file = '$out_file';
                             movement_file = '$movement_file';
                             varargin={};
                             if exist(movement_file,'file')
                                 fd = table2array(readtable(movement_file));
                                 fd = fd(:,str2num('$movement_order'));
                                 if strcmp('$movement_radians','False')
                                     fd(:,1:3) = pi/180*fd(:,1:3);
                                 end
                                 varargin={'fd',fd};
                             end
                             [~,~,outlier_tables] = outliers(fmri_data(in_file), varargin{:});
                             outlier_ind = find(any(outlier_tables.outlier_regressor_matrix_corr,2))
                             csvwrite(out_file, outlier_ind);
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
        return result.runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        return outputs
