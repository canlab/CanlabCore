function underlay = canlab_get_underlay_image(varargin)
% Get the name of an anatomical image to use as the underlay for orthviews and other displays
%
% No inputs: Use the default underlay
%
% Or enter an argument with one of these strings:
%
% 'spm2'                      'spm2_single_subj_T1_scalped.img'
% 'colin'                     'SPM8_colin27T1_seg.img'
% 'keuken'                    'keuken_2014_enhanced_for_underlay.img'
% 'icbm2009c', 'fmriprep20'   'fmriprep20_template.nii.gz'
% 'icbm2009c_0.5mm'           'tpl-MNI152NLin2009bAsym_res-1_T1w.nii.gz'
% 'mni152_1mm'                'MNI152NLin6Asym_T1_1mm.nii.gz'
% 'mni152_withskull'          'MNI152NLin6Asym_T1_1mm.nii.gz'
% 'mni152'                    'spm152.nii'
% ...or your custom underlay filename.
%
% Examples:
% underlay = canlab_get_underlay_image
% underlay = canlab_get_underlay_image('icbm2009c_0.5mm')

current_default = 'fmriprep20_template.nii';

if nargin == 0
    underlaystr = current_default;
else
    underlaystr = varargin{1};
    if isempty(underlaystr), underlaystr = current_default; end
end

switch underlaystr
    case 'spm2'
        underlayfilename = 'spm2_single_subj_T1_scalped.img';

    case 'colin'
        underlayfilename = 'SPM8_colin27T1_seg.img';

    case 'keuken'
        underlayfilename = 'keuken_2014_enhanced_for_underlay.img';

    case {'icbm2009c', 'fmriprep20'}
        underlayfilename = 'fmriprep20_template.nii.gz';

    case 'icbm2009c_0.5mm'
        underlayfilename = 'tpl-MNI152NLin2009bAsym_res-1_T1w.nii.gz';

    case {'mni152_1mm', 'mni152_withskull'}
        underlayfilename = 'MNI152NLin6Asym_T1_1mm.nii.gz';

    case 'mni152'
        underlayfilename = 'spm152.nii';

    otherwise
        underlayfilename = underlaystr;

        if iscell(underlayfilename)
            underlayfilename = underlayfilename{1};
            warning('underlay must be a string.')
        end

        if ~ischar(underlayfilename)
            error('underlay must be a string.')
        end
        
end

% Get full path to image
underlay = which(underlayfilename);

% Try to find it if it's missing

% Blanks
if ~exist(underlay, 'file')
    underlayfilename = deblank(underlayfilename);
    underlay = which(underlayfilename);
end

% Maybe curr dir?
if ~exist(underlay, 'file')
    [~, name, ext] = fileparts(underlayfilename);
    underlay = which([name ext]);
end

% Maybe strip .gz if missing this version?
if ~exist(underlay, 'file')
    [~, name] = fileparts(underlayfilename);
    underlay = which(name);
end

% maybe add .gz if missing non-gz version?
if ~exist(underlay, 'file')
    [~, name] = fileparts(underlayfilename);
    underlay = which([name '.gz']);
end

if isempty(underlay) || ~exist(underlay, 'file')
    disp('Cannot find underlay image:')
    disp(underlayfilename)

    disp('No valid underlay image specified, and I cannot find the default one.');
    disp(['The default is:' current_default ' and should be on your path.'])

    return
end

end % function
