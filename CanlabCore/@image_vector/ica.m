function icadat = ica(fmridat_obj, varargin)
% Spatial ICA of an fmri_data object
%   - icadat = ica(fmridat_obj, [number of ICs to save])
%   - icadat is also an fmri_data object, with .dat field voxels x components
%
% :Notes:
%   - icasig = W * mixedsig
%   - icasig = icadat.dat' = W * fmridat_obj.dat'
%
% A is scaled version of fmridat_obj.dat' * icadat.dat
% 
% A and W are stored in additional_info field of icadat

nic = 30;
if length(varargin) > 0, nic = varargin{1}; end

neig =  size(fmridat_obj.dat, 2);

[icasig, A, W] = icatb_fastICA(double(fmridat_obj.dat'), 'lastEig', neig, ...
    'numOfIC', nic, 'stabilization', 'on', 'verbose', 'on');

icadat = fmri_data;

icadat.dat = icasig';
icadat.mask = fmridat_obj.mask;
icadat.volInfo = fmridat_obj.volInfo;
icadat.removed_voxels = fmridat_obj.removed_voxels;
icadat.removed_images = 0;

icadat.additional_info{1} = A;
icadat.additional_info{2} = W;
icadat.history{1} = 'Spatial ICA of fmri_data object.';
icadat.history{2} = 'Stored A and W matrices in cells 1,2 of additional_info';
icadat.source_notes = 'fastICA algorithm on images (spatial)';

% Show on orthviews

tmpicadat = icadat;
for i = 1:size(tmpicadat.dat, 2)
    d = tmpicadat.dat(:, i);
    tmpicadat.dat(d < prctile(d, 85) & d > prctile(d, 15), i) = 0;
    
end

orthviews(tmpicadat);
spm_orthviews_white_background

end


% m1 = prctile(icadat.dat(:), 5);
% m2 = prctile(icadat.dat(:), 95);
%
% spm_orthviews('window', 1:min(24, size(icadat.dat, 2)), [m1 m2]);

