function icadat = ica(fmridat_obj, varargin)
% Spatial ICA of an fmri_data object
%   - icadat = ica(fmridat_obj, [number of ICs to save])
%   - icadat is also an fmri_data object, with .dat field voxels x components
%
% :Notes:
%   - Spatial component maps are saved in columns of icadat.dat, or rows of icadat.dat'
%   - icasig = W * mixedsig
%   - icasig = icadat.dat' = W * fmridat_obj.dat'
%
% W is the separation matrix, the inverse of the mixing matrix
% A is the mixing matrix, the inverse of the separation matrix
% 
% A and W are stored in additional_info field of icadat in cells {1} and {2}, respectively
% A = ica_obj.additional_info{1};
% W = ica_obj.additional_info{2};
%
% Model:
% D = A * S, where:
%   D = [n x v] dataset, images x voxels in fmridat_obj.dat'
%   A = [n x k] mixing matrix (the inverse of the separation matrix if n = k)

%   W = [k x n] separation matrix, (the inverse of the mixing matrix if n = k)
%               W' is a matrix whose columns are image series, e.g., time
%               series. They reflect the expression of each component
%               across images and can be used to analyze differences across
%               image groups (conditions, time, etc.)
%   S = [k x v] components matrix, with k components. S = icadat.dat'
%               S = W * D;
%
% - Each image in D is a mixture of maps S defined by a row of A
% - Each column of A defines the expression of one map in S across images (e.g., time, or other image series)
%
% Reconstructed data:
% D_recon = A * S;
% figure; plot(D(:), D_recon(:), 'k.')
%   A is related to D * S' (fmridat_obj.dat' * icadat.dat)';
%
% Dual regression / back-reconstruction:
% B = pinv(S') * D';  % dual regression step 1, each data image n regressed on k spatial maps S; B = k x n
% B = pinv(icadat.dat) * fmridat_obj.dat;
%   Note: fmridat_obj can be a different image object here, e.g., an independent dataset
%   If fmridat_obj is the same dataset used to derive ICA components S, B is very similar to W
%
% S_hat = pinv(B') * fmridat_obj.dat';  % dual regression step 2
%   S_obj = fmridat_obj;
%   S_obj.dat = S_hat';
%   montage(get_wh_image(S_obj, 1:4));

figure; plot(B(:), W(:), 'k.')

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

