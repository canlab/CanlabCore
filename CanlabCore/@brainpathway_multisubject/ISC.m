function bs = ISC(bs)
% Calculate brain-wide inter-subject correlations in patterns across parcel means
%
% :Usage:
% ::
%
%     bs = ISC(bs)
%
% For objects: Type methods(object_name) for a list of special commands
%              Type help object_name.method_name for help on specific
%              methods.
%
% ..
%     Author and copyright information:
%
%     Copyright (C)2020 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Inputs:
%
%   **bs:**
%        A brainpathway_multisubject object
%
%
% :Optional Inputs:
%   **'nofactor':**
%        Omit factor analysis
%
%
% :Outputs:
%
%   **bs:**
%        A brainpathway_multisubject object
%        Output is added to relevant object properties:
%
%         bs.connectivity.regions.isc = intersubject correlation matrix;
%         bs.connectivity.regions.isc_subj_unusualness = 1/mean isc with others for each subject. 1/isc = [0 Inf] (bounded at zero);
%         bs.connectivity.regions.isc_outliers = subjects with isc_subj_unusualness > 2 SD above the mean
%         bs.connectivity.regions.isc_numfactors = estimated number of factors from elbow of PCA eigenvalues
%         bs.connectivity.regions.isc_Lambda = Subject factor loadings, can be correlated with other variables
%         bs.connectivity.regions.isc_Lambda_descrip = 'Lambdas are n_subjects x numfactors factor loadings, describing dimensions of subject similarity here';
%
%
% :Examples:
% ::
%
%    % None yet - example from HCP data
%
% :References:
%
% :See also:
%   - methods(brainpathway_multisubject)
%

% ..
%    Programmers' notes:
%    12/2020 Created by Tor Wager
% ..


fprintf('Flattening conn matrices. ')
mat = bs.flatten_conn_matrices('replacenans');

% Find NaNs and replace with subject mean (so regions will have min=zero influence on isc)
% Done now in flatten_conn_matrices
% fprintf('Imputing missing values. ')
%
% whnan = isnan(mat);
% subjmean = nanmean(mat');
%
% for i = 1:size(mat, 1)
% 	mat(i, whnan(i, :)) = subjmean(i);
% end

fprintf('Calculating ISC. ')
isc = corrcoef(double(mat'));

create_figure('region isc');
imagesc(isc); colorbar
drawnow

%whnan = all(isnan(isc), 1);
% anynans = sum(any(whnan, 2));

% Subject unusualness (potential outliers)
% --------------------------
meanisc = mean(isc) - (1./size(isc, 1)); % subtract the 1 on the diagonal
isc_subj_unusualness = 1 ./ meanisc;

isc_outliers = isc_subj_unusualness > (mean(isc_subj_unusualness) + 2 * std(isc_subj_unusualness));

figure; plot( 1 ./ meanisc, 'bo', 'MarkerFaceColor', [.3 .3 1]);


% Factor analysis
% --------------------------
fprintf('Factor analysis ')

% Choose # dimensions: where acceleration is greatest, i.e., derivative of
% eigenvalues changes fastest and the curve flattens at the "elbow"
[eigvec, eigval] = pcacov(isc);
gg = gradient(gradient(eigval));
[~, isc_numfactors] = max(gg);

fprintf('(%d) factors. ', isc_numfactors)

Lambda = factoran(isc, isc_numfactors, 'Xtype', 'cov');

if isc_numfactors > 1
    
    create_figure('region isc');
    scatterhist(Lambda(:, 1), Lambda(:, 2));
    xlabel('Factor 1')
    ylabel('Factor 2')
    
end

fprintf('\n');
disp('Added results to bs.connectivity.regions.isc');

bs.connectivity.regions.isc = isc;
bs.connectivity.regions.isc_subj_unusualness = isc_subj_unusualness;
bs.connectivity.regions.isc_outliers = isc_outliers;
bs.connectivity.regions.isc_numfactors = isc_numfactors;
bs.connectivity.regions.isc_Lambda = Lambda;
bs.connectivity.regions.isc_Lambda_descrip = 'Lambdas are n_subjects x numfactors factor loadings, describing dimensions of subject similarity here';

% bs.connectivity.regions.isc = intersubject correlation matrix;
% bs.connectivity.regions.isc_subj_unusualness = 1/mean isc with others for each subject. 1/isc = [0 Inf] (bounded at zero);
% bs.connectivity.regions.isc_outliers = subjects with isc_subj_unusualness > 2 SD above the mean
% bs.connectivity.regions.isc_numfactors = estimated number of factors from elbow of PCA eigenvalues
% bs.connectivity.regions.isc_Lambda = Subject factor loadings, can be correlated with other variables
% bs.connectivity.regions.isc_Lambda_descrip = 'Lambdas are n_subjects x numfactors factor loadings, describing dimensions of subject similarity here';

end
