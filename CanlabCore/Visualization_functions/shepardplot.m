function [Y,obs,imp,stress] = shepardplot(D, ntest, varargin)
% :Usage:
% ::
%
%    [Y,obs,imp,stress] = shepardplot(D, [ntest], [k dims to save])
%
% :Outputs:
%
%   **Y:**
%        is stimulus locations (cols are dims, rows are stimuli)
%        = eigenvectors of scalar product matrix from cmdscale
%
%   **obs:**
%        is vector of observed distances
%
%   **imp:**
%        is vector of implied distances
%
%   squareform(obs) or (imp) = full matrix of distances
%
% ..
%    by Tor Wager, 8/1/04
%    Modified April 2010: ntest can be a vector of which dims to run
% ..

criterion = 'sstress';

% classical MDS
[Y] = cmdscale(D);
if nargin < 2 || isempty(ntest), ntest = size(Y,2); end

if isscalar(ntest), ntest = 1:ntest; end

% plot stress by number of dimensions
% for i = 1:size(Y,2)
%     s(i) = stress(D,Y,i);
% end

obs = squareform(D);

str = sprintf('Testing dim %03d',0); fprintf(1,str);
warning off
stress = zeros(1, max(ntest));

for i = ntest
    fprintf(1,'\b\b\b%03d',i);
    [Y,stress(i)] = mdscale(obs,i,'Criterion',criterion);
end
warning on
erase_string(str);
disp(' ')
disp(' ')

% get number of dims to save: Elbow
% g = gradient(s)
% ng = length(g);
% for i = 1:ng-1
% gd(i) = sum(g(1:i)) ./ sum(g(i+1:end));
% end

f1 = create_figure('Shepardplot', 1, 2);

plot(ntest, stress(ntest),'k-o','LineWidth',2,'MarkerFaceColor',[.5 .5 .5])
xlabel('Number of dimensions'); ylabel('Stress')
title('Stress by dimensions in model','FontSize',18)

% vector of observed distances - for diagnostics and clustering
% dist = D - tril(D);
% dist = dist';
% dist = dist(:); dist(dist == 0) = [];
% obs = dist;
if length(varargin) > 0
    k = varargin{1};
else
    k = [];
end

while isempty(k)
    k = input('Choose n NDMS dims: ');
end
disp(['Dimensions: ' num2str(k)]);

plot_vertical_line(k);

if k>5
    nrand = 250;
else 
    nrand = 100;
end
maxiter = 5000;
disp(['Running mdscale with ' num2str(nrand) ' starting configurations, max iter = ' num2str(maxiter)]);
[Y,s,err] = mdscale(obs,k,'Replicates', nrand,'Options',statset('MaxIter',maxiter),'Criterion',criterion);

% implied distances are the Euclidean distances between stimulus
% coordinates (sc).  Stim coordinates are the rows of eigenvectors (V).
% distances are given by pdist, or the outer product of an appropriately
% scaled matrix.

%implied = Y(:,1:k) * Y(:,1:k)';         % implied crossproducts (scalar products),
% implied = squareform(pdist(Y(:,1:k)));  % implied distances
%
% dist = implied - tril(implied);
% dist = dist';
% dist = dist(:); dist(dist == 0) = [];
% imp = dist;

imp = pdist(Y);

subplot(1,2,2); set(gca,'FontSize',16)
hold on;
%plot(imp,obs,'ko','MarkerFaceColor',[.5 .5 .5]);
%xlabel('Implied Euclidean distances'); ylabel('Observed distances')
%title(sprintf('Distances in %2.0f dimensions',k),'FontSize',18)

[dum,ord] = sortrows([err(:) obs(:)]);
plot(obs,imp,'ko','MarkerFaceColor',[.5 .5 .5]);
hold on; plot(obs(ord),err(ord),'k.-');
xlabel('Dissimilarities'); ylabel('Distances/Disparities')


return



function erase_string(str1)
fprintf(1,repmat('\b',1,length(str1))); % erase string
return
