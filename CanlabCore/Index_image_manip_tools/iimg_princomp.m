function [volInfo,clusters] = iimg_princomp(maskname,image_names) 
% :Usage:
% ::
%
%     [volInfo,clusters] = iimg_princomp(maskname,image_names) 


[dat,volInfo] = iimg_get_data(maskname,image_names); % Get data and principal components

[eigvec, compscore, eigval] = princomp(dat,'econ');
figure;plot(eigval,'ko-');

nc = input('how many comps to save? ');
clusters = cell(1,nc);

volInfo.eigval = eigval;
volInfo.eigvec = eigvec(:,1:nc);
volInfo.score = compscore(:,1:nc);

%% loop thru components
xbar = mean(dat);

for i = 1:nc
    k = i;

    %% flip component
    comp = eigvec(:,k); cscore = compscore(:,k);
    X = comp; X(:,end+1) = 1;
    b = pinv(X) * xbar';
    if b(1) < 0,
        comp = -comp; cscore = -cscore;
    end

    volInfo.eigvec(:,k) = comp;
    volInfo.score(:,k) = cscore;

    %% correlation with component
    x = corrcoef([cscore dat]);
    x = x(1,2:end)';

    volInfo.corr_with_comps(:,k) = x;

    %% get high voxel loadings
    sigvox = (abs(x) > .5) .* sign(x);
    
    avgregion = dat * sigvox;
    voldata = iimg_reconstruct_3dvol(sigvox,volInfo);
    cl = mask2clusters(voldata,volInfo.mat);
    cluster_orthviews(cl,'bivalent');

    clusters{k} = cl;
    volInfo.component_names{k} = input('Name this component: ','s');
    volInfo.high_loadings(:,k) = sigvox;
    volInfo.averages(:,k) = avgregion;
end
%% end loop

return

