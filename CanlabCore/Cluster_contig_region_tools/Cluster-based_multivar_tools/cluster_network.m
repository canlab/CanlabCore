function [subcl,out,pccl] = cluster_network(cl,k,beh,varargin)
% [subcl,out,pccl] = cluster_network(cl,k,beh,varargin)
%
% uses clusters(i).BARPLOT.dat (preferred) or .data as data, with column k
% beh is n x m matrix of m covariates, (n subjects), double centered
%
% tor wager
%
% if varargin, saves data in file with string = varargin{1}
%
% subcl is clusters cell array with one clusters structure per class -
% separates positively and negatively weighted clusters.
%
% pccl is cl cell array for each significant PC
% Z field in each cluster has PC weights , for imaging
% Z values are actually correlation between timeseries and PC, which is
% a normalization of the PC weights (eigenvector weights)
% try montage_clusters([],pccl,[2 2])

doresid = 1;
dorobust = 0;   % nonfunctional!  change in prplot

if length(varargin) > 0
    dosave = 1;
    savestr = varargin{1};
    try, eval(['mkdir ' savestr]),catch,end
    cd(savestr)
else
    dosave = 0;
    savestr = '-';
end

plotflag = 1;

dat = []; dopk = 0;

for i = 1:length(cl)
    
    if ~isfield(cl(i).BARPLOT,'dat') & ~isfield(cl(i).BARPLOT,'data')
        dopk = 1;
        % not implemented
    end
    
    try
        dat = [dat cl(i).BARPLOT.dat(:,k)];
    catch
        dat = [dat cl(i).BARPLOT.data(:,k)];   
    end
end

% double-center
%
dat = scale(dat,1); dat = scale(dat',1)';

% ---------------------------------------------
% principal components
% ---------------------------------------------

%out = rapca(dat);
%out.pcomps = out.T;
%out.weights = out.P;
%out.eigval = out.L;
fprintf(1,'pca_npm.')
[pc,stats] = pca_npm(dat,500);
out = stats; out.T = stats.score(:,stats.wh); out.P = pc; out.L =stats.eigval(stats.wh);
figure('Color','w');plot(out.eigval,'ro-','LineWidth',2); hold on; plot(out.thresh,'ks-','LineWidth',2)
legend({'Eigenvalues' 'Upper 95% permuted'})

out.k = k;
out.dat = dat;
out.beh = beh;

for i = 1:size(out.T,2), 
    pccl{i} = cl;
    for j = 1:length(cl)
        pccl{i}(j).Z = ones(size(pccl{i}(j).Z)) .* out.wcor(j,i);
    end
    montage_clusters([],pccl{i},[2 2])
end


% ---------------------------------------------
% classification of regions
% ---------------------------------------------

maxclusters = length(out.L);
        
%tmp = out.P * diag(out.L);   % put original scaling back so that early eigs are weighted more heavily
%out.class = docluster(tmp',maxclusters,plotflag);
        
%out.outliers = find(out.od > mean(out.od) + 1.5*std(out.od) | out.od < mean(out.od) - 1.5*std(out.od));

disp('cluster_network:')
%fprintf(1,'Input clusters:\t%3.0f\nGroups:\t%3.0f\n',length(cl),max(out.class))
fprintf(1,'Input clusters:\t%3.0f\nGroups:\t%3.0f\n',length(cl),size(out.class,2))
fprintf(1,'Data column:\t%3.0f\n',k)
fprintf(1,['Eigenvalues:\t']) 
    for i = 1:min(8,length(out.eigval))
        if out.sig(i),sstr = '*';,else,sstr='';,end
        fprintf(1,'%3.2f%s\t',out.eigval(i),sstr)
    end
fprintf(1,'\n')
    
fprintf(1,['Var. explained:\t' repmat('%3.2f\t',1,length(out.L)) '\n'],out.expl(out.wh))
fprintf(1,['p (nonparametric):\t' repmat('%3.6f\t',1,length(out.L)) '\n'],out.p(out.wh))

    
c = tril(corrcoef(dat));c=c(:);c(c==1 | c == 0) = [];
disp(sprintf('PC1stats\tmean\tstd\tmin\tmax\n'))
disp(sprintf('%s\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t',savestr,mean(c),std(c),min(c),max(c)))
    
    % ---------------------------------------------
    % correlations between components and behavior
    % ---------------------------------------------
    
    for j = 1:size(beh,2),
        figure('Color','w');[r,str,sig,ry,rx,h] = prplot(out.T,beh,j);
        set(gcf,'Position',[38    48   618   495])
        title(['Partial correlations btwn beh.' num2str(j) ' and PCs'])
        out.CORREL.rcomp{j} = r;
        out.CORREL.compsig{j} = sig;
        out.CORREL.rcomp_names{j} = ['Corr. of beh. ' num2str(j) ' with PCs'];
        
        out.CORREL.raw_str = 'Rows are components, cols are beh vectors';
        for k = 1:size(out.T,2)
            r = corrcoef(out.T(:,k),beh(:,j));
            out.CORREL.raw_compr(k,j) = r(1,2);
        end
        
        if dosave & sig,
            saveas(gcf,[savestr 'BEH' num2str(j) '_PCs.fig'])
            saveas(gcf,[savestr 'BEH' num2str(j) '_PCs.tif'])
        end

    end
   
    %
    %for i = 1:max(out.class)
    %    out.classdata(:,i) = mean(dat(:,out.class==i),2);
    %end
            
        
    % ---------------------------------------------     
    % correlations between avgs and behavior
    % ---------------------------------------------
    
    for j = 1:size(beh,2),
        figure('Color','w');[r,str,sig,ry,rx,h] = prplot(out.classdata,beh,j);
        set(gcf,'Position',[957.0000   64.0000  618.0000  495.0000])
        title(['Partial Correlations btwn beh.' num2str(j) ' and avg data of each class'])
        out.CORREL.beh = beh;
        out.CORREL.rclass{j} = r;
        out.CORREL.classsig{j} = sig;
        out.CORREL.rclass_names{j} = ['Corr. of beh. ' num2str(j) ' with class avgs'];
        
        for k = 1:size(out.classdata,2)
            r = corrcoef(out.classdata(:,k),beh(:,j));
            out.CORREL.raw_classr(k,j) = r(1,2);
        end
            
        if dosave,
            saveas(gcf,[savestr 'BEH' num2str(j) '_class.fig'])
            saveas(gcf,[savestr 'BEH' num2str(j) '_class.tif'])
        end
    end

    
    % ---------------------------------------------     
    % table
    % ---------------------------------------------
    
    fprintf(1,['Correlations\t'])
    for i = 1:size(beh,2)
        for j = 1:size(out.T,2)
            fprintf(1,'BEH%3.0fCOMP%3.0f\t',i,j)
        end
    end
    fprintf(1,'\n')
    fprintf(1,['r: component scores\t\n'])
    for i = 1:size(beh,2)
        fprintf('Partial r with Beh %3.0f\t',i)
        for j = 1:size(out.T,2)
            if out.CORREL.compsig{i}(j), rstr='*';,else,rstr='';,end
            fprintf(1,'%3.2f%s\t',out.CORREL.rcomp{i}(j),rstr)
        end
        fprintf(1,'\n')
    end
    fprintf(1,'\n')
    
    [dummy,dummy,dummy,dummy,rcrit]=r2z(.5,size(out.T,1),.05);
    out.CORREL.rcrit = rcrit;
    for i = 1:size(beh,2)
        fprintf(1,'Raw r with Beh %3.0f\t',i)
        for j = 1:size(out.T,2)
            if out.CORREL.raw_compr(j,i) >= rcrit, rstr='*';,else,rstr='';,end
            fprintf(1,'%3.2f%s\t',out.CORREL.raw_compr(j,i),rstr)
        end
        fprintf(1,'\n')
    end
    
    fprintf(1,['r: class averages\t\n'])
    for i = 1:size(beh,2)
        fprintf('Partial r with Beh %3.0f\t',i)
        for j = 1:size(out.class,2)
            if out.CORREL.classsig{i}(j), rstr='*';,else,rstr='';,end
            fprintf(1,'%3.2f%s\t',out.CORREL.rclass{i}(j),rstr)
        end
        fprintf(1,'\n')
    end
    
    for i = 1:size(beh,2)
        fprintf(1,'Raw r with Beh %3.0f\t',i)
        for j = 1:size(out.class,2)
            if out.CORREL.raw_classr(j,i) >= rcrit, rstr='*';,else,rstr='';,end
            fprintf(1,'%3.2f%s\t',out.CORREL.raw_classr(j,i),rstr)
        end
        fprintf(1,'\n')
    end
    
    %%% ***  class data predicting beh
    %    for i = 1:size(beh,2)
    %    fprintf(1,'Raw r with Beh %3.0f\t',i)
    %    for j = 1:size(out.class,2)
    %        if out.CORREL.raw_classr(j,i) >= rcrit, rstr='*';,else,rstr='';,end
    %        fprintf(1,'%3.2f%s\t',out.CORREL.raw_classr(j,i),rstr)
    %    end
   %     fprintf(1,'\n')
   %end
    
    fprintf(1,'\n')
    fprintf(1,['prop. regions in class\t\n'])
    tmp = sum(out.class)./size(out.class,1);
    fprintf(1,repmat('%3.2f\t',1,length(tmp)),tmp)
    fprintf(1,'\n')
    fprintf(1,['r: avg class weights\n'])
    for i = 1:length(out.wh)
        for j = 1:size(out.class,2)
            fprintf(1,'%3.2f\t',mean(out.pc(:,i) .* out.class(:,j)));
        end
    end
        
    fprintf(1,'\n')
    fprintf(1,'\n')
    
    
    % ---------------------------------------------     
    % save clusters and get figures for each
    % ---------------------------------------------
    mycols = {'r' 'b' 'g' 'y' 'c' 'm' 'k' 'k' 'k' 'k'}; str = [];
    for i = 1:size(out.class,2)
        %subcl{i} = cl(out.class==i);
        subcl{i} = cl(find(out.class(:,i)));
        str = [str ',subcl{' num2str(i) '}'];
    end
    str = ['montage_clusters([]' str ',[mycols(1:size(out.class,2)) {''k''} {''k''}]);'];
    eval(str)
    set(gcf,'Position',[266          73         711        1048])
    if dosave,
            saveas(gcf,[savestr 'class_montage.fig'])
            saveas(gcf,[savestr 'class_montage.fig'])
    end
        
    % ---------------------------------------------     
    % residual table
    % ---------------------------------------------
    wh_pos = zeros(size(out.T,1),size(beh,2)); wh_neg = wh_pos;
    
    fprintf(1,['Cluster\tx\ty\tz\tVox\t'])
    for i=1:length(out.L), fprintf(1,'%s\t',['R-C' num2str(i)]);,end
    for i = 1:size(out.class,2),fprintf(1,'Class_%3.0f\t',i),end
 
    fprintf(1,'Comm.\t')
    for i=1:size(beh,2), fprintf(1,'%s\t',['Res-BEH' num2str(i)]);,end
    fprintf(1,'\n')
    
    for i = 1:length(cl)
        
        if isfield(cl,'shorttitle'),mstr=[cl(i).shorttitle '(cl. ' num2str(i) ')'];
        elseif isfield(cl,'BAstr'), mstr = cl(i).BAstr;, else, mstr = ['CL' num2str(i)];,end
        
        fprintf(1,'%s\t%3.0f\t%3.0f\t%3.0f\t%3.0f\t',mstr,cl(i).mm_center(1),cl(i).mm_center(2),cl(i).mm_center(3),cl(i).numVox)
        fprintf(1,repmat('%3.2f\t',1,length(out.L)),out.wcor(i,:))
        for j = 1:size(out.class,2),fprintf(1,'%3.0f\t',out.class(i,j)),end
        
        if doresid
            
        X = [out.T ones(size(out.T,1),1)]; y = dat(:,i);
        rd = y - X * (pinv(X) * y);
        out.residstr = 'Residual activation after removing principal components';
        out.residuals(:,i) = rd;
        out.communality(i) = 1 - (var(rd) ./ var(y));
        fprintf(1,'%3.2f\t',out.communality(i));
        
        for j = 1:size(beh,2)
            figure('Color','w');[r,str,sig,ry,rx,h] = prplot(rd,beh,j);
            set(gcf,'Position',[1068         701         524         410])
            title(['RESIDUAL r: BEH' num2str(j) ' cl ' num2str(i) ' for data col. ' num2str(k)])
            
            if dosave,
                saveas(gcf,[savestr 'BEH' num2str(j) 'cl' num2str(i) '_residcor.fig'])
                saveas(gcf,[savestr 'BEH' num2str(j) 'cl' num2str(i) '_residcor.tif'])
            end
            
            if sig, 
                rstr='*';,
                if r > 0,wh_pos(i,j) = 1;,elseif r < 0,wh_neg(i,j) = 1;,else, error('Uh-oh!'),end
            else,rstr='';,close,
            end
            fprintf(1,'%3.2f%s\t',r,rstr);
        end

        end
    
        fprintf(1,'\n')
    end

    if doresid
    % ---------------------------------------------     
    % show areas with residual correlations
    % ---------------------------------------------

    for i = 1:size(beh,2)
        mycols = {[1 .5 0] [.2 .6 .5]}; str = []; 
        
        clpos{i} = cl(find(wh_pos(:,i)));
        clneg{i} = cl(find(wh_neg(:,i)));
        
        if ~isempty(clpos{i})
            str = [str ',clpos{' num2str(i) '}'];
        else
            mycols = mycols(2);
        end
                
        if ~isempty(clneg{i})
            str = [str ',clneg{' num2str(i) '}'];
        end

        if ~isempty(clpos{i}) |  ~isempty(clneg{i})
            str = ['montage_clusters([]' str ',mycols);'];
            eval(str)
            set(gcf,'Position',[266          73         711        1048])
    
            if dosave,
                saveas(gcf,[savestr 'BEH' num2str(i) '_resid_montage.fig'])
                saveas(gcf,[savestr 'BEH' num2str(i) '_resid_montage.tif'])
            end
        end
    end
    
    end % if doresid
    
    if dosave, eval(['save ' savestr '_clusters subcl out']),end
    

    
return
    
    
    
    
    
    function class = docluster(a,maxclusters,doplot)

    Y = pdist(a','euclid');     % transpose so the voxels are observations, eigenvectors the variables
    Z = linkage(Y,'complete');
    if maxclusters > 1
        class = cluster(Z,maxclusters)';
        
        if doplot, 
            dendrogram(Z,0); title('Dendrogram for clustering')
        end
    
    else
        class = ones(1,size(a,2));
    end
    
    
    
return