function cluster_manova(clusters,fnames,varargin)
% :Usage:
% ::
%
%     cluster_manova(clusters,fnames,verbose)
%
% :Inputs:
%
%   **clusters:**
%        is output of clusters2database, with all fields
%        from database
%
%   **fnames:**
%        is cell array of strings with names to test
%        e.g., {'Rule' 'Task'}
%
%   uses stats toolbox
%
%   **verbose:**
%        optional, produces more output and tests
%
% ..
%    tor wager
% ..

if length(varargin) > 0, verb = 1; else, verb = 0;, end

mycol = {'bo' 'ys' 'g^' 'cv' 'rd' 'ys'};
warning off     % avoid openGL warnings on my machine
fprintf(1,'\nNote: df will not equal sums of individual class memberships if classes are fuzzy; \nMANOVA run separately on each class')
fprintf(1,'\n-------------------------------------------------\n* Cluster_Manova.m :: Table of Clusters \n-------------------------------------------------\n')
clusters = cluster_table(clusters);

con = [1 -1]; 


for i = 1:length(clusters)
    
    %x = [clusters(i).x clusters(i).y clusters(i).z];
    fprintf(1,'\n-------------------------------------------------\n* Cluster %3.0f :: %s\n-------------------------------------------------\n',i,clusters(i).BAstring)
    
    im = {};, clear allcent
    for j = 1:length(fnames)
        % get centers and distances
        eval(['im{j} = clusters(i).' fnames{j} ';'])
        cent = mean(im{j});
        allcent(j,:) = cent;
        if isempty(cent), cent = [NaN NaN NaN];,end
    end

    % differences
    for j = 1:length(con), im{j} = im{j} .* con(j);,end
    x = cat(3,im{:}); x = sum(x,3); % differences
    
    % distances between group centers
    fprintf(1,'\nDistances between class centers\n')
    for j = 1:length(fnames),fprintf(1,'%s\t',fnames{j}),end,fprintf(1,'\n')
    mydist = squareform(pdist(allcent,'euclid'));
   
    
    %fprintf(1,'%s\t%3.0f\t%3.0f\t%3.2f\t%3.4f\t%3.0f%\t',fnames{j},stats.dfB,stats.dfW,stats.lambda,p,round(100*missrate));
    fprintf(1,'Name\tnum_peaks\tdfB\tdfW\tWilk''s\tp\tmissclass. rate\t\n')
    
    
    for j = 1:length(fnames)
        
        d = 0;
        
        eval(['group = clusters(i).' fnames{j} ';'])
        if max(group) == 1 | any(group==0), group = group+1;,end
        
        % run manova if there are enough observations
        % print output line
        nvox = eval(['sum(clusters(i).' fnames{j} ');']);
        [d,stats] = manova(x,group,fnames{j},nvox);
        
        % print centers and distances
        
        fprintf(1,'%s\t%3.2f\t%3.2f\t%3.2f\t',fnames{j},allcent(j,1),allcent(j,2),allcent(j,3))
        fprintf(1,[repmat('%3.2f\t',1,size(allcent,1)) '\n'],mydist(j,:)');
        
        
            if d > 0 & verb
                % stuff to print if it's significant
                
                figure('Color','w'); hold on;
                %fprintf(1,'\n%s Dimension Weights: %3.2f\t%3.2f\t%3.2f\t', ...
                %    fnames{j},stats.eigenvec(1,1),stats.eigenvec(2,1),stats.eigenvec(3,1))
                
                fprintf(1,'\n%s Group Centers',fnames{j})
                fprintf(1,'\nGroup\tx\ty\tz\t')
                for k = 1:max(group), 
                    tmp = mean(x(group==k,:));
                    fprintf(1,'\n%3.0f\t%3.0f\t%3.0f\t%3.0f\t',k,tmp(1),tmp(2),tmp(3))
                    plot3(x(group==k,1),x(group==k,2),x(group==k,3),mycol{k})
                end
                addbrain
                camzoom(.7)
                title(['Cluster ' num2str(i) ' ' fnames{j} ': blue o = no, yellow square = yes'])
                set(gca,'FontSize',18)
                
            end 
            
            if d > 0    % confidence volume stuff
                
                if sum(group==2) > 12,
                    results = confidence_volume(x(group==2,:),mycol{j}(1));
                end
                
                hold on;
                plot3(x(group==2,1),x(group==2,2),x(group==2,3),mycol{j},'LineWidth',2);
                
            end
            
    end
        
    if verb            
    % add ALL cluster centers and pairwise manovas
    all_pairwise(clusters(i),fnames)
    end
    
    fprintf(1,'\n')
    
end

warning on

return



function [missrate,c] = get_missclassrate(group,class)

% misclassification rate with confusion matrix
for i = 1:max(group)
    for j = 1:max(class)
        c(i,j) = sum(group==i & class==j);
    end
end
missrate = 1 - trace(c) ./ sum(c(:));

return




function all_pairwise(clusters,fnames)

xyz = [clusters.x clusters.y clusters.z];

fprintf(1,'\n\nCenters for each class\n')

for i = 1:length(fnames),
    eval(['im(:,i) = clusters.' fnames{i} ';'])
    cent = mean(xyz(find(im(:,i)),:),1);
    
    if isempty(cent), fprintf(1,'No peaks of type %s in cluster\n', fnames{i}),cent = [NaN NaN NaN];
    else, fprintf(1,'%s\t%3.2f\t%3.2f\t%3.2f\t\n',fnames{i},cent(1),cent(2),cent(3))
    end
    allcent(i,:) = cent;
    
end

% distances between group centers
fprintf(1,'\nDistances between class centers\n')
for i = 1:length(fnames),fprintf(1,'%s\t',fnames{i}),end,fprintf(1,'\n')
mydist = squareform(pdist(allcent,'euclid'));
fprintf(1,[repmat('%3.2f\t',1,size(allcent,1)) '\n'],mydist');

% now do pairwise manovas
fprintf(1,'\nPairwise MANOVAs\n')
fprintf(1,'\nName\tdfB\tdfW\tWilk''s\tp\tmissclass. rate\t\n')

for i = 1:size(im,2)
    for j = (i+1):size(im,2)
        
        % build group identities, eliminate points in both groups
        group = im(:,i); group(find(im(:,j))) = 2; group(im(:,i) & im(:,j)) = 0;
        xyztmp = xyz; xyztmp(group == 0,:) = []; group(group == 0) = [];
        
        [d,stats] = manova(xyztmp,group,[fnames{i} ' vs. ' fnames{j}]);
           
    end
end


return




function [d,stats,missrate] = manova(x,group,myname,nvox)

d = 0;

        for k = 1:max(group), numg(k) = sum(group==k);,end
        
        if any(numg < 3), fprintf(1,'%s\t%3.0f\t Too few in cell to classify\t\t\t\t\t\t\t\t',myname,nvox);,
        elseif numg(1) == length(group), fprintf(1,'%s\t%3.0f\t No activations for one group in cluster\t\t\t\t\t\t\t\t',myname,nvox);,
            
        else
            class = classify(x,x,group);
            [d,p,stats] = manova1(x,group);
            
            %whos group, whos x, sum(group==1),sum(group==2), stats
            
            [missrate,c] = get_missclassrate(group,class);
        
            fprintf(1,'\n%s\t%3.0f\t%3.0f\t%3.0f\t%3.2f\t%3.4f\t%3.0f%%\t',myname,nvox,stats.dfB,stats.dfW,stats.lambda,p,round(100*missrate));
            fprintf(1,'%3.2f\t%3.2f\t%3.2f\t', ...
                    stats.eigenvec(1,1),stats.eigenvec(2,1),stats.eigenvec(3,1))
        end
        
 return
        
    
    
