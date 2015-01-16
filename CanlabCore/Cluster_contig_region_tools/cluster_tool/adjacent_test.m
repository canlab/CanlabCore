function [ind]=adjacent_test(XYZ,varargin)

% USAGE:
%
% [ind]=adjacent_test(XYZ)
% [ind]=adjacent_test(...,criterion)
% [ind]=adjacent_test(...,criterion,min_adj)
% [ind]=adjacent_test(...,criterion,min_adj,adjacency)
% [ind]=adjacent_test(...,criterion,min_adj,adjacency,erosions)
%
% Note that inputs must be specified in the correct order.
%
% Returns a vector of same length as the non-tripleton dimension of XYZ,
% with numeric values indicating which cluster each voxel belongs too.
% Thus, assuming each column refers to a voxel coordinate,
% XYZ(:,ind==1) is all of the voxels in one cluster, and
% XYZ(:,ind==2) is another.
%
% adjacency is a string, which must be 'surface' (6 connectivity scheme),
% 'edge' (18 connectivity scheme), or 'corner' (26 connectivity scheme).
% Default is 'edge'.
%
% criterion (must be scalar) is the number of voxels a voxel must be
% adjacent too for it to be considered a 'core' cluster voxel (all voxels
% adjacent to 'core' voxels will be included in a cluster, but those that
% are adjacent to voxels that are part of a cluster but are not themselves
% 'core' voxels will not be included in the cluster). Setting criteria
% above the connectivity value associated with a given adjacency criteria
% will result in no clusters. Default is 0.
%
% min_adj (must be scalar) is the minimum number of voxels that a voxel
% must be adjacent to for it to be included in a cluster. Voxels adjacent
% to less voxels than min_adj will recieve a value of 0 in ind. Setting
% criteria above the connectivity value associated with a given adjacency
% criteria will result in no clusters. Default is 0.
%
%
% erosions is the integer number of times that the input image will be
% eroded according to the criterion and adjacency inputs before 'core'
% cluster voxels are defined. Default is 1.
%
%
% Note that:
% ind=adjacent_test(XYZ,'edge'); %or ind=adjacent_test(XYZ,'edge',0,0,1);
% ind2=spm_clusters(XYZ);
% isequal(ind,ind2)
%
% ans =
%
%      1
%
%
%
%

persistent runN
runN = 1;

if isempty(varargin)
    criterion=0;
    min_adj=0;
    adjacency='edge';
    er = 1;
elseif length(varargin)<2
    criterion=varargin{1};
    min_adj=0;
    adjacency='edge';
    er = 1;
elseif length(varargin)<3
    criterion=varargin{1};
    min_adj=varargin{2};
    adjacency='edge';
    er = 1;
elseif length(varargin)<4
    criterion=varargin{1};
    min_adj=varargin{2};
    adjacency=varargin{3};
    er = 1;
else
    criterion=varargin{1};
    min_adj=varargin{2};
    adjacency=varargin{3};
    er = varargin{4};
end

if ~strcmp(adjacency,'corner')&&~strcmp(adjacency,'edge')&&~strcmp(adjacency,'surface')
    error('Adjacency input incorrectly specified')
end

if size(XYZ,1)~=3
    if size(XYZ,2)~=3
        error('XYZ matrix must have a tripleton dimension!')
    else
        XYZ=XYZ';
    end
end

if ~isequal(XYZ,round(XYZ))||any(XYZ(:)<1)
    error('XYZ input contains invalid or duplicate coordinates!')
end

if strcmp(adjacency,'edge')&&~criterion&&~min_adj&&er==1
    try ind=spm_clusters(XYZ);return,catch end
end



ind = voxel_test([],XYZ,adjacency,min_adj,criterion);
for i = 2:er
    ind = voxel_test(logical(ind),XYZ,adjacency,min_adj,criterion);
end

k=1;
while k<=max(ind)
    if ~any(ind==k)
        ind(ind>k)=ind(ind>k)-1;
        k=k-1;
    end
    k=k+1;
end



disp('Attaching stray voxels to existing clusters')

lost = [];
for i = 1:er
    attach = zeros(size(ind));
    for k=1:length(ind)

        if ~rem(k,10000)
            disp([int2str(((k+(i-1)*length(ind))/(er*length(ind)))*100) '% complete'])
        end

        if ~ind(k)
            
            adj = find_adjacent(k,XYZ,adjacency,[]);

            if ( sum(adj) > min_adj && ( sum(adj)>criterion || i==er ) ) && any(ind(adj))
                attach(k) = 1;
            end
        end
    end

    for k = 1:length(ind)
        if attach(k)

            adj = find_adjacent(k,XYZ,adjacency,[]);

            if any(ind(adj))
                if min(nonzeros(ind(adj)))==max(ind(adj))
                    ind(k) = max(ind(adj));
                elseif i < er
                    ind(k) = attach_lost(ind,adj);
                else
                    lost(end+1) = k;
                end
            end
        end
    end
end

disp('done')

for i = 1:length(lost)
    ind(lost(i)) = attach_lost(ind,find_adjacent(lost(i),XYZ,adjacency,[]));
end

end




function ind = voxel_test(index,XYZ,adjacency,min_adj,criterion)

persistent runN

if isempty(runN)
    runN = 1;
else
    runN = runN + 1;
end

ind = zeros(1,size(XYZ,2));

for k=1:length(ind)

    if ~rem(k,3000)
        disp(['Erosion #' num2str(runN) ' is ' int2str((k/length(ind))*100) '% complete'])
    end

    adj = find_adjacent(k,XYZ,adjacency,index);    
    
    if sum(adj) < min_adj + 1 || sum(adj) < criterion + 1
        continue
    end

    if any(ind(adj))
        ind(k) = min(nonzeros(ind(adj)));
        if min(nonzeros(ind(adj)))~=max(ind(adj))
            cls = unique(nonzeros(ind(adj)));
            for i = 1:length(cls)
                ind(ind==cls(i)) = min(cls);
            end
        end
    else
        ind(k) = max(ind) + 1;
    end
end

end


function [adj] = find_adjacent(k,XYZ,adjacency,index)

if isempty(index)
    index = true(1,size(XYZ,2));
end

adjacent=zeros(size(XYZ));
for dim = 1:3
    a = ( (XYZ(dim,k)==XYZ(dim,:)+1) | (XYZ(dim,k)==XYZ(dim,:)-1) ) & index;
    b = (XYZ(dim,k)==XYZ(dim,:)) & index;
    if ~isempty(a)
        adjacent(dim,a)=1;
    end
    if ~isempty(b)
        adjacent(dim,b)=2;
    end
end

% slower: but agrees in faces...
% tic
% tor_adj = [XYZ(1, :) - XYZ(1,k); XYZ(2, :) - XYZ(2,k); XYZ(3, :) - XYZ(3,k)];
% tor_adj = abs(tor_adj);
% tor_adj(tor_adj > 1) = Inf;  % not actually adjacent
% tor_adj = sum(tor_adj);
% tor_adj(~index) = Inf;
% toc

% 1 = share face, 2 = share edge, 3 = share corner

% in "adjacent" var:
% adjacent faces differ by only one coord
% adjacent edges differ by 1 or 2 coords
% adjacent corners differ by 1-3 coords

% these methods also include "self" voxel
if strcmp(adjacency,'corner')
    adj = all(adjacent);
elseif strcmp(adjacency,'edge')
    adj = all(adjacent) & any(adjacent - 1);
elseif strcmp(adjacency,'surface')
    adj = all(adjacent) & sum(adjacent) > 4;  % tor fixed bug 12/11/13
end

end


function [lind] = attach_lost(ind,adj)

cls=unique(nonzeros(ind(adj)));
for m=1:length(cls)
    num(m)=sum((ind(adj))==cls(m));
end
if sum(num==max(num))==1
    lind=cls(num==max(num));
else
    a=find(num==max(num));
    for m=1:length(a)
        l(m)=sum(ind==cls(a(m)));
    end
    if sum(l==max(l))==1
%         warning(['Voxel at ' num2str(XYZ(1,lost(i))) ',' num2str(XYZ(2,lost(i))) ',' num2str(XYZ(3,lost(i))) ' shares an equal number of adjacent voxels with more than one seperate cluster. It is being associated with the largest of them.']);
        lind=cls(a(l==max(l)));
    else
        warning(['Voxel shares an equal number of adjacent voxels with more than one seperate cluster of identical size. It is being randomly associated with one of them.']);
        r=randperm(length(a));
        lind=cls(a(r(1)));
    end
end

end