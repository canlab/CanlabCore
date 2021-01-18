
%   Class for  building branches of the dendrogram built by QfTree.m 
%
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
%   QF dissimilarity is described in 
%   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6586874/
%
%   QF Tree (phenogram) is described in 
%   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6586874/

classdef Branch < handle
    properties
        leftName;
        rightName;
    end
    properties(SetAccess=private)
        parentPtr=0;
        ptr;
        leftSize;
        rightSize;
        leftIds;
        leftPtr;
        rightIds;
        rightPtr;
        levelFromRoot;
        levelFromLeaf;
        qfScore;
        H;
        merge;
        key;
    end
    methods
        function this=Branch(tree, ids, qfScore, leftSize, rightSize, ...
                leftIds, rightIds, ptr)
            ptr=int16(ptr);
            this.merge=[leftIds rightIds];
            this.leftSize=leftSize;
            this.rightSize=rightSize;
            this.key=num2str(sort(this.merge));
            this.qfScore=qfScore;
            this.leftPtr=getPtr(leftIds);
            this.rightPtr=getPtr(rightIds);
            this.ptr=ptr;
            this.leftIds=leftIds;
            this.rightIds=rightIds;
            tree.set(this.key, this);
            
            function childPtr=getPtr(childIds)
                if length(childIds)==1
                    childPtr=find(ids==childIds,1);
                else
                    childKey=num2str(sort(childIds));
                    branch_=tree.get(childKey);
                    childPtr=branch_.ptr;
                    branch_.parentPtr=ptr;
                end
            end
        end
        
    end
    
    methods(Static)
        function [phyTree, branchNames, nodeQfs, branchSzs, branchQfs]=...
                PhyTree(branches, numLeaves)
            N=length(branches);
            branchSzs=zeros(1, N);
            branchQfs=zeros(1, N);
            nodeQfs=zeros(1, N+numLeaves);
            top=[];
            for i=1:N
                if branches{i}.parentPtr==0
                    traverse(branches{i}, 1);
                    top(end+1)=i;
                end
                branchSzs(i)=branches{i}.leftSize+branches{i}.rightSize;
            end
            treeLevels=zeros(1,N);
            branchNames=cell(1,N);
            for i=1:N
                br=branches{i};
                treeLevels(i)=br.levelFromRoot;
                branchNames{i}=num2str(br.merge);
                branchQfs(i)=br.qfScore;
                nodeQfs(br.leftPtr)=br.qfScore;
                nodeQfs(br.rightPtr)=br.qfScore;                
            end
            treeLevels=unique(treeLevels);
            mx=max(treeLevels);
            for i=1:N
                branches{i}.levelFromLeaf=(mx-branches{i}.levelFromRoot)+1;
            end
            phyTree=zeros(N,2);
            found=false(1, numLeaves);
            for i=1:N
                phyTree(i, 1)=branches{i}.leftPtr;
                if phyTree(i,1)<=numLeaves
                    found(phyTree(i,1))=true;
                end
                phyTree(i, 2)=branches{i}.rightPtr;
                if phyTree(i,2)<=numLeaves
                    found(phyTree(i,2))=true;
                end
                
            end
            if length(top)>1
                phyTree(end, 1)=numLeaves+top(1);
                phyTree(end, 2)=numLeaves+top(2);
                branchNames{end+1}='root';
            else
                lost=find(~found);
                if length(lost)==1
                    phyTree(end, 1)=numLeaves+top(1);
                    phyTree(end, 2)=lost(1);
                    branchNames{end+1}='root';
                end
            end
            function traverse(branch, total)
                branch.levelFromRoot=total;
                idx=branch.leftPtr-numLeaves;
                if idx>0
                    traverse(branches{idx}, total+1);
                end
                idx=branch.rightPtr-numLeaves;
                if idx>0
                    traverse(branches{idx}, total+1);
                end
            end
        end
    end
end
