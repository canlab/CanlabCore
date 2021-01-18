%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef Clusters < handle
    properties(SetAccess=private)
        ids;
        clues;
        sizes;
        orderIdxs;
        fcsIdxs;
        N;
    end
    properties
        colors;
        selColors;
    end
    
    methods
        function this=Clusters(clues, cluColorFactory)
            ids=unique(clues);
            zeroIdx=find(ids==0,1);
            if ~isempty(zeroIdx)
                ids(zeroIdx)=[];
            end
            sizes=histc(clues, ids);
            orderIdxs=MatBasics.GetSizeRankings(sizes);
            N=length(ids);
            if nargin>1
                if isempty(zeroIdx)
                    colors=zeros(N, 3);
                    selColors=zeros(N, 3);
                    ii=0;
                else
                    colors=zeros(N+1, 3);
                    selColors=zeros(N+1, 3);
                    colors(1,:)=[0 0 0];
                    selColors(1,:)=[0 0 0];
                    ii=1;
                end
                for i=1:N
                    colors(i+ii, :)=cluColorFactory.get(orderIdxs(i), false);
                    selColors(i+ii, :)=cluColorFactory.get(orderIdxs(i), true);
                end
                this.colors=colors;
                this.selColors=selColors;
            end
            fcsIdxs=cell(1, N);
            for i=1:N
                id=ids(i);
                fcsIdxs{i}=find(clues==id);
            end
            this.fcsIdxs=fcsIdxs;
            this.ids=ids;
            this.clues=clues;
            this.sizes=sizes;
            this.orderIdxs=orderIdxs;
            this.N=N;
        end
        
        function f=fMeasure(this, thisI, that, thatI)
            same=intersect(this.fcsIdxs{thisI}, that.fcsIdxs{thatI});
            f=Clusters.F_measure(length(same), this.sizes(thisI), ...
                that.sizes(thatI));
        end
    end
    
    methods(Static)
        function f=F_measure(nTruePos, sizeTrainingSet, sizeTestSet)
            precision=nTruePos/sizeTestSet;
            recall=nTruePos/sizeTrainingSet;
            f=(2*precision*recall)/(precision+recall);
            if isnan(f)
                f=0;
            end
        end
       
         function f=F_measures(nsTruePos, sizeTrainingSet, sizesTestSet)
            precision=nsTruePos./sizesTestSet;
            recall=nsTruePos./sizeTrainingSet;
            f=(2*precision.*recall)./(precision+recall);
            if any(isnan(f))
                f(isnan(f))=0;
            end
         end
         
         function p=Precision(nTruePos, sizesTestSet)
            p=nTruePos./sizesTestSet;
            if any(isnan(p))
                p(isnan(p))=0;
            end
         end
         
         function r=Recall(nTruePos, sizeTrainingSet)
            r=nTruePos./sizeTrainingSet;
            if any(isnan(r))
                r(isnan(r))=0;
            end
         end
         
         function [clrs, selClrs]=GrayColors(clustSizes)
             numClusts=length(clustSizes);
             clrs=zeros(numClusts+1,3);
             selClrs=zeros(numClusts+2,3);
             [~,I]=sort(clustSizes, 'descend');
             grayGap=.8/numClusts;
             for idx_=2:numClusts
                 cluIdx=I(idx_);
                 %cluIdx=idx_;
                 clr=.04+(idx_*grayGap);
                 rank=mod(idx_, 3);
                 if rank==0
                     greenFactor=1;
                     blueFactor=1.1;
                 elseif rank==1
                     greenFactor=1.15;
                     blueFactor=1;
                 else
                     blueFactor=.98;
                     greenFactor=1.03;
                 end
                 clrs(cluIdx,:)=[clr greenFactor*clr clr ];
                 selClrs(cluIdx,:)=.9*clrs(cluIdx,:);
                 %clrs(cluIdx,:)=[clr clr*greenFactor blueFactor*clr ];
             end
         end
         
         function F=MeasureF(this, that)
            nThis=length(this.clues);
            F=0;
            for i=1:this.N
                a=zeros(1, that.N);
                for j=1:that.N
                    a(j)=this.fMeasure(i, that, j);
                end
                F=F+(this.sizes(i)/nThis*max(a));
            end
        end
    end
end