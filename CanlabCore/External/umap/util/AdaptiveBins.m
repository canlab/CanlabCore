%   Class for probability binning described at
%   https://www.ncbi.nlm.nih.gov/pubmed/11598946
%
%   Bioinformatics inventors
%   Original:  Roederer M1, Moore W, Treister A, Hardy RR, Herzenberg LA.
%   Incorporation into QF matching:  Darya Orlova <dyorlova@gmail.com>
%
%   Software Developers: Connor Meehan <connor.gw.meehan@gmail.com>
%           Stephen Meehan <swmeehan@stanford.edu> 
%
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

classdef AdaptiveBins<handle
    properties(SetAccess=private)
        means;
        dists;
        teachPtrs;
        studPtrs;
        teachSize;
        studSize;
        maxDist=0;
        binSize;
    end 
    
    properties(Constant)
        MAX_SIZE=25000;
        MIN_SIZE=30;
    end
    
    methods
        function this=AdaptiveBins(teachData, studData, minBinSize, ...
                sizeIsFixedSplit, teachStudCacheFile, studTeachCacheFile)
            if nargin<6
                studTeachCacheFile=[];
                if nargin<5
                    teachStudCacheFile=[];
                    if nargin<4
                        sizeIsFixedSplit=false;
                        if nargin<3
                            minBinSize=[];
                        end
                    end
                end
            end
            if ~isempty(teachStudCacheFile)
                teachStudCacheFile=[teachStudCacheFile 'pb_v3_' num2str(minBinSize) ...
                    '_' num2str(sizeIsFixedSplit) '.mat'];
            end
            if ~isempty(studTeachCacheFile)
                studTeachCacheFile=[studTeachCacheFile 'pb_v3_' num2str(minBinSize) ...
                        '_' num2str(sizeIsFixedSplit) '.mat'];
            end
            if exist(teachStudCacheFile, 'file')
                load(teachStudCacheFile, 'means', 'teachPtrs', 'studPtrs', 'binSize'); 
                this.means=means;
                this.teachPtrs=teachPtrs;
                this.studPtrs=studPtrs;
                this.binSize=binSize;
            elseif exist(studTeachCacheFile, 'file')
                load(studTeachCacheFile, 'means', 'teachPtrs', 'studPtrs', 'binSize'); 
                this.means=means;
                this.teachPtrs=studPtrs;
                this.studPtrs=teachPtrs;
                this.binSize=binSize;
            else
                [this.means, this.teachPtrs, this.studPtrs, testW1, ...
                    testW2, this.binSize]=AdaptiveBins.Create(teachData,...
                    studData, minBinSize, sizeIsFixedSplit, ...
                    teachStudCacheFile);
                assert(round(sum(testW1),11)==1);
                assert(round(sum(testW2),11)==1);
                %disp('AdaptiveBins.Create CALLED!');
            end
            this.teachSize=size(teachData, 1);
            this.studSize=size(studData, 1);
            if size(this.means, 1)<AdaptiveBins.MAX_SIZE
                this.dists=MatBasics.PDist2Self(this.means);
                this.maxDist=max(max(this.dists));
            end
        end
        
        function [f,p,r]=fMeasure(this, teachChoices, studChoices)
            hPtrs=this.teachPtrs(teachChoices);
            fPtrs=this.studPtrs(studChoices);
            u=unique([hPtrs fPtrs]);
            teachCnts=histc(hPtrs, u);
            studCnts=histc(fPtrs, u);
            [f, p, r]=AdaptiveBins.F_measure(teachCnts,studCnts);
        end
        
        function e=emd(this, teachChoices, studChoices)
            hPtrs=this.teachPtrs(teachChoices);
            fPtrs=this.studPtrs(studChoices);
            u=unique([hPtrs fPtrs]);
            tBinCnts=histc(hPtrs, u);
            sBinCnts=histc(fPtrs, u);
            tSize=this.teachSize;%=sum(teachChoices);
            sSize=this.studSize;%=sum(studChoices);
            tW=tBinCnts/tSize;
            sW=sBinCnts/sSize;
            M=this.means(u, :);
            [~,e]=emd_flow(M, M, tW, sW, @gdf);
        end
        
        function [meansOrDistances, teachWeights, studWeights]=...
                weigh(this, teachChoices, studChoices, weighBySampleSize)
            if weighBySampleSize
                tSize=this.teachSize;%=sum(teachChoices);
                sSize=this.studSize;%=sum(studChoices);
            else
                tSize=sum(teachChoices);
                sSize=sum(studChoices);
            end
            hPtrs=this.teachPtrs(teachChoices);
            fPtrs=this.studPtrs(studChoices);
            u=unique([hPtrs fPtrs]);
            tBinCnts=histc(hPtrs, u);
            sBinCnts=histc(fPtrs, u);
            if AdaptiveBins.DEBUG_LEVEL>0
                f=AdaptiveBins.F_measure(tBinCnts, sBinCnts);
            end
            teachWeights=tBinCnts/tSize;
            studWeights=sBinCnts/sSize;
            if QfHiDM.DEBUG_LEVEL>0
                if ~weighBySampleSize
                    assert(round(sum(teachWeights),11)==1);
                    assert(round(sum(studWeights),11)==1);
                end
            end
            if isempty(this.dists)
                meansOrDistances=this.means(u, :);
            else
                meansOrDistances=this.dists(u, u);
                if QfHiDM.DEBUG_LEVEL>0
                    test=MatBasics.PDist2Self(this.means(u, :));
                    if ~isequal(round(test,8), round(meansOrDistances,8))
                        disp('floating point insensitity at 8th digit?...');
                    end
                end
            end
        end
        
        
    end
    
    methods(Static)
        function [f, precision, recall]=F_measure(teachBinCnts, studBinCnts)
            tSize=sum(teachBinCnts);
            sSize=sum(studBinCnts);
            teachBins=studBinCnts>teachBinCnts&teachBinCnts>0;
            studBins=studBinCnts<=teachBinCnts&studBinCnts>0;
            truePositive=sum(teachBinCnts(teachBins))+sum(studBinCnts(studBins));
            
            % this statement on 5 lines below does nothing much to help
            %truePositive = floor( truePositive + ...
            %    sum((studBinCnts(teachBins)-teachBinCnts(teachBins))...
            %    .* (studBinCnts(teachBins)/sSize)) +...
            %  sum((teachBinCnts(studBins)-studBinCnts(studBins))...
            %    .*  (teachBinCnts(studBins)/tSize)) );
            precision=truePositive/sSize;
            recall=truePositive/tSize;
            f=(2*precision*recall)/(precision+recall);
            if isnan(f)
                f=0;
            end
        end
        
        function [f, p, r ]=F_measureOld(teachCnts, studCnts)
            tSize=sum(teachCnts);
            sSize=sum(studCnts);
            cnt=tSize+sSize;
            tRatio=tSize/cnt;
            sRatio=sSize/cnt;
            isOverlap=teachCnts>0 & studCnts>0;
            truePositive=sum(studCnts(isOverlap));
            recalled=sum(teachCnts(isOverlap));
            percentOverlap=(truePositive+recalled)/(sSize+tSize);
            tRecallRatio=recalled/tSize;
            sRecallRatio=truePositive/sSize;
            ratio=(tRecallRatio*tRatio)+(sRecallRatio*sRatio);
            mnCnt=min([tSize sSize]);
            finalTruePositive=ratio*mnCnt;
            f=Clusters.F_measure(finalTruePositive, tSize, sSize);
            p=finalTruePositive/sSize;
            r=finalTruePositive/tSize;
        end
        
        function [means, teachPtrs, studPtrs, teachWeights, studWeights, ...
                minBinSize]=Create(teachData, studData, minBinSize, ...
                doFixedSplits, cacheFile, trackRow)
            if nargin < 4
                doFixedSplits=false;
            end
            doNotDuplicate=isequal(studData, teachData);
            if doNotDuplicate
                data=teachData;
                studData=[];
                [N, m] = size(data);
                minN=N;
            else
                data=[teachData;studData];
                [N, m] = size(data);
                minN=min([size(teachData, 1) size(studData,1)]);
            end
            currentBins=cell(1,1);
            firstEntry=cell(1);
            firstEntry{1}=data; %first set is 1 bin-->the original data
            currentBins{1}=firstEntry;
            currentBinPtrs{1}={(1:N)'};
            if nargin<3 || isempty(minBinSize)
                %minBinSize=floor(2*log(N));
                minBinSize=floor(2*log(minN));
            else
                minBinSize=floor(minBinSize);
            end
            %fprintf('The minimum probability bin size is %d\n', minBinSize);
            i=0;
            while true
                i=i+1;
                if doFixedSplits && i >minBinSize
                    break;
                end
                lastBins=currentBins{i};
                lastBinPtrs=currentBinPtrs{i};
                eventsPerBin=size(lastBins{1}, 1);
                if eventsPerBin==1
                    break;
                end
                if eventsPerBin/2<minBinSize
                    if ~doFixedSplits
                        break;
                    end
                end
                nextBins=cell(1,2^i);
                nextBinPtrs=cell(1,2^i);
                binIdx=1;
                numOfBins=length(lastBins);
                for j=1:numOfBins
                    bin=lastBins{j};
                    binPtrs=lastBinPtrs{j};
                    if size(bin,1) == 1
                        nextBins{binIdx} = bin;
                        nextBinPtrs{binIdx}=binPtrs;
                    else
                        variance=var(bin,1);
                        [~,maxCol]=max(variance);
                        [~,ptrs] = sort(bin(:,maxCol));
                        splitPtr = ceil(length(ptrs)/2);
                        %This choice of the splitting_index gives the bin with smaller numbers
                        %the extra data point, if there is one. I made this choice arbitrarily.
                        nextBins{binIdx}=bin(ptrs(1:splitPtr),:);
                        nextBinPtrs{binIdx}=binPtrs(ptrs(1:splitPtr),:);
                        binIdx=binIdx+1;
                        nextBins{binIdx}=bin(ptrs(splitPtr + 1:length(ptrs)),:);
                        nextBinPtrs{binIdx}=binPtrs(ptrs(splitPtr+1:length(ptrs)),:);
                    end
                    binIdx=binIdx+1;
                end
                currentBins{i + 1} = nextBins(1:binIdx-1);
                currentBinPtrs{i + 1} = nextBinPtrs(1:binIdx-1); 
            end
            N1 = size(teachData,1);
            teachPtrs=zeros(1, N1);
            finalBins = currentBins{end};
            finalBinPtrs = currentBinPtrs{end};
            numBins=length(finalBins);
            means=zeros(numBins,m);
            if nargout>3
                N2 = size(studData,1);
                studPtrs=zeros(1, N2);
                teachWeights = zeros(numBins,1);
                studWeights = zeros(numBins,1);
                for i=1:numBins
                    bins=finalBins{i};
                    means(i,:)=mean(bins,1);
                    binPtrs=finalBinPtrs{i};
                    p=binPtrs<=N1;
                    teachPtrs(binPtrs(p))=i;
                    studPtrs(binPtrs(~p)-N1)=i;
                    teachWeights(i)=sum(p)/N1;
                    studWeights(i)=sum(~p)/N2;
                end
            else
                if doNotDuplicate
                    for i=1:numBins
                        bins=finalBins{i};
                        means(i,:)=mean(bins,1);
                        binPtrs=finalBinPtrs{i};
                        p=binPtrs<=N1;
                        teachPtrs(binPtrs(p))=i;
                    end
                else
                    N2 = size(studData,1);
                    studPtrs=zeros(1, N2);
                    for i=1:numBins
                        bins=finalBins{i};
                        means(i,:)=mean(bins,1);
                        binPtrs=finalBinPtrs{i};
                        p=binPtrs<=N1;
                        teachPtrs(binPtrs(p))=i;
                        studPtrs(binPtrs(~p)-N1)=i;
                    end
                end
            end
           if AdaptiveBins.DEBUG_LEVEL>0
               teachCnt=floor(mean(teachWeights(teachWeights>0))*N1);
               studCnt=floor(mean(studWeights(studWeights>0))*N2);
               fprintf('this=%d/%d %s%% and stud=%d/%d %s%%\n', ...
                   teachCnt, N1, num2str(round(teachCnt/N1*100, 3)),...
                   studCnt, N2, num2str(round(studCnt/N2*100, 3)));
           end
           if doNotDuplicate
               studPtrs=teachPtrs;
               if nargout>3
                   studWeights=teachWeights;
               end
           end
           if nargin>4 && ~isempty(cacheFile)
               binSize=minBinSize;
               save(cacheFile, 'means', 'teachPtrs', 'studPtrs', 'binSize'); 
           end
        end
       
        function lvl=DEBUG_LEVEL
            lvl=0;
            %lvl=1;
        end
        
        function emd=Emd(tData, sData, bins)
            try
                tN=size(tData, 1);
                sN=size(sData, 1);
                mx=max([tN sN]);
                minEvents=floor(2*log(mx));
                numberOfBins=floor(mx/minEvents);
                if numberOfBins<32
                    bins=5;
                elseif numberOfBins<64
                    bins=6;
                elseif numberOfBins<128
                    bins=7;
                elseif nargin<3
                    bins=8;
                end
                [tM, ~, ~, tW]=AdaptiveBins.Create(tData, tData, bins, true);
                [sM, ~, ~, sW]=AdaptiveBins.Create(sData, sData, bins, true);
                [~, emd]=emd_flow(tM, sM, tW, sW, @gdf);
            catch
                emd=nan;
            end
        end
        
         function bins=MeansWeightsPtrs(data)
            MIN_BINS=8192;
            MIN_EVENTS_PER_BIN=5;
            MAX_EVENTS_PER_BIN=34;
            N=size(data, 1);
            eventsPerBin=floor(2*log(N));
            numberOfBins=floor(N/eventsPerBin);
            if numberOfBins<MIN_BINS
                eventsPerBin=floor(N/MIN_BINS);
            end
            if eventsPerBin<MIN_EVENTS_PER_BIN
                eventsPerBin=MIN_EVENTS_PER_BIN;
            elseif eventsPerBin>MAX_EVENTS_PER_BIN
                eventsPerBin=MAX_EVENTS_PER_BIN;
            end
            if numberOfBins>2^14 %16384
                eventsPerBin=MAX_EVENTS_PER_BIN;
            end
            [bins.means, bins.ptrs, ~, bins.weights]=...
                AdaptiveBins.Create(data, data, eventsPerBin, false);
            %binObj.ptrs=AdaptiveBins.Fit(data, binObj, binObj.ptrs);
         end
        
         function bins=Fit(data, binObj, ptrs)
             space='CityBlock';
             [dists, bins]=pdist2(binObj.means, data, space, ...
                 'Smallest', 1);
             if nargin>2
                 difs=abs(bins-ptrs);
                 lgcl=difs>0;                 
                 R=size(binObj.means,1);
                 if R>100
                     mn=500;
                 else
                     mn=R/2;
                 end
                 [dists4, bins4]=pdist2(binObj.means, data(lgcl,:), ...
                     space, 'Smallest', mn);%top 100
                 N=sum(lgcl);
                 idxs=find(lgcl);
                 trueBinDists=zeros(N, 2);
                 for i=1:N
                     idx=idxs(i);
                     d=dists4(:, i);
                     b=bins4(:, i);
                     ptr=ptrs(idx);
                     D=pdist2(data(idx,:), binObj.means(ptr,:), space);
                     idx2=find(b==ptr, 1);
                     if ~isempty(idx2)
                         trueBinDists(i,1)=idx2;
                     end
                     trueBinDists(i,2)=D;
                 end
             end
         end
    end
    
end
