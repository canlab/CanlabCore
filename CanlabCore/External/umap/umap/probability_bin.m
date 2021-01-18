%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Funded by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
function [means, teachPtrs, studPtrs, teachWeights, studWeights, ...
    minBinSize]=probability_bin(teachData, studData, minBinSize, ...
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
    currentBinPtrs{1}={[1:N]'};
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
       
        
