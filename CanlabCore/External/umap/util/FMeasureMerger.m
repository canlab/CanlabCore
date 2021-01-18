classdef FMeasureMerger < handle
    properties(Constant)
        USE_ORIGINAL_IDEA=false;
    end
    
    methods(Static)
        function [bestIds, fMeasure]=GetBest(matchIdPerRow, ...
                mergeIdPerRow, matchId, mergeIds, adaptiveBins, transpose)
            if isempty(adaptiveBins)
                [matchSize, fMeasures4MergeIds,mergeSizes,truePositives]...
                    =FMeasureMerger.Compute(matchIdPerRow, ...
                    mergeIdPerRow, matchId, mergeIds, transpose);
                if FMeasureMerger.USE_ORIGINAL_IDEA
                    [bestIds, fMeasure]=...
                        FMeasureMerger.FindBestUsingOriginalIdea(...
                        matchSize, mergeIds, fMeasures4MergeIds, mergeSizes, ...
                        truePositives);
                else
                    [bestIds, fMeasure]=FMeasureMerger.FindBest(...
                    matchSize, mergeIds, fMeasures4MergeIds, mergeSizes, ...
                    truePositives);
                end
            else
                [matchChoices, fMeasures4MergeIds, allMergeChoices] ...
                    =FMeasureMerger.ComputeBins(matchIdPerRow, ...
                    mergeIdPerRow, matchId, mergeIds, adaptiveBins, transpose);
                [bestIds, fMeasure]=FMeasureMerger.FindBestBins(...
                    mergeIds, fMeasures4MergeIds, matchChoices, ...
                    allMergeChoices, adaptiveBins, transpose);
            end
        end
        

        function [matchChoices, fMeasures4MergeIds, allMergeChoices]=...
                ComputeBins(idPerRow1, idPerRow2, ...
                matchId, mergeIds, adaptiveBins, transpose)
            if ~transpose
                matchChoices=MatBasics.LookForIds(idPerRow1, matchId);
            else
                matchChoices=MatBasics.LookForIds(idPerRow2, matchId);
            end
            nMergeIds=length(mergeIds);
            allMergeChoices=cell(1, nMergeIds);
            fMeasures4MergeIds=zeros(1,nMergeIds);
            for i=1:nMergeIds
                if ~transpose
                    mergeChoices=MatBasics.LookForIds(idPerRow2, mergeIds(i));
                    fMeasures4MergeIds(i)=adaptiveBins.fMeasure(matchChoices, mergeChoices);
                else
                    mergeChoices=MatBasics.LookForIds(idPerRow1, mergeIds(i));
                    fMeasures4MergeIds(i)=adaptiveBins.fMeasure(mergeChoices, matchChoices);
                end
                allMergeChoices{i}=mergeChoices;
                if isnan(fMeasures4MergeIds(i))
                    fMeasures4MergeIds(i)=0;
                end
            end
        end

        function [bestIds, fMeasure]=FindBestBins(mergeIds,...
                fMeasures4MergeIds, matchChoices, allMergeChoices, ...
                adaptiveBins, transpose)
            [maxF, idx]=max(fMeasures4MergeIds);
            bestIds=mergeIds(idx);
            fMeasure=-100;% force at least ONE merger 
            nMergeCandidates=length(fMeasures4MergeIds);
            nextFms=zeros(1, nMergeCandidates);
            mergeChoices=allMergeChoices{idx};
            mergeIds(idx)=nan;
            while ~all(isnan(mergeIds))
                for i=1:nMergeCandidates
                    if isnan(mergeIds(i))
                        nextFms(i)=-1;
                    else
                        if ~transpose
                            nextFms(i)=adaptiveBins.fMeasure(...
                                matchChoices, mergeChoices|allMergeChoices{i});
                        else
                            nextFms(i)=adaptiveBins.fMeasure(...
                                mergeChoices|allMergeChoices{i}, matchChoices);
                        end
                    end
                end
                [f, idx]=max(nextFms);
                if f>fMeasure
                    fMeasure=f;
                    mergeChoices=mergeChoices|allMergeChoices{idx};
                    bestIds(end+1)=mergeIds(idx);
                else
                    break; % nothing changed!
                end
                mergeIds(idx)=nan;
            end
            if QfHiDM.F_MEASURE_MERGE_FAST==-1
                FMeasureMerger.NoteIfNoMergerImprovement(fMeasure, maxF, ...
                    bestIds(1));
            end
            bestIds=sort(bestIds);
        end
        
        function [matchSize, fMeasures4MergeIds, mergeSizes, truePositives]=...
                Compute(idPerRow1, idPerRow2, matchId, mergeIds, transpose)
            if ~transpose
                matchChoices=MatBasics.LookForIds(idPerRow1, matchId);
            else
                matchChoices=MatBasics.LookForIds(idPerRow2, matchId);
            end
            matchSize=sum(matchChoices);
            nMergeIds=length(mergeIds);
            fMeasures4MergeIds=zeros(1,nMergeIds);
            truePositives=zeros(1,nMergeIds);
            mergeSizes=zeros(1,nMergeIds);
            for i=1:nMergeIds
                if ~transpose
                    mergeChoices=MatBasics.LookForIds(idPerRow2, mergeIds(i));
                else
                    mergeChoices=MatBasics.LookForIds(idPerRow1, mergeIds(i));
                end
                mergeSizes(i)=sum(mergeChoices);
                truePositives(i)=sum(matchChoices&mergeChoices);
                p=truePositives(i)/mergeSizes(i);
                r=truePositives(i)/matchSize;
                fMeasures4MergeIds(i)=(2*p*r)/(p+r);
                if isnan(fMeasures4MergeIds(i))
                    fMeasures4MergeIds(i)=0;
                end
            end
        end
        
        
        
%This method was necessary over the
%FindBestUsingOriginalIdea when it was found that
%the next highest F-measure does not ALWAYS 
%lead to the best F-measure after merging.
%This was found with the Neutrophil scenario 
%when trying to find first best merger of 2 after 
%starting with a merge candidate with max F-measure 
%that is based on a precision of 212/584 and matchSize of 300
%and the next 2 candidates in decreasing F-measure 
%order is precision 8/105 and 2/54.
%The best F-measure when combining the 2nd candidate
%is 2/54 and not 8/105

        function [bestIds, fMeasure]=FindBest(matchSize, mergeIds, ...
            fMeasures4MergeIds, mergeSizes, truePositives)
            [maxF, idx]=max(fMeasures4MergeIds);
            bestIds=mergeIds(idx);
            fMeasure=-100;% force at least ONE merger 
            truePositive=truePositives(idx);
            mergeSize=mergeSizes(idx);
            truePositives(idx)=nan;
            while ~all(isnan(truePositives))
                tp=truePositive+truePositives;
                p=tp./(mergeSizes+mergeSize);
                r=tp./matchSize;
                [f, idx]=max((2*p.*r)./(p+r));
                if f>fMeasure
                    fMeasure=f;
                    truePositive=truePositive+truePositives(idx);
                    mergeSize=mergeSize+mergeSizes(idx);
                    bestIds(end+1)=mergeIds(idx);
                else
                    break; % nothing changed!
                end
                truePositives(idx)=nan;
            end
            if QfHiDM.F_MEASURE_MERGE_FAST==-1
                FMeasureMerger.NoteIfNoMergerImprovement(fMeasure, maxF,...
                    bestIds(1));
            end
            bestIds=sort(bestIds);
        end
        
        function [bestIds, fMeasure]=FindBestUsingOriginalIdea(matchSize,...
                mergeIds, fMeasures4MergeIds, mergeSizes, truePositives)
            [fmBest2Worst, I]=sort(fMeasures4MergeIds, 'descend');
            idx=I(1);
            bestIds=mergeIds(idx);
            fMeasure=-100;% force at least ONE merger 
            truePositive=truePositives(idx);
            mergeSize=mergeSizes(idx);
            truePositives(idx)=nan;
            tp=truePositive+truePositives;
            p=tp./(mergeSizes+mergeSize);
            r=tp./matchSize;
            [f, idx]=max((2*p.*r)./(p+r));
            if f>fMeasure
                fMeasure=f;
                truePositive=truePositive+truePositives(idx);
                truePositives(idx)=nan;
                mergeSize=mergeSize+mergeSizes(idx);
                bestIds(end+1)=mergeIds(idx);
                nMergeCandidates=length(fMeasures4MergeIds);
                for i=2:nMergeCandidates
                    if fMeasure<.5 && fmBest2Worst(i)<.5
                        %NOT going to get better!
                        break;
                    end
                    idx=I(i);
                    if isnan(truePositives(idx))
                        continue;
                    end
                    c=mergeSizes(idx);
                    tp=truePositive+truePositives(idx);
                    p=tp/(c+mergeSize);
                    r=tp/matchSize;
                    f=(2*p*r)/(p+r);
                    if f>fMeasure
                        fMeasure=f;
                        truePositive=tp;
                        mergeSize=c+mergeSize;
                        bestIds(end+1)=mergeIds(idx);
                    end
                end
            end
            
            if QfHiDM.F_MEASURE_MERGE_FAST==-1
                FMeasureMerger.NoteIfNoMergerImprovement(fMeasure,...
                    fmBest2Worst(1), bestIds(1));
            end
            bestIds=sort(bestIds);
        end
        
        function NoteIfNoMergerImprovement(fMeasure, maxF, bestId)
            if fMeasure<maxF
                fprintf(['No improvement by merging:  Best merger (fm=%s) is LESS than '...
                    'ID=%d (fm=%s)!\n'], ...
                    String.encodeRounded(fMeasure, 3), bestId, ...
                    String.encodeRounded(maxF(1),3));
            end
        end
        
        function [dif, newF, oldF]=DetectChange(matchSize, ...
                 currentTruePositives, currentMergeSize, newTruePositives, newMergeSize)
             p=currentTruePositives/currentMergeSize;
             r=currentTruePositives/matchSize;
             oldF=(2*p*r)/(p+r);
             p=newTruePositives/currentMergeSize;
             r=newTruePositives/matchSize;
             individualF=(2*p*r)/(p+r);
             tp=currentTruePositives+newTruePositives;
             p=tp/(newMergeSize+currentMergeSize);
             r=tp/matchSize;
             newF=(2*p*r)/(p+r);
             dif=newF-oldF;
             if dif<0
                 word='WORSER';
             else
                 word='better';
             end
             fprintf('%s F-measures:  old=%s, individual=%s, merged=%s, change=%s\n',...
                 word, ...
                 String.encodeRounded(oldF, 3), ...
                 String.encodeRounded(individualF, 3), ...
                 String.encodeRounded(newF, 3), ...
                 String.encodeRounded(dif, 3));
        end
    end
end
    