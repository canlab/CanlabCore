%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%


classdef MatBasics
    methods(Static)
        function order=Sort1D(nums, canon)
            order=[];
            N=length(canon);
            for i=1:N
                idx=find(canon(i)==nums, 1);
                if ~isempty(idx)
                    order(end+1)=nums(idx);
                    nums(idx)=[];
                end
            end
            if ~isempty(nums)
                order=[order nums];
            end
        end
        
        function data=StretchXyData(data, perc)
            if size(data,1)>2
                mn=min(data);
                mx=max(data);
                r=range(data);
                s=r*perc;
                data=[data; mn(1)-s(1) mn(2)-s(2); mx(1)+s(1) mx(2)+s(2)];
            end
        end
        
        function cnt=CountNonZeroPerColumn(ids)
            N=size(ids,2);
            cnt=0;
            for i=1:N
                cnt=cnt+sum(ids(:, i)>0);
            end
        end

        function idxs=IndexOf(v1, v2)
            idxs=[];
            N=length(v1);
            for i=1:N
                num=v1(i);
                idx=find(v2==num,1);
                if ~isempty(idx)
                    idxs(end+1)=idx;
                end
            end
        end
        function rank=RankForStdDev(data, topRank)
            avg=mean(data);
            dev=std(data);
            devs=abs(data-avg)/dev;
            l=devs<=1;
            d=range(data(l))/topRank;
            mn=min(data(l));
            rank=zeros(1, length(data))+topRank;
            rank(l)=floor((data(l)-mn)/d)+1;
        end
        % Produces distance. 
        %   size of dist is Nx1 where N=size(matrix1,1)=size(matrix2,1)
        %   size(matrix1,2) must equal size(matrix2, 2)
        %   distance is cityblock, cosine or euclidean.  
        % EQUIVALENT (but faster) to 
        %    for i=1:size(matrix1,1)
        %       pdist2(matrix1(i,:), matrix2(i,:), distance)
        %    end
        function dist=Distance(matrix1, matrix2, distance)
            if strcmpi(distance, 'cityblock')
                dist=MatBasics.CityBlock(matrix1, matrix2);
            elseif strcmpi(distance, 'cosine')
                dist=MatBasics.Cosine(matrix1, matrix2);
            elseif strcmpi(distance, 'euclidean')
                dist=MatBasics.Euclidean(matrix1, matrix2);
            else
                assert(false, [ distance ' is not a supported distance ']);
            end
        end

        % Produces euclidean distance. 
        %   size of dist is Nx1 where N=size(matrix1,1)=size(matrix2,1)
        %   size(matrix1,2) must equal size(matrix2, 2)
        % EQUIVALENT (but faster) to 
        %    for i=1:size(matrix1,1)
        %       pdist2(matrix1(i,:), matrix2(i,:), 'euclidean')
        %    end
        function dist=Euclidean(matrix1, matrix2)
            dist=sqrt(sum((matrix1-matrix2).^2, 2));
        end

        % Produces cityblock distance. 
        %   size of dist is Nx1 where N=size(matrix1,1)=size(matrix2,1)
        %   size(matrix1,2) must equal size(matrix2, 2)
        % EQUIVALENT (but faster) to 
        %   for i=1:size(matrix1,1)
        %       pdist2(matrix1(i,:), matrix2(i,:), 'cityblock')
        %    end
        
        function dist=CityBlock(r1, r2)
            dist=sum(abs(r1-r2), 2);
        end

        % Produces cosine distance. 
        %   size of dist is Nx1 where N=size(matrix1,1)=size(matrix2,1)
        %   size(matrix1,2) must equal size(matrix2, 2)
        % EQUIVALENT (but faster) to 
        %   for i=1:size(matrix1,1)
        %       pdist2(matrix1(i,:), matrix2(i,:), 'cosine')
        %    end
        
        function dist=Cosine(matrix1, matrix2)
            numerator=sum(matrix1.*matrix2, 2);
            denominator=sqrt(sum(matrix1.^2,2)).*sqrt(sum(matrix2.^2, 2));
            dist=1-(numerator./denominator);
            
        
        end
        
        function dist=DevDist(avg, dev)
            dist=pdist2(avg, avg+dev);
        end
        
        function xy=FlipY(grid, M)
            [X, Y]=ind2sub([M M], grid);
            M1=M+1;
            Y=M1-Y;
            xy=sub2ind([M M], X, Y);
        end
        function sdus=SduDist(self,other)
            N=size(self,2);
            if N~=size(other,2)
                sdus=[];
                return;
            end
            selfM=mean(self);
            selfS=std(self);
            otherM=mean(other);
            otherS=std(other);
            sdus=zeros(1,N);
            for i=1:N
                dist=abs(selfM(i)-otherM(i));
                if selfS(i)<otherS(i)
                    sdus(i)=dist/selfS(i);
                else
                    sdus(i)=dist/otherS(i);
                end
            end
        end
        
        function sdus=SduDist2(self,otherM, otherS)
            N=size(self,2);
            if N~=length(otherM)
                sdus=[];
                return;
            end
            selfM=mean(self);
            selfS=std(self);
            sdus=zeros(1,N);
            for i=1:N
                dist=abs(selfM(i)-otherM(i));
                if selfS(i)<otherS(i)
                    sdus(i)=dist/selfS(i);
                else
                    sdus(i)=dist/otherS(i);
                end
            end
        end
        
        function out=Trim(in, less)
            if nargin<2
                less=1;
            end
            out=[];
            R=size(in);
            for r=1:R-less
                out(end+1,:)=in(r,:);
            end
        end
        function charIds=ToStrings(ids)
            N=length(ids);
            charIds=cell(1, N);
            for i=1:N
                charIds{i}=num2str(ids(i));
            end            
        end
        
        function [mnRows, mnCols]=Min(M)
            [R, C]=size(M);
            if R==1
                mnCols=ones(1,C);
            else
                [~,mnCols]=min(M);
            end
            if C==1
                mnRows=ones(1, R);
            else
                [~,mnRows]=min(M');
            end
        end
        function D=Dist(arg1, arg2)
            [R1,C]=size(arg1);
            [R2, C2]=size(arg2);
            D=zeros(R1, R2);
            if C ~= C2
                error('arg1 has %d columns and arg2 has %d', C, C2);
            end
            for r1=1:R1
                left=arg1(r1, :);
                for r2=1:R2
                    right=arg2(r2,:);
                    total=0;
                    for c=1:C
                        diff=right(c)-left(c );
                        total=total+(diff*diff);
                    end
                    D(r1, r2)=sqrt(total);
                end
            end
        end
        function b=UniqueSpacedBy(mtrx, spacing)
            if nargin<2
                spacing=1000;
            end
            u=unique(mtrx);
            [~,b]=ismember(mtrx, u);
            b=(b-1)*spacing;
        end

        function b=UniqueOnLogScale(mtrx)
            u=unique(mtrx);
            u=[0 u'];
            [~,b]=ismember(mtrx, u);
            b=(b-1).^(1+(.22*b));
        end

        function tmr=DoOften(fnc, secs, startRightAway)
            tmr=timer;
            tmr.StartDelay=secs;
            tmr.Period=1;
            tmr.ExecutionMode='fixedRate';
            tmr.TimerFcn=fnc;
            if startRightAway
                start(tmr);
            end
        end
        
        function DoLater(fnc, delaySecs)
            tmr=timer;
            tmr.StartDelay=delaySecs;
            tmr.TimerFcn=fnc;
            start(tmr);
        end
        
        function point=getPercentPoint(data, percent)
            if isempty(data)
                point=[];
                return;
            end
            data=sort(data);
            N=length(data);
            P=floor(N*percent);
            if P<1
                P=1;
            end
            point=data(P);
        end
        
        function [pairs, unpaired]=FindBestPairs(matrix, removeIdentity)
            if nargin<2
                removeIdentity=0;
                if nargin<1
                    %TESTING
                    matrix=[0 1 2 3 .5;1 0 1 2 .5; 2 1 0 4 .5; 3 2 4 0 .5;.5 .5 .5 .5 0];
                    matrix=[0 1 2 3 .5;1 0 1 2 3; 2 1 0 4 .5; 3 2 4 0 .5;.5 3 .5 .5 0];
                    matrix=[0 1.1 2 3 .5;1.1 0 1 2 3; 2 1 0 4 .5; 3 2 4 0 .5;.5 3 .5 .5 0];
                    matrix=[0 1.1 2 3 .5;1.1 0 1 2 3; 2 1 0 4 5; 3 2 4 0 .5;.5 3 5 .5 0];
                    removeIdentity=100;
                end
            end
            if removeIdentity>0
                [R,C]=size(matrix);
                if R==C
                    for col=1:R
                        matrix(col,col)=removeIdentity;
                    end
                end
            end
            [minForCol, cntMinForCol, isMinForColMinForRow]=...
                QfHiDM.FindMinRowsForCol(matrix,false);
            N=length(minForCol);
            done=false(1, N);
            pairs=[];
            unpaired=[];
            for col=1:N
                if done(col)
                    continue;
                end
                if ~isMinForColMinForRow(col)
                    unpaired(end+1)=col;
                    continue;
                end
                nums=minForCol{col};
                N2=length(nums);
                if N2==1
                    pairs(end+1,:)=[col nums];
                    done(nums)=true;
                else
                    [mn, I]=min(matrix(nums, col));
                    pairs(end+1,:)=[col nums(I)];
                    done(nums(I))=true;
                    for idx=1:N2
                        if idx~=I
                            unpaired(end+1)=nums(idx);
                        end
                    end
                end
            end
            unpaired=unique(unpaired);
        end
        
        function [u, szs]=GetUniqueIdsAndCntIfGt0(ids)
            u=unique(ids);
            u=u(u>0);
            szs=histc(ids(:), u);
        end
        function found=LookForIds(ids, lookFor)
            [R, C]=size(ids);
            if C>0 && R>0
                found=ismember(ids(:, 1), lookFor);
                for i=2:C
                    found=found|ismember(ids(:,i), lookFor);
                end
            else
                found=[];
            end
        end
        
        function bfm=BestPossibleF_measures(tIds, tLookFor, sIds, sIdSets)
            tCnt=MatBasics.CountIds(tIds, tLookFor);
            sCnts=MatBasics.CountIdSets(sIds, sIdSets);
            sizesRecalled=min(tCnt, sCnts);
            bfm=Clusters.F_measures(sizesRecalled, tCnt, sCnts);
        end

        function cnt=CountIds(ids, lookFor)
            if iscell(lookFor)
                lookFor=lookFor{1};
            end
            [R, C]=size(ids);
            if C>0 && R>0
                cnt=sum(ismember(ids(:, 1), lookFor));
                for i=2:C
                    cnt=cnt+sum(ismember(ids(:,i), lookFor));
                end
            else
                cnt=0;
            end
        end

        function cnts=CountIdSetsSlow(ids, idSets)
            [R, C]=size(ids);
            N=length(idSets);
            cnts=zeros(1, N);
            if C==0 || R==0
                return;
            end
            for i=1:N
                if iscell(idSets)
                    lookFor=idSets{i};
                else
                    lookFor=idSets(i);
                end
                cnts(i)=sum(ismember(ids(:, 1), lookFor));
                for j=2:C
                    cnts(i)=cnts(i)+sum(ismember(ids(:,j), lookFor));
                end                
            end
        end
       
        function cnts=CountIdSets(ids, idSets)
            [R, C]=size(ids);
            if C==0 || R==0
                return;
            end
            u=[];
            for i=1:C
                u=[u unique(ids)];
            end
            u=unique(u);
            nU=length(u);
            ttl=zeros(1, nU);
            for i=1:nU
                id=u(i);
                ttl(i)=sum(sum(ids==id));
            end
            N=length(idSets);
            cnts=zeros(1, N);
            if iscell(idSets)
                for i=1:N
                    cnts(i)=sum(ttl(ismember(u, idSets{i})));
                end
            else
                for i=1:N
                    cnts(i)=sum(ttl(ismember(u, idSets(i))));
                end
            end
        end
        
        function d = PDist2Self(a)
            d=pdist2(a, a);
            R=size(d);
            for i=1:R
                d(i,i)=0;
            end
        end

        function [gb, sGb]=GigaByteMatrix(num)
            gb=num*num*8;
            sGb=String.encodeRounded(gb/(2^30), 2, true);
        end
        
        function [D, d_max, A_IJ]=Qf(distances, h, f)
            d_max=max(max(distances));
            A_IJ=1-distances/d_max;
            [H, F]=meshgrid(h-f, h-f);
            D=sqrt(sum(sum(A_IJ.*H.*F)));
            if isnan(D)
                D=QfHiDM.MAX_QF_DISTANCE;%close to max QF distance
            end
            
        end
        
        function [D, d_max, A_IJ]=QfNonZero(a, h, f)
            l=h~=0 | f~=0;
            if ~any(l)
                D=0;
                return;
            end
            [R, C]=size(a);
            if R==C
                d_max=max(max(a));
            end
            if ~all(l)
                h=h(l);
                f=f(l);
                nonZero=sum(l);
                if R==C
                    if nonZero/R<.12
                        old=a;
                    end
                    a=a(l,:);
                    a=a(:,l);
                else
                    a=a(l, :);
                end
            end
            if R~=C
                a=pdist2(a, a);
                d_max=max(max(a));
            end
            A_IJ=1-a/d_max;
            [H, F]=meshgrid(h-f, h-f);
            D=sqrt(sum(sum(A_IJ.*H.*F)));
            if isnan(D)
                D=QfHiDM.MAX_QF_DISTANCE;%close to max QF distance
            end
        end
        
        %slower than MatLab's pdist2 but works with complex numbers
        function d = PDist2(a,b)
            if nargin<2
                b=a;%self distance
            end
            if (size(a,1) ~= size(b,1))
                error('A and B should be of same dimensionality');
            end
            aa=sum(a.*a,1); bb=sum(b.*b,1); ab=a'*b;
            d = sqrt(abs(repmat(aa',[1 size(bb,2)]) + ...
                repmat(bb,[size(aa,2) 1]) - 2*ab));
            if nargin<2
                [R,C]=size(d);
                for i=1:R
                    d(i,i)=0;
                end
            end
        end
        
        function [devUnits, tFreq, sFreq,tMeans, sMeans]=GetMeanUnits(...
                tData, sData, tIdPerRow,tIdSet, sIdPerRow, sIdSet)
            devUnits=[];
            tFreq=[];
            sFreq=[];
            tMeans=[];
            sMeans=[];
            if isempty(tData)||isempty(sData)
                return;
            end
            tChoices=MatBasics.LookForIds(tIdPerRow, tIdSet);
            sChoices=MatBasics.LookForIds(sIdPerRow, sIdSet);
            if MatBasics.DEBUG('meanunits')
                [~, C]=size(tIdPerRow);
                if C==1
                    assert(isequal(tChoices, ismember(tIdPerRow, tIdSet)) )
                end
                [~, C]=size(sIdPerRow);
                if C==1
                    assert(isequal(ismember(sIdPerRow, sIdSet), sChoices))
                end
            end
            tData_=tData(tChoices, :);
            sData_=sData(sChoices, :); 
            if size(sData_,1)<2
                return;
            end
            sMeans=mean(sData_);
            if size(tData_,1)<2
                return;
            end
            tMeans=mean(tData_);
            tStds=std(tData_);
            devUnits=abs(sMeans-tMeans)./tStds;
            if any(devUnits>4)
                try
                    [tMeans;sMeans;tStds;devUnits]
                    fprintf('%d parameters are over 4 dev unit away\n',...
                        sum(devUnits>4));
                    disp('herp');
                catch ex
                    ex.getReport
                end
            end
            if nargout>1
                tFreq=sum(tChoices)/size(tData,1);
                sFreq=sum(sChoices)/size(sData,1);
            end
        end
        
        function [f,p,r]=F_measure(tData, sData, tIdPerRow,tIdSet, sIdPerRow, ...
                sIdSet, areEqual)
            f=-1;
            if isempty(tData)||isempty(sData)
                return;
            end
            if nargin<7 
                if ~isequal(tData, sData)
                    return;
                end
            else
                if ~areEqual
                    return;
                end
            end
            tChoices=MatBasics.LookForIds(tIdPerRow, tIdSet);
            sChoices=MatBasics.LookForIds(sIdPerRow, sIdSet);
            sizeRecalled=sum(tChoices&sChoices);
            sizeTrainingSet=sum(tChoices);
            sizeTestSet=sum(sChoices);
            f=Clusters.F_measure(sizeRecalled, sizeTrainingSet, sizeTestSet);
            p=Clusters.Precision(sizeRecalled, sizeTestSet);
            r=Clusters.Recall(sizeRecalled, sizeTrainingSet);
        end
        
        
        function mergers=NumCombos(nMerger)
            mergers=1;
            lastComboSize=nMerger-1;
            for comboSize=2:lastComboSize
                mergers=mergers+nchoosek(nMerger, comboSize);
            end
        end

        function [tO, sO]=Overlap(tData, sData, tIdPerRow,tIdSet, ...
                sIdPerRow, sIdSet, areEqual)
            sO=2;
            tO=2;
            if isempty(tData)||isempty(sData)
                return;
            end
            if nargin<7 
                if ~isequal(tData, sData)
                    return;
                end
            else
                if ~areEqual
                    return;
                end
            end
            tChoices=MatBasics.LookForIds(tIdPerRow, tIdSet);
            sChoices=MatBasics.LookForIds(sIdPerRow, sIdSet);
            sizeRecalled=sum(tChoices&sChoices);
            sizeTrainingSet=sum(tChoices);
            sizeTestSet=sum(sChoices);
            sO=sizeRecalled/sizeTrainingSet;
            tO=sizeRecalled/sizeTestSet;
        end

        function [madUnits, tFreq, sFreq,tMdns, sMdns, tIsMaxSd]=GetMdnUnits(...
                tData, sData, tIdPerRow,tIdSet, sIdPerRow, sIdSet)
            madUnits=[];
            tFreq=[];
            sFreq=[];
            tMdns=[];
            sMdns=[];
            tIsMaxSd=[];
            if isempty(tData)||isempty(sData)
                return;
            end
            tChoices=MatBasics.LookForIds(tIdPerRow, tIdSet);
            sChoices=MatBasics.LookForIds(sIdPerRow, sIdSet);
            tData_=tData(tChoices, :);
            sData_=sData(sChoices, :);
            if size(sData_,1)<2
                return;
            end
            sMdns=median(sData_);
            tIsMaxSd=true(1, length(sMdns));
            if size(tData_,1)<2
                return;
            end
            tMdns=median(tData_);
            tMads=mad(tData_);
            madUnits=abs(sMdns-tMdns)./tMads;
            if any(madUnits>4)
                [tMdns;sMdns;tMads;madUnits]
                fprintf('%d parameters are over 4 dev unit away\n',...
                    sum(madUnits>4));
                disp('herp');
            end
            if nargout>1
                tFreq=sum(tChoices)/size(tData,1);
                sFreq=sum(sChoices)/size(sData,1);
            end
        end
        
        
        function outData=LogReal(fcn, inData, shiftDoNotRefit)
            if nargin<3
                shiftDoNotRefit=true;
            end
            [~,C]=size(inData);
            outData=inData;
            for i=1:C
                avoidNonReal(i)
            end
            outData=feval(fcn, outData);
            
            function avoidNonReal(col)
                %log converts 0s to inf and negatives to complex
                %thus 0s and negs must be re-fit BELOW the lowest unconverted
                %number between 0 and 1 or re-fit between .01 and 1
                %BEFORE log10 conversion
                nonRealIndex=inData(:, col)<=0;
                if ~any(nonRealIndex) 
                    return;
                end
                nonReal=outData(nonRealIndex, col);
                maxNonReal=max(nonReal);
                sum(nonRealIndex) %display count to console
                real=outData(~nonRealIndex, col);
                minNegLog=min(real(real<=1));
                if isempty(minNegLog)
                    if maxNonReal==0
                        maxRefit=1; % log10 of 1 is 0
                    else
                        maxRefit=.99;
                    end
                    minRefit=0.01;
                else
                    maxRefit=.99*minNegLog;
                    minRefit=.01*minNegLog;
                end
                if shiftDoNotRefit
                    mn=min(nonReal);
                    plus=minRefit-mn;
                    outData(:, col)=outData(:,col)+plus;
                else                    
                    refitRegion=maxRefit-minRefit;
                    minNonReal=min(nonReal);
                    nonRealRange=maxNonReal-minNonReal;
                    if  nonRealRange==0
                        nonRealRatio=.5;
                    else
                        nonRealRatio=(nonReal-minNonReal)/nonRealRange;
                    end
                    nonRealAdjusted=minRefit+(nonRealRatio*refitRegion);
                    outData(nonRealIndex, col)=nonRealAdjusted;
                end
            end
        end
        
        function [tMdns, tFreq, tMads]=GetMdns(tData, tIdPerRow, tIds)
            N=length(tIds);
            [~,C]=size(tData);
            tFreq=zeros(1,C);
            tMdns=zeros(N, C);
            tMads=zeros(N, C);
            for i=1:N
                tChoices=MatBasics.LookForIds(tIdPerRow, tIds(i));
                tData_=tData(tChoices, :);
                if size(tData_,1)<2
                    continue;
                end
                tMdns(i,:)=median(tData_);
                tMads(i,:)=mad(tData_);
                tFreq(i)=sum(tChoices)/size(tData,1);
            end
        end
        
        function [maxPerc, tPercs, sPercs]=GetMaxPercOverlap(tData, ...
                sData, tIdPerRow, tIdSet, sIdPerRow, sIdSet)
            maxPerc=[];
            tPercs=[];
            sPercs=[];
            if isempty(tData)||isempty(sData)
                return;
            end
            tChoices=MatBasics.LookForIds(tIdPerRow, tIdSet);
            tData_=tData(tChoices, :);
            sChoices=MatBasics.LookForIds(sIdPerRow, sIdSet);
            sData_=sData(sChoices, :);
            [maxPerc, tPercs, sPercs]=MatBasics.GetMaxPercOverlap2(tData_, sData_);
        end
        
        function [maxPerc, tPercs, sPercs]=GetMaxPercOverlap2(...
                tData_, sData_)
            [R,tC]=size(tData_);
            if R < 1
                return;
            end
            [R,sC]=size(sData_);
            if R < 1
                return;
            end
            if sC ~= tC
                return;
            end
            tPercs=zeros(1,sC);
            sPercs=zeros(1,sC);
            for i=1:sC
                tPercs(i)=MatBasics.PercentInRange(tData_(:, i), sData_(:, i));
                sPercs(i)=MatBasics.PercentInRange(sData_(:,i), tData_(:,i));
            end
            maxPerc=zeros(1,sC);
            l=tPercs>sPercs;
            maxPerc(l)=tPercs(l);
            maxPerc(~l)=sPercs(~l);
        end
        
        function perc=PercentInRange(first,second)
            mx=max(second);
            mi=min(second);
            perc=sum(first >= mi & first <= mx)/length(first)*100;
        end
        
        function [devUnits, tFreq, sFreq, tMdns, sMdns, tIsMaxSd]=...
                GetMad1DevUnits(tData, sData, tIdPerRow, ...
                tIdSet, sIdPerRow, sIdSet)
            devUnits=[];
            tFreq=[];
            sFreq=[];
            tMdns=[];
            sMdns=[];
            tIsMaxSd=[];
            if isempty(tData)||isempty(sData)
                return;
            end
            tChoices=MatBasics.LookForIds(tIdPerRow, tIdSet);
            sChoices=MatBasics.LookForIds(sIdPerRow, sIdSet);
            tData_=tData(tChoices, :);
            sData_=sData(sChoices, :);
            if size(sData_,1)<2
                return;
            end
            sMdns=median(sData_);
            tIsMaxSd=true(1, length(sMdns));
            sMads=mad(sData_,1);
            if size(tData_,1)<2
                return;
            end
            tMdns=median(tData_);
            tMads=mad(tData_,1);
            N=length(sMdns);
            devUnits=zeros(1, N);
            for i=1:N
                if tMads(i)>=sMads(i)
                    devUnits(i)=abs(sMdns(i)-tMdns(i))/tMads(i);
                else
                    tIsMaxSd(i)=false;
                    devUnits(i)=abs(sMdns(i)-tMdns(i))/sMads(i);
                end
            end
            if nargout>1
                tFreq=sum(tChoices)/size(tData,1);
                sFreq=sum(sChoices)/size(sData,1);
            end
        end
        
        function [devUnits, tFreq, sFreq, tMdns, sMdns, tIsMaxSd]=...
                GetMadDevUnits(tData, sData, tIdPerRow, ...
                tIdSet, sIdPerRow, sIdSet)
            devUnits=[];
            tFreq=[];
            sFreq=[];
            tMdns=[];
            sMdns=[];
            tIsMaxSd=[];
            if isempty(tData)||isempty(sData)
                return;
            end
            tChoices=MatBasics.LookForIds(tIdPerRow, tIdSet);
            sChoices=MatBasics.LookForIds(sIdPerRow, sIdSet);
            tData_=tData(tChoices, :);
            sData_=sData(sChoices, :);
            if size(sData_,1)<2
                return;
            end
            sMdns=median(sData_);
            tIsMaxSd=true(1, length(sMdns));
            sMads=mad(sData_);
            if size(tData_,1)<2
                return;
            end
            tMdns=median(tData_);
            tMads=mad(tData_);
            N=length(sMdns);
            devUnits=zeros(1, N);
            for i=1:N
                if tMads(i)>=sMads(i)
                    devUnits(i)=abs(sMdns(i)-tMdns(i))/tMads(i);
                else
                    tIsMaxSd(i)=false;
                    devUnits(i)=abs(sMdns(i)-tMdns(i))/sMads(i);
                end
            end
            if nargout>1
                tFreq=sum(tChoices)/size(tData,1);
                sFreq=sum(sChoices)/size(sData,1);
            end
        end
        
        function [mins, maxs]=GetMinsMaxs(data, trim)
            if nargin<2
                trim=.15;
            end
            trimUp=1+trim;
            trimDown=1-trim;
            mins=min(data);
            maxs=max(data);
            D=size(data, 2);
            for i=1:D
                if mins(i)<0
                    mins(i)=mins(i)*trimUp;
                else
                    mins(i)=mins(i)*trimDown;
                end
                if maxs(i)<0
                    maxs(i)=maxs(i)*trimDown;
                else
                    maxs(i)=maxs(i)*trimUp;
                end
            end
        end

        
        function [devUnits, tFreq, sFreq, tMdns, sMdns, tIsMaxSd]=...
                GetStdDevUnits(tData, sData, tIdPerRow, ...
                tIdSet, sIdPerRow, sIdSet)
            devUnits=[];
            tFreq=[];
            sFreq=[];
            tMdns=[];
            sMdns=[];
            tIsMaxSd=[];
            if isempty(tData)||isempty(sData)
                return;
            end
            tChoices=MatBasics.LookForIds(tIdPerRow, tIdSet);
            sChoices=MatBasics.LookForIds(sIdPerRow, sIdSet);
            tData_=tData(tChoices, :);
            sData_=sData(sChoices, :);
            if size(sData_,1)<2
                return;
            end
            sMdns=median(sData_);
            tIsMaxSd=true(1, length(sMdns));
            sStds=std(sData_);
            if size(tData_,1)<2
                return;
            end
            tMdns=median(tData_);
            tStds=std(tData_);
            N=length(sMdns);
            devUnits=zeros(1, N);
            if sIdSet==1799
                disp('huh');
            end
            for i=1:N
                if tStds(i)>=sStds(i)
                    devUnits(i)=abs(sMdns(i)-tMdns(i))/tStds(i);
                else
                    tIsMaxSd(i)=false;
                    devUnits(i)=abs(sMdns(i)-tMdns(i))/sStds(i);
                end
            end
            if nargout>1
                tFreq=sum(tChoices)/size(tData,1);
                sFreq=sum(sChoices)/size(sData,1);
            end
        end
        
        function rec=SortMatrix(matrix)
            if nargin<1
                matrix=[500:-100:100; 10:3:22; 1066:33:1198];
                matrix(2,3)=9982;
            end
            rec.matrix=matrix;
            rec.vector=matrix(:);
            [rec.sorted, rec.sortedIdxs]=sort(rec.vector);
            rec.N=length(rec.vector);
        end
        
        function [matrixValue, x, y, idx]=PokeSortedMatrix(rec, i)
            [R, C]=size(rec.matrix);
            idx=rec.sortedIdxs(i);
            [x, y]=ind2sub([R C], idx);
            matrixValue=rec.matrix(x,y);
        end
        
        function rank=Rank(score, scoreMin, scoreMax, maxRank)
            range=scoreMax-scoreMin;
            rank=round((score-scoreMin)/range*maxRank);
            if rank<1
                rank=1;
            elseif rank>maxRank
                rank=maxRank;
            end
        end
        
        function i=GetIdxOfMore(q, item, maxSize)
            size=length(q);
            if size<maxSize
                i=length(q)+1; % add to end
                return;
            end
            [mx, i]=max(q);
            if mx<=item
                i=-1;
            end
        end
        function lvl=HandleIdx(a, lvl)
            if lvl>length(a)
                lvl=mod(lvl, length(a));
                if lvl==0
                    lvl=length(a);
                end
            end
        end
        function PrintStackTrace(msg, alert)
            try
                msgID='MatBasics:printTrace';
                exc= MException(msgID, msg);
                throw(exc);
            catch ex
                disp(ex.getReport);
                if nargin>1 
                    if alert==1
                        msgBox(msg);
                    else
                        msgbox(msg);
                    end
                end
            end
        end
        
        function groups=FindClusterGroups(clusterIds, oldPtrs, newPtrs)
            N1=length(clusterIds);
            newClusters=unique(newPtrs);
            N2=length(newClusters);
            
            if N2>0 && newClusters(1)==0
                N2=N2-1;
            end
            groups=cell(1, N2);
            for i=1:N1
                clue=clusterIds(i);
                idxs=unique(newPtrs(find(oldPtrs==clue)));
                assert( length(idxs)==1 );
                groups{idxs(1)}=[clue groups{idxs(1)}];
            end
        end
        
        function SetMatrixAsProp(props, name, matrix)
            [R, C]=size(matrix);
            A={};
            for r=1:R
                A{end+1}=num2str(matrix(r,:));
            end
            props.setAll(name, A);
        end
        function matrix=GetMatrixAsProp(props, name)
            matrix=[];
            A=props.getAll(name);
            R=length(A);
            if iscell(A)
                for r=1:R
                    matrix(r,:)=str2num(A{r});
                end
            else
                for r=1:R
                    matrix(r,:)=str2num(A(r));
                end
            end
        end
        
        function n=Mod(i,N)
            n=mod(i, N);
            if n==0
                n=N;
            end
        end
        function idxs=Limit(N, LIMIT)
            if N>LIMIT
                i=floor(N/LIMIT);
                if i==1
                    idxs=true(N,1);
                    idxs(1:LIMIT-N)=false;
                else
                    idxs=false(N,1);
                    idxs(1:i:LIMIT*i)=true;
                end
            else
                idxs=ones(N);
            end
        end
        
        function perc=Num2Perc(num, lower, upper)
            perc=(num-lower) / (upper-lower);
        end
        
        function num=Perc2Num(perc, lower, upper)
            num=roundTo(lower+perc*(upper-lower), .1);
        end
        
        function [polyX, polyY]=Polygon(x, y)
            try
                k=convhull(x, y);
                [polyX, polyY]=poly2cw(x(k), y(k));
            catch ex
                %disp(ex);
                polyX=[];
                polyY=[];
            end
        end
        
        function sizeOrderIdxs=GetSizeOrder(cellOfMatrices, dimension, start)
            if nargin<3
                start=1;
            end
            offset=start-1;
            N=length(cellOfMatrices);
            us=zeros(1,N-offset);
            for i=start:N
                us(i-offset)=size(cellOfMatrices{i}, dimension);
            end
            [~,b]=sort(0-us);
            sizeOrderIdxs=zeros(1,N-offset);
            for i=1:N-offset
                idx=b(i);
                sizeOrderIdxs(idx)=i;
            end
        end
        
        function [out,ok]=CheckBounds(a, N)
            n=length(a);
            ok=true;
            out=[];
            for i=1:n
                if a(i)<1||a(i)>N
                    ok=false;
                else
                    out(end+1)=a(i);
                end
            end
        end
        
        function idxs=GetSizeRankings(sizes)
            N=length(sizes);
            [~,b]=sort(0-sizes);
            idxs=zeros(1,N);
            for i=1:N
                idxs(b(i))=i;
            end
        end
        
        function c=Matrix(javaCollection)
            N=javaCollection.size;
            c=zeros(1,N);
            it=javaCollection.iterator;
            for i=1:N
                c(i)=it.next;
            end
        end
        
        function [left, right, bottom, top]=getBorder(xy, ...
                minWidths, minSteps)
            if nargin==1
                minWidths=max(xy)-min(xy)/20;
                minSteps=[20, 20];
            end
            [left,right]=MatBasics.GetSides(xy,minWidths(1), ...
                minSteps(1), 2, 1);
            [bottom,top]=MatBasics.GetSides(xy,minWidths(1), ...
                minSteps(1), 1,2);
        end
               
        function [side1, side2]=GetSides(xy, minWidth, minSteps,...
                p1, p2)
            Min=min(xy(:, p1));
            Max=max(xy(:, p1));
            N=ceil((Max-Min)/minWidth);
            if N>size(xy,1)/3
                N=floor(size(xy,1)/3);
                stepWidth=ceil( (Max-Min)/N );
            elseif N<minSteps
                N=minSteps;
                stepWidth=ceil( (Max-Min)/N );
            else
                stepWidth=minWidth;
            end
            
            side1=[];
            side2=[];
            for i=1:N
                u=Min+(i*stepWidth);
                l=u-stepWidth;
                ll=xy(:, p1)>=l;
                lu=xy(:, p1)<u;
                d=ll+lu;
                [mn, iMn]=min(xy(d==2, p2));
                if ~isempty(mn)
                    zone=xy(d==2,:);
                    [mx, iMx]=max(zone(:, p2));
                    side1(end+1, p1)=zone(iMn, p1);%u+(u-l)/2;
                    side2(end+1, p1)=zone(iMx, p1);%u+(u-l)/2;
                    side1(end,p2)=mn;
                    side2(end,p2)=mx;
                end
            end
        end
        function ok=isEqual(m1,m2)
            ok=all(size(m1)==size(m2));
            if ok
                ok=isempty(find(m1~=m2,1));
            end
        end
        function ok=isEqualNoOrder(m1,m2)
            n1=MatBasics.TotalCount(m1);
            n2=MatBasics.TotalCount(m2);
            ok=n1==n2;
            if ok
                u1=unique(m1(:));
                u2=unique(m2(:));
                n1=length(u1);
                n2=length(u2);
                ok=n1==n2;
                if ok
                    for i=1:n1
                        ok=any(u1==u2(i));
                        if ~ok
                            break;
                        end
                    end
                end
            end
        end
        
        function n=TotalCount(m)
            if isempty(m)
                n=0;
            end
            s=size(m);
            n=s(1);
            D=length(s);
            for i=2:D
                n=n*s(i);
            end
        end
            
        function out=Num2Str(mat)
            R=size(mat);
            str=num2str(mat);
            out='';
            for i=1:R
                out=[out str(i,:) '#'];
            end
        end
        
        function nums=Str2Num(mat)
            strs=strsplit(mat, '#');
            N=length(strs);
            if N>0
                nums=str2num(strs{1});
                for i=2:N
                    if ~isempty(strs{i})
                        nums(end+1,:)=str2num(strs{i});
                    end
                end
            else
                nums=[];
            end
        end
        
        function str=toString(mat, delimiter)
            %UNTITLED2 Summary of this function goes here
            %   Detailed explanation goes here
            if nargin<2
                delimiter=', ';
            end
            rows=size(mat,1);
            cols=size(mat,2);
            str='';
            n=0;
            for i = 1:rows
                for j =1:cols
                    if n>0
                        str=strcat(str, delimiter);
                    end
                    n=n+1;
                    v=mat(i,j);
                    if isnumeric(v)
                        v=num2str(v);
                    elseif iscell(v)
                        v=cell2mat(v);
                    end
                    str=strcat(str, v);
                end
            end
            
        end
        
        function html=toRoundedTable(mat, rounded, highlight, names)
            %UNTITLED2 Summary of this function goes here
            %   Detailed explanation goes here
            if nargin<3
                highlight=[];
            end
            rows=size(mat,1);
            cols=size(mat,2);
            html='<table>';
            n=0;
            all=1;
            if nargin>3
                html=[html '<thead><tr>'];
                for i = 1:cols
                    if ~isempty(find(all==highlight))
                        html=[html '<th color="red">' String.ToHtml(names{i}) '</th>'];
                    else
                        html=[html '<th>' String.ToHtml(names{i}) '</th>'];
                    end
                    all=all+1;
                end
                html=[html '</tr></thead>'];
            end
            all=1;
            for i = 1:rows
                html=[html '<tr>'];
                for j =1:cols
                    n=n+1;
                    v=mat(i,j);
                    v=String.encodeRounded(v, rounded, true);
                    if ~isempty(find(all==highlight))
                        v=['<b>' v '</b>'];
                    end
                    html=[html '<td align="right">' v '</td>'];
                    all=all+1;
                end
                html=[html '</tr>'];                
            end
        end
        
        function [x, y, maximum]=MaxIn2D(data)
            [colMax,rowIdxs]=max(data);
            [maximum,colIdx]=max(colMax);
            x=rowIdxs(colIdx);
            y=colIdx;
            
        end
        
        function m=Cell2Mat(c)
            N=length(c);
            m=zeros(1, N);
            for i=1:N
                m(i)=c{i};
            end
        end
        
        function c=Mat2Cell(m)
            N=length(m);
            c=cell(1, N);
            for i=1:N
                c{i}=m(i);
            end
        end

        function [lowIdx, lowValue]=NextMin(a,min)
            N=length(a);
            lowValue=0;
            lowIdx=0;
            for i=1:N
                if a(i)>min
                    if a(i)<lowValue
                        lowValue=a(i);
                        lowIdx=i;
                    end
                end
            end
            
        end
        function counts=CountMinEdge(matrix,mins)
            edge=[];
            N=size(matrix,2);
            for i=1:N
                edge(:,end+1)=matrix(:,i)<mins(i);                
            end
            counts=sum(edge);
        end
        
        function edgeRows=FindMinOrAbove(matrix, mins)
            [M, N]=size(matrix);
            edge=zeros(M,N);
            for i=1:N
                edge(:,i)=matrix(:,i)>=mins(i);                
            end
            edgeRows=sum(edge,2)==N;
        end

        function [edgeRows, reduced]=FindMaxOrBelow(matrix, maxs, idxs)
            if nargin<3
                [M, N]=size(matrix);
                idxs=1:N;
            else
                N=length(idxs);
                M=size(matrix,1);
            end
            try
                edge=bsxfun(@le, matrix(:, idxs), maxs);
                edgeRows=sum(edge,2)==N;
                reduced=sum(edgeRows)<M;
            catch e
                e.getReport
                edgeRows=[];
                reduced=[];
            end
        end

        function edgeRows=FindMinEdge2(m, t)
            edgeRows=(m(:,1)<t(1))+(m(:,2)<t(2))~=0;
        end
        
        function edgeRows=FindMinEdge(matrix, mins)
            edge=[];
            N=size(matrix,2);
            for i=1:N
                edge(:,end+1)=matrix(:,i)<mins(i);                
            end
            edgeRows=sum(edge,2)>0;
        end
        
        function data=PileUpAtMinimum(data, mins)
            N=size(data,2);
            for i=1:N
                l=data(:,i)<mins(i);
                if find(l,1)>0
                    %data(l, i)=mins(i);
%                    disp('hmm');
                end
            end            
        end

        function data=PileUpAtMaximum(data, maxs)
            N=size(data,2);
            for i=1:N
                l=data(:,i)>maxs(i);
                if find(l,1)>0
                    data(l, i)=maxs(i);
                end
            end            
        end

        function edgeRows=FindOnScale(matrix, mins, maxs)
            N=size(matrix,2);
            edgeRows=sum([bsxfun(@ge, matrix, mins) bsxfun(@le, matrix, maxs)],2)==2*N;
        end

        function letters=toLetters(numbers, sort1st)
            if nargin==1
                sort1st=true;
            end
            N=length(numbers);
            letters='';
            if sort1st && length(numbers)>1
                numbers=sort(numbers);
            end
            for j=1:N
                l=numbers(j);
                if j>1
                    letters=strcat(letters, ', ', String.toLetter(l));
                else
                    letters=String.toLetter(l);
                end
            end
        end
        
        function [pos, xLim, yLim]=AxesInfo(ax)
            pos=getpixelposition(ax);
            axlim=axis(ax); % Get the axis limits [xlim ylim (zlim)]
            xLim=axlim(1:2);
            yLim=axlim(3:4);
        end
        
        function [figPixelsX, figPixelsY]=AxesXy2FigPixels(dataX, dataY,...
                axOrPos, xLim, yLim)
            if nargin<4
                [pos, xLim, yLim]=MatBasics.AxesInfo(axOrPos);
            else
                pos=axOrPos;
            end
            W=diff(xLim);
            H=diff(yLim);
            figPixelsX=(pos(3)*(dataX-xLim(1))/W);
            figPixelsY=(pos(4)*(dataY-yLim(1))/H);
        end
        
        
        function varargout = dsxy2figxy(varargin)
            % dsxy2figxy -- Transform point or position from data space
            % coordinates into normalized figure coordinates
            % Transforms [x y] or [x y width height] vectors from data space
            % coordinates to normalized figure coordinates in order to locate
            % annotation objects within a figure. These objects are: arrow,
            % doublearrow, textarrow, ellipse line, rectangle, textbox
            %
            % Syntax:
            %    [figx figy] = dsxy2figxy([x1 y1],[x2 y2])  % GCA is used
            %    figpos      = dsxy2figxy([x1 y1 width height])
            %    [figx figy] = dsxy2figxy(axes_handle, [x1 y1],[x2 y2])
            %    figpos      = dsxy2figxy(axes_handle, [x1 y1 width height])
            %
            % Usage: Obtain a position on a plot in data space and
            %        apply this function to locate an annotation there, e.g.,
            %   [axx axy] = ginput(2); (input is in data space)
            %   [figx figy] = dsxy2figxy(gca, axx, axy);  (now in figure space)
            %   har = annotation('textarrow',figx,figy);
            %   set(har,'String',['(' num2str(axx(2)) ',' num2str(axy(2)) ')'])
            %
            %   Copyright 2006-2009 The MathWorks, Inc.
            
            % Obtain arguments (limited argument checking is done)
            % Determine if axes handle is specified
            if length(varargin{1}) == 1 && ishandle(varargin{1}) ...
                    && strcmp(get(varargin{1},'type'),'axes')
                hAx = varargin{1};
                varargin = varargin(2:end); % Remove arg 1 (axes handle)
            else
                hAx = gca;
            end;
            
            % Remaining args are either two point locations or a position vector
            if length(varargin) == 1        % Assume a 4-element position vector
                pos = varargin{1};
            else
                [x,y] = deal(varargin{:});  % Assume two pairs (start, end points)
            end
            
            % Get limits
            axun = get(hAx,'Units');
            set(hAx,'Units','normalized');  % Make axes units normalized
            axpos = get(hAx,'Position');    % Get axes position
            axlim = axis(hAx);              % Get the axis limits [xlim ylim (zlim)]
            axwidth = diff(axlim(1:2));
            axheight = diff(axlim(3:4));
            
            % Transform from data space coordinates to normalized figure coordinates
            if exist('x','var')     % Transform a and return pair of points
                varargout{1} = (x - axlim(1)) * axpos(3) / axwidth + axpos(1);
                varargout{2} = (y - axlim(3)) * axpos(4) / axheight + axpos(2);
            else                    % Transform and return a position rectangle
                pos(1) = (pos(1) - axlim(1)) / axwidth * axpos(3) + axpos(1);
                pos(2) = (pos(2) - axlim(3)) / axheight * axpos(4) + axpos(2);
                pos(3) = pos(3) * axpos(3) / axwidth;
                pos(4) = pos(4) * axpos(4 )/ axheight;
                varargout{1} = pos;
            end
            
            % Restore axes units
            set(hAx,'Units',axun)
        end
        
        function dataXy=figXyToData(obj, ax, whToo)
            u=get(obj, 'units');
            set(obj, 'units', 'normalized');
            objP=get(obj, 'position');
            set(obj, 'units', u);
            u=get(ax, 'units');
            set(ax, 'units', 'normalized');
            axP=get(ax, 'position');
            X=adjust(1, xlim(ax));
            Y=adjust(2, ylim(ax));
            
            if nargin>2 && whToo
                W=adjustWh(3, xlim(ax));
                H=adjustWh(4, ylim(ax));
                dataXy=[X Y W H];
            else
                dataXy=[X Y];
            end
            set(ax, 'units', u);
            function out=adjust(idx, lim)
                p=abs( objP(idx)-axP(idx) );
                sz=p/axP(idx+2);
                fr=sz*(lim(2)-lim(1));
                if objP(idx)>axP(idx)
                    out=lim(1)+fr;
                else
                    out=lim(1)-fr;
                end
            end 
            function out=adjustWh(idx, lim)
                out=objP(idx)/axP(idx)*(lim(2)-lim(1));
            end 
 
        end
        
        function figXy=dataXyToFig(ax, xy)
            u=get(ax, 'units');
            set(ax, 'units', 'normalized');
            axP=get(ax, 'position');
            X=adjust(1, xlim(ax));
            Y=adjust(2, ylim(ax));
            figXy=[X Y];
            set(ax, 'units', u);
            function out=adjust(idx, lim)
                p=abs(xy(idx));
                r=(lim(2)-lim(1));
                fr=(p-lim(1))/r;
                if xy(idx)>=0
                    out=axP(idx)+fr*axP(idx+2);
                else
                    out=axP(idx)-fr*axP(idx+2);
                end
            end 
        end
        
        function figXy=dataWhToFig(ax, wh)
            u=get(ax, 'units');
            set(ax, 'units', 'normalized');
            axP=get(ax, 'position');
            X=adjust(1, xlim(ax));
            Y=adjust(2, ylim(ax));
            figXy=[X Y];
            set(ax, 'units', u);
            function out=adjust(idx, lim)
                out=wh(idx)/(lim(2)-lim(1))*axP(idx+2);
            end 
        end
        
        function idxBoolVector=NthTrueInBoolVector(boolVector, nTh)
            N=length(boolVector);
            count=0;
            for i=1:N
                if boolVector(i)
                    count=count+1;
                    if count==nTh
                        idxBoolVector=i;
                        return;
                    end
                end
            end
            idxBoolVector=0;
        end
        
        function str=Encode(m)
            [R, C]=size(m);
            if R==1 && C==1
                str=num2str(m(1,1));
                return;
            end
            str='[';
            for r=1:R
                if r>1
                    str=[str ';'];
                end
                for c=1:C
                    if c>1
                        str=[str ','];
                    end
                    str=[str num2str( m(r,c), 15)];
                end
            end
            str=[str ']'];
        end
        
        function  idxs=SmoothLow(ln)
            idxs=[1];
            N=length(ln);
            pI=1;
            pV=ln(pI);
            nI=2;
            while nI<N
                nV=ln(nI);
                if nV<pV
                    pI=nI;
                else
                    if nV>pV
                        [~, mI]=min(ln(nI+1:N));
                        pI=nI+mI;
                    else
                        nI=nI+1;
                        continue;
                    end
                end
                pV=ln(pI);
                idxs(end+1)=pI;
                nI=pI+1;
            end
            idxs(end+1)=N;
        end
        
        function  idxs=SmoothLowExact(ln)
            idxs=[1];
            N=length(ln);
            pI=1;
            pV=ln(pI);
            nI=2;
            la=1;
            while nI<N
                nV=ln(nI);
                if nV<pV
                    pI=nI;
                else
                    good=false;
                    if nV>pV
                        for la=nI+1:N
                            if ln(la)<nV
                                pI=la;
                                good=true;
                                break;
                            end
                            % avoid "gulleys" into which other clusters can dip
                            if la>nI+6
                                pI=la;
                                good=true;
                                break;
                            end
                        end
                    end
                    if ~good
                        nI=nI+1;
                        continue;
                    end
                end
                pV=ln(pI);
                idxs(end+1)=pI;
                nI=pI+1;
            end
            idxs(end+1)=N;
        end
        
        function  idxs=SmoothHigh(ln)
            idxs=[1];
            N=length(ln);
            pI=1;
            pV=ln(pI);
            nI=2;
            while nI<N
                nV=ln(nI);
                if nV>pV
                    pI=nI;
                else
                    if nV<pV
                        [~, mI]=max(ln(nI+1:N));
                        pI=nI+mI;
                    else
                        nI=nI+1;
                        continue;
                    end
                end
                pV=ln(pI);
                idxs(end+1)=pI;
                nI=pI+1;
            end
            idxs(end+1)=N;
        end
    
        function  idxs=SmoothHighExact(ln)
            idxs=[1];
            N=length(ln);
            pI=1;
            pV=ln(pI);
            nI=2;
            la=1;
            while nI<N
                nV=ln(nI);
                if nV>pV
                    pI=nI;
                else
                    good=false;
                    if nV<pV
                        for la=nI+1:N
                            if ln(la)>nV
                                pI=la;
                                good=true;
                                break;
                            end
                            % avoid "gulleys" into which other clusters can dip
                            if la>nI+6
                                pI=la;
                                good=true;
                                break;
                            end
                        end
                    end
                    if ~good
                        %pI=nI;
                        %pV=ln(pI);
                        %idxs(end+1)=pI;
                        nI=nI+1;
                        continue;
                    end
                end
                pV=ln(pI);
                idxs(end+1)=pI;
                nI=pI+1;
                
            end
            idxs(end+1)=N;
        end
        
        function missing=getMissing(left, right)
            missing=java.util.ArrayList;            
            N1=length(left);
            N2=length(right);
            for i=1:N1
                found=false;
                for j=1:N2
                    if left(i)==right(j)
                        found=true;
                    end
                end
                if ~found
                    missing.add(left(i));
                end
            end
        end
        function [added,subtracted]=diff(old, new)
            subtracted=MatBasics.getMissing(old,new);
            added=MatBasics.getMissing(new,old);
        end
        
        function out=ToStr(in)
            out='';
            if iscell(in)
                for i=1:length(in)
                    out=[out ' ' MatBasics.ToStr(in{i}) ' '];
                end
            elseif isnumeric(in)
                out=num2str(in);
            elseif ischar(in)
                out=in;
            end            
        end
        
        function lines=ToStrs(mat)
            
            N=size(mat);
            R=N(1);
            C=N(2);
            lines=cell(1,R);
            
            for i=1:R
                ln=[];
                
                for j=1:C
                    if j<C
                        ln=[ln num2str(mat(i,j)) ','];
                    else
                        ln=[ln num2str(mat(i,j))];
                    end
                end
                lines{i}=sprintf('%s\n', ln);
            end
        end
        function idx=closest(ss,pr)
            [~, idx]=min(abs(ss-pr));
        end
        
        function point=MidRight(xy)
            ym=mean(xy(:,2));
            ys=std(xy(:,2));
            if ys==0
                xm=max(xy(:,1));
            else
                l=xy(:,2)<ym+ys & xy(:,2)>ym-ys;
                xm=max(xy(l,1));
            end
            point=[xm,ym];
        end
        
        function idxs=UpToPercent(cellA, percent, startIdx)
            allSizes=cellfun(@(n)size(n,1), cellA);
            sizes=allSizes(startIdx+1:end);
            [ss, ii]=sort(sizes, 'descend');
            ssp=cumsum(ss)/sum(ss);
            idx=MatBasics.closest(ssp, percent);
            idxs=zeros(1,idx);
            for i=1:idx
                idxs(i)=ii(i);
            end
        end
        
        function [point, i]=Closest(xy, X, Y)
            try
                xy2=[abs(xy(:,2)-Y), abs(xy(:,1)-X)];
                xy3=xy2(:,1)+xy2(:,2);
                [~,i]=min(xy3);
                point=xy(i,:);
            catch ex
                point=[];
                i=0;
                disp(ex);
            end
        end
        
        function set=AllArrays(array, set, legal)
            if nargin<3
                hasLegal=false;
                if nargin<2
                    set=java.util.LinkedHashSet;
                end
            else
                hasLegal=true;
            end
            set.add(num2str(sort(array)));
            N=length(array)-1;
            for i=1:N
                v=nchoosek(array, i);
                N2=length(v);
                if i>1 && hasLegal
                    for j=1:N2
                        v2=v(j,:);
                        if all(ismember(v2, legal{v2(1)}));
                            set.add(java.lang.String(num2str(sort(v2))));
                        end
                    end
                    
                else
                    for j=1:N2
                        set.add(java.lang.String(num2str(sort(v(j,:)))));
                    end
                end
            end
        end
        function set=AllArrays2(array, set, legal)
            set.add(num2str(sort(array)));
            N=length(array)-1;
            for i=1:N
                v=nchoosek(array, i);
                N2=length(v);
                if i>1
                    for j=1:N2
                        v2=v(j,:);
                        N3=length(v2);
                        bad=false;
                        for k=1:N3
                            num3=v2(k);
                            legals=legal{num3};
                            if ~any(ismember(v2, legals(2:end)));
                                bad=true;
                                break;
                            end
                        end
                        if ~bad
                            set.add(java.lang.String(num2str(sort(v2))));
                        end
                    end
                else
                    for j=1:N2
                        set.add(java.lang.String(num2str(sort(v(j,:)))));
                    end
                end
            end
        end
        function [neighbors,neighborsAndNans, selfIdx]=NeighborHoodInd(ind,grid,radius)
            if nargin<3
                radius=1;
            end
            [x,y]=ind2sub(size(grid), ind);
            if nargout>1
                [neighbors,neighborsAndNans,selfIdx]=...
                    MatBasics.NeighborHood(x,y,grid, radius);
            else
                neighbors=MatBasics.NeighborHood(x,y,grid, radius);
            end
        end
        
        function [xIdxs,yIdxs, bins]=NeighborIdxs(x,y,lX, lY)
            xIdxs=[];
            yIdxs=[];
            for i=-1:1
                nx=x+i;
                if nx>0 && nx<=lX
                    for j=-1:1
                        if i==0&&j==0
                            continue;
                        end
                        ny=y+j;
                        if ny>0 && ny<=lY
                            xIdxs(end+1)=nx;
                            yIdxs(end+1)=ny;
                        end
                    end
                end
            end
            bins=sub2ind([lX, lY], xIdxs, yIdxs);
        end
        
        function [neighbors,neighborsAndNans, selfIdx]=...
                NeighborHood(x,y,grid, radius)
            if nargin<4 || isempty(radius)
                radius=1;
            end
            [lX, lY]=size(grid);
            neighbors=nan(3,3);
            for i=-radius:radius
                nx=x+i;
                if nx>0 && nx<=lX
                    for j=-radius:radius
                        ny=y+j;
                        if ny>0 && ny<=lY
                            neighbors(radius+1+i,radius+1+j)=grid(nx,ny);
                        end
                    end
                end
            end
            if nargout>1
                selfIdx=5;
                if any(isnan(neighbors(:)))
                    neighborsAndNans=neighbors;
                    neighbors=[];
                    for i=1:3
                        rr=neighborsAndNans(i,:);
                        if any(~isnan(rr))
                            neighbors(end+1,:)=rr(~isnan(rr));
                        end
                        if i==2
                            if isnan(rr(3))
                                selfIdx=length(neighbors(:));
                            else
                                selfIdx=length(neighbors(:))-1;
                            end
                        end
                    end
                else
                    neighborsAndNans=[];
                end
            end
        end
        
        function out=debugN(x,y,grid)
            s={'southEast', 'east', 'northEast', 'south', 'self', 'north', 'southWest', 'west', 'northWest'};
            ind=1;
            out={};
            for i=1:-1:-1
                for j=1:-1:-1
                    disp(s{ind});
                    ind=ind+1;
                    r=MatBasics.Neighborhood(x+j,y+i,grid);
                    disp(r);
                    out{end+1}=r;
                end
            end
        end
        function s=sprintf(mat, fmt, prefix)
            if nargin<3
                prefix='';
            else
                prefix=sprintf('%s\t', prefix);
            end
            s='';
            [rows,cols]=size(mat);
            for r=1:rows
                s=[s prefix];
                for c=1:cols
                    if c==cols
                        s=[s sprintf(fmt, mat(r,c))];
                    else
                        s=[s sprintf([fmt '\t'], mat(r,c))];
                    end
                end
                s=sprintf('%s\n', s);
            end
        end
        
        function html=ToHtml(matrix, columnLabels, rowLabels, ...
                boldThreshold, bestPairsMax, rounded)
            [rows,cols]=size(matrix);
            if nargin<6
                rounded=6;
                if nargin<5
                    bestPairsMax=[];
                    if nargin<4
                        boldThreshold=[];
                        if nargin<3
                            rowLabels=[];
                            if nargin<2
                                columnLabels=[];
                            end
                        end
                    end
                end
            end
            if isempty(columnLabels)
                columnLabels=cell(1,cols);
                for col=1:cols
                    columnLabels{col}=num2str(col);
                end
            end
            if isempty(rowLabels)
                rowLabels=columnLabels;
            end
            hasBestPairs=~isempty(bestPairsMax);
            if hasBestPairs
                [pairs, unpaired]=MatBasics.FindBestPairs(matrix, bestPairsMax);
                
            end
            html='<table cellpadding="3" border="1"><tr><th></th>';
            for col=1:cols
                html=[html '<th><font color="blue">' columnLabels{col} '</font></th>'];
            end
            html=[html '</tr>'];
            
            for row=1:rows
                html=[html '<tr><td><font color="blue">' rowLabels{row} '</font></td>'];
                for col=1:cols
                    v=matrix(row,col);
                    num=String.encodeRounded(v,rounded);
                    if ~isempty(boldThreshold) && v>=boldThreshold
                        html=[html '<td align="right"><b>' num '</b></td>'];
                    else
                        if hasBestPairs 
                            if  (ismember([row col], pairs, 'rows') || ...
                                    ismember([col row], pairs, 'rows'))
                                html=[html '<td align="right"><b>' num '</b></td>'];
                            else
                                if ~isempty(find(unpaired==row,1))
                                    row_=matrix(row,:);
                                    row_(row)=bestPairsMax;
                                    [~, I]=min(row_);
                                    if I==col
                                        html=[html '<td align="right"><font color="red">' num '</font></td>'];
                                    else
                                        html=[html '<td align="right">' num '</td>'];
                                    end
                                else
                                    html=[html '<td align="right">' num '</td>'];
                                end
                            end
                        else
                            html=[html '<td align="right">' num '</td>'];
                        end
                    end
                end
                html=[html '</tr>'];
            end
            html=[html '</table>'];
        end
        
        function html=ToHtmlIntegers(matrix, columnLabels, rowLabels, ...
                notes1, notes2)
            [rows,cols]=size(matrix);
            if nargin<4
                noteCols=0;
                nNotes=0;
                if nargin<3
                    rowLabels=[];
                    if nargin<2
                        columnLabels=[];
                    end
                end
            else
                nNotes=length(notes1);
                noteCols=2;
                assert(nNotes==length(notes2));
            end
            if isempty(columnLabels)
                columnLabels=cell(1,cols);
                for col=1:cols
                    columnLabels{col}=num2str(col);
                end
            end
            if isempty(rowLabels)
                rowLabels=columnLabels;
            end
            html='<table cellpadding="3" border="1"><tr><th></th>';
            app=BasicMap.Global;
            for col=1:cols+noteCols
                html=[html '<th width=40px><font color="blue">' ...
                    app.smallStart columnLabels{col} app.smallEnd ...
                    '</font></th>'];
            end
            html=[html '</tr>'];
            mx=max(matrix);
            for row=1:rows
                html=[html '<tr><td><font color="blue">' rowLabels{row} '</font></td>'];
                for col=1:cols
                    v=matrix(row,col);
                    num=String.encodeInteger(v);
                    
                    if v==mx(col)
                        html=[html '<td align="right"><b>' num '</b></td>'];
                    else
                        html=[html '<td align="right">' num '</td>'];
                    end                    
                end
                if nargin>3
                    if row>nNotes
                        html=[html '<td colspan="2"></td>'];
                    else
                        html=[html '<td>' notes1{row} '</td><td>' notes2{row} '</td>'];
                    end
                end
                html=[html '</tr>'];
            end
            html=[html '</table>'];
        end
        
        function [matrix, columnLabels]=ReadMatrix(file)
            matrix=[];
            columnLabels={};
            fid = fopen(file);
            if (fid<0)
                return;
            end
            N=0;
            try
                mustBeOneColumn=false;
                line=[];
                lines={};
                while N<2 || String.StartsWith(line, '<') || isempty(line)
                    line = fgetl(fid);
                    if ischar(line)
                        line=strtrim(line);
                        columnLabels = textscan(line, '%s', 'Delimiter', '\t');
                        N = length(columnLabels{1});
                        mustBeOneColumn=N==1;
                        lines{end+1}=line;
                    else
                        break;
                    end
                end
                matrix=zeros(N,N);
                if mustBeOneColumn
                    numbers = textscan(lines{end}, '%f', 'Delimiter', '\t');
                    columnLabels = textscan(lines{end-1}, '%s', 'Delimiter', '\t');
                    matrix(1,:) = numbers{1}';
                else
                    for i = 1:N
                        line = fgetl(fid);
                        numbers = textscan(line, '%f', 'Delimiter', '\t');
                        matrix(i,:) = numbers{1}';
                    end
                end
                fclose(fid);
                columnLabels=columnLabels{1};
            catch ex
                fclose(fid);
                matrix=[];
                columnLabels={};
                rethrow(ex);
            end
        end
        
        
        function [slopeIdxVec, peaks, clusters, slopeIdxMat ]=...
                GetClusters(matValues, isBackground, peakSigTest)
            %[slopeIdxVec, peaks, clusters, slopeIdxMat ]=...
            %    GetClusters(matValues, isBackground, peakSigTest)
            %matValues -> 2D matrix of any v alues
            %isBackground (optional) -> vector of true/false to exclude 
            %cells in matrix
            %peakSigTest (optional) -> vector of values that peaks in 
            %matValues must be greater than to be considered a cluster
            
            sz=size(matValues);
            assert(length(sz)==2);
            assert(sz(1)==sz(2));
            M_=sz(1);
            MM=M_^2;
            %The 3d matrix neighborValues finds most dense of 8 neighbors
            %in the 2D neighborhood by having MxMx9 cubes that are
            %rotated so that the statement max(F,[],3) sees ONLY
            %each lattice point's 8 neighbors plus self which is 0
            %the 1-9 indices in matrix space is
            %   1=southEast, 2=east, 3=northEast
            %   4=south 5=none/self 6=north
            %   7=southWest, 8=west, 9=northWest 
            %In graph plot() space is
            %   1=northEast, 2=north, 3=northWest
            %   4=east 5=none/self 6=west
            %   7=southEast, 8=south , 9=southWest
            
            neighborValueCube=zeros([M_ M_ 9]);
            matrixIndiceCube=reshape(1:MM,[M_ M_]);
            neighborIdxs=zeros([M_ M_ 9]);
            neighborIdx=1;
            for col=-1:1
                for row=-1:1
                    neighborValueCube(:,:,neighborIdx)=circshift(matValues,[row col]);  %matching up densities of neighbors
                    neighborIdxs(:,:,neighborIdx)=circshift(matrixIndiceCube,[row col]);  %matching up corresponding indices of neighbors
                    neighborIdx=neighborIdx+1;
                end
            end
            
            neighborValueCube(end,:,[1 4 7])=nan; %don't let density at one edge switch to opposite edge
            neighborValueCube(:,end,[1 2 3])=nan;
            neighborValueCube(1,:,[3 6 9])=nan;
            neighborValueCube(:,1,[7 8 9])=nan;
            
            %% (07/2011 RJ) Normalize by length of distance vector to actually select steepest gradient (not in paper)
            nvcCenter = cat(3, neighborValueCube(:,:,5), neighborValueCube(:,:,5), neighborValueCube(:,:,5), neighborValueCube(:,:,5));
            neighborValueCube(:,:,[1 3 7 9]) = nvcCenter + (neighborValueCube(:,:,[1 3 7 9]) - nvcCenter)/ sqrt(2);
            [~, maxNeighIdxs]=max(neighborValueCube,[],3); %find 1-9 neighbor index of where the max density is
            slopeIdxMat=zeros(M_,M_);
            for i=1:9
                if i~=5
                    a=maxNeighIdxs==i;  %finds max neighbors in ith direction that have 
                    thisP=neighborIdxs(:,:,i); %sufficiently large slope grid points indices in ith direction
                    slopeIdxMat(a)=thisP(a);  %makes association pointers
                    %disp( sprintf('neighbor# %d, %d, %d!', i, sum(a(:)), sum(thisP(:))));
                end
            end
            slopeIdxMat(matValues==0)=0;
            slopeIdxVec=reshape(slopeIdxMat, [1 MM]);
            if nargin>1 && ~isempty(isBackground)
                assert(length(isBackground)==MM);
                slopeIdxVec(isBackground)=-1;
            end
            unassigned=find(slopeIdxVec==0); %the indices of gridpoints with no association pointers
            pointed_to=slopeIdxVec(slopeIdxVec>0); %the indices of gridpoints with association pointers into them
            peaks=intersect(unassigned,pointed_to);  %the indices of pathends: they have pointers into them but no pointers out of them            
            if nargin>2
                clusters=slopeIdxVec;
                if isempty(peakSigTest)
                    sigPeaks=true(1, length(peaks));
                else
                    assert(length(peakSigTest)==MM);
                    f=reshape(matValues, [1 MM]);
                    sigPeaks=f(peaks) >= peakSigTest(peaks);  %these pathends' densities are significant
                end
                
                clusters(peaks(sigPeaks))= -1 - peaks(sigPeaks);  %assign these pathends pointers to dummy state representing a cluster
                insigPeaks=peaks(~sigPeaks);  %pathends that have insignificant densities
                isInsig=false(1,MM);
                isInsig(insigPeaks)=true;  %keeps track of gridpoints that are on paths to insig_pathends
                
                while ~isempty(insigPeaks)
                    to_insigs=ismember(clusters,insigPeaks);  %finds gridpoints that point into insig_pathends
                    insigPeaks=find(to_insigs);  %updates insig_pathends to be the gridpoints that point into them
                    isInsig(to_insigs)=true;  %adds gridpoints that are on paths to insig_pathends
                end
                
                clusters(isInsig)=-1;  %sends to background all gridpoints on path to an insignificant pathend
                
            end
        end
        
        function RunLater(cb, secs)
            tmr=timer;
            tmr.StartDelay=secs;
            tmr.TimerFcn=@callback;
            start(tmr);
            function callback(h,e)
                try
                    feval(cb, h, e);
                catch ex
                    ex.getReport
                end
            end
        end
        
        function [avgDif, outDif]=CompareDistances2D(dist1, dist2, ...
                out1, out2)
            dif=abs(dist1-dist2);
            avgDif=median(dif(:));
            thisDifN=sum(out1);
            outDif=out1==out2;
            outDif=sum(outDif);
            outDif=outDif/thisDifN;
        end
        
        function [dist, avg, mx, out, outNum, unq2Vec]=...
                ComputeDistance2D(data, mins, maxs,  unq2Vec, M)
            if nargin<5
                M=64;
            end
            B=Binning.ToGridXy(data, M, mins, maxs);
            if nargin>3
                B=B(unq2Vec,:);
            else
                [B, unq2Vec, vec2Unq]=unique(B, 'rows');
            end
            dist=pdist2(B,B);
            avg=median(dist(:));
            mx=max(dist(:));
            dev=mad(dist(:));
            ee=pdist2(B,B, 'Euclidean', 'Smallest', 2);
            mn=ee(2,:);
            out=mn>avg+1.3;
            outNum=sum(out);
        end
        
        function inds=IdentiyMatrixHalfInds(matrixWidth, inds)
            if nargin<2
                subs=nchoosek([1:matrixWidth], 2);
            else
                subs=nchoosek(inds, 2);
            end
            inds=sub2ind([matrixWidth matrixWidth], subs(:,1), subs(:,2));
        end
        
        function ok=DEBUG(context)
            ok=false;
            if strcmpi(context, 'meanunits')
                ok=true;
            end
        end
        
        function [dataInTemplatePosition, reducedParams]=...
                ReOrg(data, columns, ...
                templateColumns, allow1EmptyLikeFMO, mustBeSameColumnCount)
            if nargin<5
                mustBeSameColumnCount=true;
                if nargin<4
                    allow1EmptyLikeFMO=false;
                end
            end
            reducedParams=[];
            [R,C]=size(data);
            N=length(columns);
            if N~=C || N~=length(templateColumns)
                if mustBeSameColumnCount || N<length(templateColumns)
                    dataInTemplatePosition=[];
                    return;
                end
                [allFound, templateIdxs]=StringArray.Find(...
                    templateColumns, columns, allow1EmptyLikeFMO);
                if ~allFound
                    dataInTemplatePosition=[];
                else
                    dataInTemplatePosition=data(:, templateIdxs);
                    reducedParams=templateIdxs;
                end
            else
                [allFound, templateIdxs]=StringArray.Find(columns, ...
                    templateColumns, allow1EmptyLikeFMO);
                if ~allFound
                    dataInTemplatePosition=[];
                else
                    dataInTemplatePosition=zeros(R,C);
                    dataInTemplatePosition(:,templateIdxs)=data;
                    [~, reducedParams]=StringArray.Find(...
                        templateColumns, columns, allow1EmptyLikeFMO);
                end
            end
            
        end
        
        function ScaleBins(ax, nBins, data)
            lbls{1}=get(ax, 'xticklabel');
            lbls{2}=get(ax, 'yticklabel');
            lbls{3}=get(ax, 'zticklabel');
            for i=1:3
                mn=min(data(:,i));
                mx=max(data(:,i));
                r=mx-mn;
                pt=r/nBins;
                l=lbls{i};
                old=l;
                N=length(l);
                for j=1:N
                    v=l{j};
                    n=str2double(v)*pt;
                    np=String.encodeRounded(mn+n, 1);
                    l{j}=np;
                end
                lbls{i}=l;
            end
            set(ax, 'xticklabel', lbls{1});
            set(ax, 'yticklabel', lbls{2});
            set(ax, 'zticklabel', lbls{3});
            
            
        end
        
        function idxs=OrderByCloseness(xy, start)
            if nargin<2
                start=[0 0];
            end
            idxs=[];
            while ~all(isnan(xy(:,1))) && ~all(isnan(xy(:,2)))
                [~,idx]=pdist2(xy, start, 'euclidean', 'smallest', 1);
                if isempty(idx) || idx==0
                    idxs=1:size(xy,1);
                    return;
                end
                idxs(end+1)=idx;
                start=xy(idx,:);
                xy(idx,:)=nan;
            end
        end
        
        function u=Unique(vector, range, preserveOrder)
            if nargin<3
                preserveOrder=true;
                if nargin<2
                    range=[];
                end
            end
            if ~isempty(range)
                v=vector(vector>=range(1) & vector<=range(2));
            else
                v=vector;
            end
            [u, I]=unique(v);
            if preserveOrder
               u=u(I);
            end
        end
        
        function counts = HistCounts(x, labels)
            if isempty(labels)
                counts = [];
                return;
            end
            if ~any(isinf(labels))
                labels(end+1) = inf;
            end
            counts = histcounts(x, labels);
        end
        
        function cnts=HistCountsTally(labels, ids)
            u=unique(labels);
            c1=MatBasics.HistCounts(labels, u);
            [c2,I]=sort(c1, 'descend');
            N=length(u);
            idxs=zeros(1, N);
            for i=1:N
                id=u(i);
                found=find(ids==id, 1);
                if ~isempty(found)
                    idxs(i)=found;
                end
            end
            nIds=length(ids);
            cnts=zeros(1,nIds);%plus 1 for not found
            
            for i=1:N
                idx=I(i);
                idx=idxs(idx);
                if idx==0
                    idx=nIds+1;
                end
                cnts(idx)=c2(i);
            end
        end
        
        function str=HistCountsText(labels, ids, names, doId)
            if nargin<4
                doId=false;
            end
            u=unique(labels);
            c1=MatBasics.HistCounts(labels, u);
            [c2,I]=sort(c1, 'descend');
            N=length(u);
            idxs=zeros(1, N);
            for i=1:N
                id=u(i);
                found=find(ids==id, 1);
                if ~isempty(found)
                    idxs(i)=found;
                end
            end
            sb=java.lang.StringBuilder;
            nIds=length(ids);
            if doId
                for i=1:N
                    if i>1
                        sb.append(',');
                    end
                    idx=I(i);
                    idx=idxs(idx);
                    if idx==0
                        name='no subset';
                        id=nan;
                        idx=nIds+1;
                    else
                        name=names{idx};
                        id=ids(idx);
                    end
                    sb.append([num2str(id) '=' name '(' ...
                        String.encodeInteger(c2(i)) ')' ]);
                    
                end
            else                
                for i=1:N
                    if i>1
                        sb.append(',');
                    end
                    idx=I(i);
                    idx=idxs(idx);
                    if idx==0
                        name='no subset';
                        idx=nIds+1;
                    else
                        name=names{idx};
                    end
                    sb.append([name '(' String.encodeInteger(c2(i)) ')' ]);
                end
            end
            str=char(sb.toString);
        end
        
        function [labels, ids, names]=TestHistDescribeCounts
            X=36;
            N=1000;
            ID=1203;
            third=floor(X/3);
            labels=randi(third, 1, N);
            for i=1:2:N
                r=mod(i,2);
                labels(i)=labels(i)+labels(i)*r;
            end
            ids=1:X;
            ids=ID+ids;
            labels=ID+labels+3;
            labels(1:5:end)=0;
            names=cell(1,X);
            for i=1:X
                names{i}=['id=' num2str(ID+i)];
            end
        end
        
    end
end
