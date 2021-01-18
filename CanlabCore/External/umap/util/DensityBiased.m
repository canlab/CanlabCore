%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef DensityBiased<handle
    properties(SetAccess=private)
        gridPoints;
        gridPtsOnDataScale;        
        M;
        mxDistance;
        xyData;
        kdTree;
        N;
        nnVec;
        eigVec;
        dnsVec;
        nearN;
        onScale;
        gridIdxs;
    end
    
    methods(Access=private)
        function this=DensityBiased(xyGridIdxs, xyData, xyMin, xyMax, M)
            this.M=M;
            [this.xyData, xyMin, xyMax]=...
                DensityBiased.NormalizeIfNotLogicle(xyData, xyMin, xyMax);
            this.kdTree=KDTreeSearcher(this.xyData);            
            [this.gridPtsOnDataScale, ptsGridScale]=...
                DensityBiased.ScaleGridPoints(M, xyMin, xyMax);
            [~, this.mxDistance]=...
                DensityBiased.BandWidth(M, this.gridPtsOnDataScale);
            [this.gridPoints, this.N]=DensityBiased.RemoveEmptyPoints( ...
                xyGridIdxs, [ptsGridScale(:,2) ptsGridScale(:,1)], M, 8);
            this.nearN=DensityBiased.DFLT_NEIGH(size(xyData,1));
            [this.gridIdxs, this.onScale]=Binning.MakeUnivariate(...
                this.xyData, xyMin, xyMax, M);
        end
        
        function calcFast(this)
            gps=this.gridPtsOnDataScale(this.gridPoints,:);
            [idxs, dists]=knnsearch(this.kdTree, gps, 'K', this.nearN);
            [R, C]=size(idxs);
            closeGps=reshape(gps, R, 1, 2);
            closeD=reshape(this.xyData(idxs, :), R, C, 2);
            finalD=closeD-closeGps;
            d1=finalD(:,:,1);
            d2=finalD(:,:,2);
            d1(dists>this.mxDistance)=nan;
            d2(dists>this.mxDistance)=nan;
            r11=nansum(d2'.^2);
            r22=nansum(d1'.^2);
            r12=nansum(d1'.*d2');
            cnts=sum(dists'<this.mxDistance);
            result=[r11./cnts;r12./cnts; r12./cnts; r22./cnts];
            dets=result(4,:).*result(1,:)-result(2,:).*result(3,:);
            densities=real(1./sqrt(dets));
            if any(isinf(densities))
                densities(isinf(densities))=0;
            end
            densities(cnts<2)=0;
            this.dnsVec=zeros(1, this.M^2);
            this.dnsVec(this.gridPoints)=densities;
            this.nnVec=zeros(1, this.M^2);
            this.nnVec(this.gridPoints)=cnts;
            eigT=(result(4,:) + result(1,:));
            eigSQ=sqrt(eigT.^2/4-dets);
            eigNum=eigT/2-eigSQ;
            eigDen=eigT/2+eigSQ;
            eigRat=eigNum./eigDen;
            this.eigVec=zeros(1, this.M^2);
            this.eigVec(this.gridPoints)=eigRat;
        end
    end
    methods(Static)
        function verbose=DEBUG
            verbose=0;
        end
        
        % xyNonGridedData is compensated and logicle values (unless
        %           scatter on linear scale)  
        % nearestNeighbors is usually 256 or higher
        % M is size of grid ... usually 256
        % xyMin is lowest compensated & transformed value
        %   possible for both X and Y parameter
        %% xyMax is highest compensated & transformed value
        %   possible for both X and Y parameter.  Always 1 (unless
        %           scatter on linear scale) 
        
        
        % In different contexts XY coordinates increase differently
        % context MatLab GUI:   X left to right; Y  bottom to top
        % context MatLab matrix: X top to bottom; Y left to right
        % context JAVA GUI:  X left to right & Y top to bottom
        % context MatLab vector:  X left to right & Y top to bottom
        function [dataScale, gridScale]=ScaleGridPoints(M, mins, maxs)
            MM=M^2;
            [cols, rows]=ind2sub([M M], 1:MM);
            gridScale=[rows;cols]';
            [Y, X]=...
                Binning.ToDataScale(gridScale, mins, maxs, M);
            dataScale=[X Y];
        end

        function [idxs, N]=RemoveEmptyPoints(eventIdxs, gridXy, M, gridDistance)
            eventL=false(1,M^2);
            eventL(unique(eventIdxs))=true;
            nonEventGridXy=gridXy(~eventL, :);
            eventGridXy=gridXy(eventL, :);
            [D, I]=pdist2(eventGridXy, nonEventGridXy, ...
                'euclidean', 'Smallest', 1);
            close_=D<gridDistance;
            closeNonEventL=false(1,M^2);
            closeNonEventL(~eventL)=close_;
            idxs=find(closeNonEventL | eventL);
            N=length(idxs);
        end
        
        function [bbwMat, nearN, nnVec, eigVec]=CreateFast(dns, events)
            bbw=DensityBiased(dns.eventBinIdxs, ...
                events, dns.mins, dns.maxs, dns.M);
            bbw.calcFast;
            bbwMat=reshape(bbw.dnsVec, [dns.M dns.M]);
            nearN=bbw.nearN;
            nnVec=bbw.nnVec;
            eigVec=bbw.eigVec;
        end
        
        function n=MIN_PLOT_SIZE
            n=50;
        end
        
        function bw=ConfirmBandWidth(dns)
            bw=dns.bandwidth;
            if dns.N <=DensityBiased.MIN_PLOT_SIZE
                if bw==-1
                    bw=0;
                end
            end
        end
        
        function [bbwMat, nearN, nnVec, eigVec]=Create(dns, events)
            bbwMat=[];
            nnVec=[];
            eigVec=[];
            nearN=0;
            if dns.bandwidth>=0
                return;
            end
            try
            bbw=CytoGate.Get.getNumeric('biasedBandwidth', 0);
            catc ex
            bbw=0;
            end
            tm=tic;
            nearN=DensityBiased.DFLT_NEIGH(dns.N);
            if bbw==5
                jd=msg('Automatic bandwidth is calculating....', ...
                    [], 'north east++', 'Extra time needed ...');
                [bbwMat, nnVec, eigVec]=DensityBiased.CalcSlow(...
                    events, dns);
            else
                [bbwMat, nearN, nnVec, eigVec]=DensityBiased.CreateFast(...
                    dns, events);
            end
            disp('Cost of biased bandwidth');
            toc(tm);
            if exist('jd', 'var')
                jd.dispose;
            end
        end
        
        function ok=USE_EIG_AS_BACKGROUND
            ok=true;
        end
        
        function [numClusts, eventClusterIds]=ClusterAnalyze(dns,...
                backgroundFactor)
            f=dns.fmatVector;
            stdErr=f/sqrt(dns.nearN);
            %stdErr=sqrt(f.^2 - (f.^2/(this.nearN-1)));
            if DensityBiased.USE_EIG_AS_BACKGROUND
                f2=f;
                f2(dns.eigVec<DensityBiased.DFLT_EIG_BACKGROUND)=0;
                fMat=reshape(f2, dns.M, dns.M);
                back=f2<=backgroundFactor*stdErr; %points designated background
                %fMat=dns.fmat;
                [~, ~, Pointers]=MatBasics.GetClusters(fMat, back, []);
            else
                back=f<=backgroundFactor*stdErr; %points designated background
                [~, ~, Pointers]=MatBasics.GetClusters(dns.fmat, back, []);
            end
            javaDbm=edu.stanford.facs.swing.Dbm(dns.M);
            javaDbm.pointers=Pointers;
            javaDbm.density=f;
            javaDbm.stdErr=stdErr;
            javaDbm.reportChangeCount=false;
            javaDbm.merge;
            javaDbm.fixClusterTear;
            Pointers=javaDbm.pointers';
            dns.pointers=Pointers; %this is to save Pointers for making vector plot later
            peaks=Merge.Peaks(Pointers);
            p=dns.pointers>0;
            while any(p)
                dns.pointers(p)=dns.pointers(dns.pointers(p));  %follows the path of all positive pointers until they all go to dummy states
                p=dns.pointers>0;
            end
            changes=Merge.ForNonConstant(dns, peaks);
            nChanges=length(changes);
            peakClues=zeros(nChanges,1);
            mergeClues=zeros(nChanges, 1);
            mxDns=zeros(nChanges, 1);
            for i=1:nChanges
                change=changes{i};
                peak=change(1);
                merge=change(2);
                mxDns(i)=change(3);
                peakClues(i)=dns.pointers(peak);
                mergeClues(i)=dns.pointers(merge);
            end
            set=java.util.TreeSet;
            for i=1:nChanges
                peakClue=peakClues(i);
                mergeClue=mergeClues(i);
                idx=find(peakClues==mergeClue,1);
                if idx>0
                    set.clear;
                    set.add(idx);
                    dnsty=mxDns(idx);
                    while idx>0
                        mergeClue=mergeClues(idx);
                        idx2=find(peakClues==mergeClue,1);
                        if isempty(idx2) || idx2==idx 
                            break;
                        end
                        if set.contains(idx2)
                            break;
                        end
                        set.add(idx2);
                        mx2=mxDns(idx2);
                        if mx2>dnsty
                            dnsty=mx2;
                            idx=idx2;
                        end
                    end
                    mergeClue=mergeClues(idx);
                end
                l=dns.pointers==peakClue;
                l2=dns.pointers==mergeClue;
                [x1, y1]=ind2sub([dns.M dns.M], find(l));
                [x2, y2]=ind2sub([dns.M dns.M], find(l2));
                DD=min(pdist2([x1;y1]', [x2;y2]', 'euclidean',...
                    'smallest', 1));
                if DD>1
                    %disp('huh');
                else
                    dns.pointers(l)=mergeClue;
                end
            end
            dns.pointers(dns.pointers==-1)=0;  %send background gridpoints to label 0
            eventClusterIds=dns.pointers(dns.eventBinIdxs); %assigns each original data point the dummy cluster number of its closest gridpoint
            clusts=unique(eventClusterIds);
            clusts=clusts(clusts~=0);  %all non-background clusters
            if DensityBiased.DEBUG>0
                for clue=clusts
                    idxs=find(dns.pointers==clue);
                    idxs=idxs(ismember(idxs, dns.eventBinIdxs));
                    N=length(idxs);
                    [mxDns, dnsI]=max(dns.fmatVector(idxs));
                    dnsI=idxs(dnsI);
                    eigDns=dns.eigVec(dnsI);
                    [x,y]=ind2sub([dns.M dns.M], dnsI);
                    fprintf('N=%d \tMAX dns=%s@%d/%d\teig=%s',  ...
                        N, String.encodeRounded(mxDns, 0), x, y, ...
                        String.encodeRounded(eigDns,2));
                    [mxEig, eigI]=max(dns.eigVec(idxs));
                    eigIdx=idxs(eigI);
                    [x,y]=ind2sub([dns.M dns.M], eigIdx);
                    fprintf('\tMAX eig=%s@%d/%d\tMEAN %s\t(%d)\n', ...
                        String.encodeRounded(mxEig, 2), x, y, ...
                        String.encodeRounded(mean(dns.eigVec(idxs)),2), clue); 
                end
            end

            for clue=clusts
                idxs=find(dns.pointers==clue);
                idxs=idxs(ismember(idxs, dns.eventBinIdxs));
                mxEig=max(dns.eigVec(idxs));
                if mxEig<DensityBiased.DFLT_EIG_MIN
                    dns.pointers(idxs)=0;
                end
            end
            dns.pointers(dns.pointers==-1)=0;  %send background gridpoints to label 0
            eventClusterIds=dns.pointers(dns.eventBinIdxs); %assigns each original data point the dummy cluster number of its closest gridpoint
            clusts=unique(eventClusterIds);
            clusts=clusts(clusts~=0);  %all non-background clusters            
            numClusts=length(clusts);
            j=1;
            for i=clusts
                eventClusterIds(eventClusterIds==i)=j;  %relabel clusters from 1 to NumClusts
                dns.pointers(dns.pointers==i)=j;  %relabel pointers
                j=j+1;
            end
            if dns.truncated
                ca=zeros(1, length(dns.onScale));
                ca(dns.onScale)=eventClusterIds;
                eventClusterIds=ca;
            end
            clustersWithoutEventBinIdxs=dns.pointers<-1;
            if sum(clustersWithoutEventBinIdxs)>0
                if Density.IsDebugging
                    fprintf('There are %d grid points with clusters but no events!', sum(clustersWithoutEventBinIdxs));
                end
                dns.pointers(clustersWithoutEventBinIdxs)=0;  %send background gridpoints that had clusters but no events
            end
            dns.computePointersWithEvents(eventClusterIds);
        end
       
        function [xyData, xyMin, xyMax]=NormalizeIfNotLogicle(...
                xyData, xyMin, xyMax)
            for axis=1:2
                if xyMax(axis)>1
                    m=xyMin(axis);
                    r=xyMax(axis)-m;
                    xyData(:,axis)=(xyData(:,axis)-m)/r;
                    xyMin(axis)=0;
                    xyMax(axis)=1;
                end
            end
        end
        
        function [dataIdxs, nDataIdxs]=...
                Search(kdTree, gridPointOnDataScale, withinDistance)
            
            % MatLab function rangesearch returns dataIdxs that 
            %   are sorted in the ascending order of the corresponding 
            %   distance values
            dataIdxs=rangesearch(kdTree, gridPointOnDataScale, withinDistance);
            nDataIdxs=length(dataIdxs);
            if nDataIdxs==1 && iscell(dataIdxs)
                dataIdxs=dataIdxs{1};
                nDataIdxs=length(dataIdxs);
            end
        end
        
        function nearN=DFLT_NEIGH(N)
            %nearN=256;
            %nearN=512;
            %nearN=576;
            nearN=640;
            %nearN=704;
            %nearN=768;
            %nearN=896;
            %nearN=1024;
            if nearN>floor(N/7.825)
                nearN=floor(N/7.825);
                if nearN<10
                    nearN=10;
                end
            end     
        end

        function n=DFLT_MAX_TWIN_PEAK_DIST
            n=1/32;
        end

        function n=DFLT_MAX_D
            n=.08;
        end

        function n=DFLT_EIG_BACKGROUND
            n=.3;
        end

        function n=DFLT_EIG_MIN
            n=.88;
        end
        function n=DFLT_MIN_D
            n=.016;
        end
        
        function [mnDistance, mxDistance]=BandWidth(M, gridPtsOnDataScale)
            mnIdx=sub2ind([M M], 1, 1);
            mxIdx=sub2ind([M M], M, M);
            maxGridD=pdist2(...
                gridPtsOnDataScale(mnIdx, :),gridPtsOnDataScale(mxIdx, :));
            mxDistance=DensityBiased.DFLT_MAX_D*maxGridD;
            mnDistance=DensityBiased.DFLT_MIN_D*maxGridD;
        end
        
        
        function result=WayneCovariance(data)
            %means=mean(data);
            q1=sum(data(:,1).^2);
            q2and3=sum(data(:,2).^2);
            q4=sum( data(:,1).*data(:,2));
            N=size(data, 1);
            result=[q1/N q4/N;q4/N q2and3/N];
        end
        

        function [densityMat, nnVec, eigVec]=CalcSlow(xyData, options)
            xyMin=options.mins; % minimum values of X and Y
            xyMax=options.maxs; % maximum values of X and Y
            M=options.M;
            xyGridIdxs=options.eventBinIdxs; %Each data point's assigned grid point
            nearestNeighbors=DensityBiased.DFLT_NEIGH(options.N);
            nnVec=zeros(1, M^2);
            eigVec=zeros(1, M^2);
            [xyData, xyMin, xyMax]=DensityBiased.NormalizeIfNotLogicle(...
                xyData, xyMin, xyMax);
            kdTree=KDTreeSearcher(xyData);            
            [gridPtsOnDataScale, ptsGridScale]=...
                DensityBiased.ScaleGridPoints(M, xyMin, xyMax);
            [~, mxDistance]=...
                DensityBiased.BandWidth(M, gridPtsOnDataScale);
            densityVec=zeros(1, M^2);
            [gridPoints, N]=DensityBiased.RemoveEmptyPoints(xyGridIdxs, ...
                [ptsGridScale(:,2) ptsGridScale(:,1)], M, 8);
            for idx=1:N
                gridPt=gridPoints(idx);
                gridPtOnDataScale=gridPtsOnDataScale(gridPt,:);
                [dataIdxs, nDataIdxs]=DensityBiased.Search(...
                    kdTree, gridPtOnDataScale, mxDistance);
                if nDataIdxs>nearestNeighbors
                    dataIdxs=dataIdxs(1:nearestNeighbors);
                end
                if length(dataIdxs)<2
                    continue;
                end
                closestEvents=xyData(dataIdxs,:);
                finalData=closestEvents - gridPtOnDataScale;
                covValue=DensityBiased.WayneCovariance(finalData);
                num=1/sqrt(det(covValue));
                if isreal(num) && ~isinf(num) && ~isnan(num)
                    densityVec(gridPt)=num;
                    nnVec(gridPt)=length(dataIdxs);
                    eigValue=eig(covValue);
                    ratioEigValue=eigValue(1)/eigValue(2);
                    eigVec(gridPt)=ratioEigValue;                    
                else
                    % do nothing
                end
            end
            densityMat=reshape(densityVec, [M M]);
        end

        function [densityMat, nnVec, eigVec2]=CalcSlowDebug(xyData, dns, debug)
            events=xyData;
            xyMin=dns.mins;
            xyMax=dns.maxs;
            M=dns.M;
            xyGridIdxs=dns.eventBinIdxs;
            nearestNeighbors=DensityBiased.DFLT_NEIGH(dns.N);
            nnVec=zeros(1, M^2);
            eigVec2=zeros(1, M^2);
            [xyData, xyMin, xyMax]=DensityBiased.NormalizeIfNotLogicle(xyData, xyMin, xyMax);
            kdTree=KDTreeSearcher(xyData);            
            [gridPtsOnDataScale, ptsGridScale]=...
                DensityBiased.ScaleGridPoints(M, xyMin, xyMax);
            [~, mxDistance]=...
                DensityBiased.BandWidth(M, gridPtsOnDataScale);
            densityVec=zeros(1, M^2);
            [gridPoints, N]=DensityBiased.RemoveEmptyPoints(xyGridIdxs, ...
                [ptsGridScale(:,2) ptsGridScale(:,1)], M, 8);
            if debug
                bbw=DensityBiased(xyGridIdxs, xyData, xyMin, xyMax, M);
                bbw.calcFast;
            end
            for idx=1:N
                gridPt=gridPoints(idx);
                gridPtOnDataScale=gridPtsOnDataScale(gridPt,:);
                [dataIdxs, nDataIdxs]=DensityBiased.Search(...
                    kdTree, gridPtOnDataScale, mxDistance);
                if nDataIdxs>nearestNeighbors
                    dataIdxs=dataIdxs(1:nearestNeighbors);
                end
                if length(dataIdxs)<2
                    continue;
                end
                closestEvents=xyData(dataIdxs,:);
                finalData=closestEvents - gridPtOnDataScale;
                covValue=DensityBiased.WayneCovariance(finalData);
                num=1/sqrt(det(covValue));
                [numXx, numYy]=ind2sub([M M], gridPt);
                if numYy==25 && numXx >26 && numXx < 32
                    numStr=String.encodeRounded(num,0);
                    disp(numStr);
                end
                if isreal(num) && ~isinf(num) && ~isnan(num)
                    densityVec(gridPt)=num;
                    nnVec(gridPt)=length(dataIdxs);
                    eigValue=eig(covValue);
                    ratioEigValue=eigValue(1)/eigValue(2);
                    eigVec2(gridPt)=ratioEigValue;
                    
                    unequal=false;
                    if debug
                        s=String.encodeRounded(densityVec(gridPt),4);
                        sVec=String.encodeRounded(bbw.dnsVec(gridPt),4);
                        if ~isequal(sVec, s)
                            fprintf(['%s does not equal %s @idx=%d, %d '...
                                'neighbors\n'], s, sVec, idx, nDataIdxs);
                            unequal=true;
                        end
                        s=String.encodeRounded(eigVec2(gridPt), 4);
                        sVec=String.encodeRounded(bbw.eigVec(gridPt),4);
                        if ~isequal(sVec, s)
                            fprintf(['eig %d does not equal %d->idx=%d, %d '...
                                'neighbors\n'], s, sVec,idx, nDataIdxs);
                            unequal=true;
                        end
                    end
                    if mod(idx+1,500)==0
                        fprintf('%d\n', idx);
                    end
                else
                    [numXx, numYy]=ind2sub([M M], gridPt);
                    str=sprintf('x=%d, y=%d', numXx, numYy);
                    if ~isreal(num)
                        disp(['not real@' str]);
                    elseif isinf(num)
                        disp(['infinitity@' str]);
                    elseif isnan(num)
                        disp(['NAN@' str]);
                    else
                        disp(['SOMETHING wrong @' str]);
                    end
                end
            end
            densityMat=reshape(densityVec, [M M]);
        end


    end
    
end
