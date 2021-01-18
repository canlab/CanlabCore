%   AUTHORSHIP
%   Software Developers:  Rachel Finck, 
%                       Connor Meehan <connor.gw.meehan@gmail.com>
%                       Stephen Meehan <swmeehan@stanford.edu> 
%   
%   ORIGINAL PUBLICATION
%   http://cgworkspace.cytogenie.org/GetDown2/demo/dbm.pdf
%
%   PUBLICATION REVISIONS
%   https://static-content.springer.com/esm/art%3A10.1038%2Fs42003-019-0467-6/MediaObjects/42003_2019_467_MOESM1_ESM.pdf
%
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

classdef Density < handle
    
    properties(Constant)
        D=2;
        CONTOUR_LIMIT=1000;
        CONTOUR_BANDWIDTH=.014;
        DbmBandwidth='DbmBandwidth'; 
        DbmBackgroundType='DbmBackgroundType'; 
        DbmBackgroundFactor='DbmBackgroundFactor'; 
        DbmSlopeIsSignificant='DbmSlopeIsSignificant';
        DbmM='DbmM';
        DbmNmin='DbmNmin';
        DETAILS={'most high', 'very high', 'high', ...
            'medium', 'low', 'very low', 'adaptive', ...
            'nearest neighbor', 'dbscan arguments'};
        
    end
    
    properties(SetAccess=private)
        opts;
        labels;
        clusterNames;
        clusterColors;
        clusterMdns;
        clusterLabels;
        backgroundFactor;
        backgroundType;
        decimals;
        notes;
        gradients;
        backGround;
        maxNbrs;
        isSlopeBig;
        isSlopeSignificant;
        debug=false;
        Z={};
        slopes=[];
        pathEnds=[];
        sigPathEnds=[];
        criticalValue=0;
        Kappa=0;
        M=0;
        N=0;
        N_min;
        mins=[];
        maxs=[];
        deltas=[];
        eventBinIdxs=[];
        fmat=[];
        wmat=[];
        phimat=[];
        rawPointers=[];
        contourH=[];
        h=[];
        L=[];
        onScale=[];
        truncated=false;
        ye=[];
        wmatVector;
        fmatVector;
        fmatVector2;
        nnVec;
        eigVec;
        nnMat; % for tool tips
        eigMat;% for tool tips
        bandWidthInput=[];
        bandwidth=[];
        xgrid;
        ygrid;
        gridAssign;
        pointersWithEvents=[];
        contourZ=[];
        pseudoZ=[];
        contourXm=[];
        contourYm=[];
        contourV=[];
        colorMap=[];
        probabilityContourPercent=0;
        isBiased=false;
        nearN=0;
    end
    
    properties
        probabilityDensity;
        pointers=[]; %with or without events
    end
    
    methods
        function yes=needToRecreate(this, options)
            yes=this.M ~= options.DbmM;
            if ~yes
                yes=this.truncated;
            end
            if ~yes
                if isfield(options,  Density.DbmBandwidth)
                    yes=this.bandWidthInput ~= options.DbmBandwidth;
                end
                if ~yes
                    if isfield(options, Density.DbmNmin)
                        yes=this.N_min~=options.DbmNmin;
                    end
                end
            end
        end
        
        function [displayData, display2Event, displayWeight]=getDisplayData(this, data)
            d=Density.ToData(this.ye(1,:), this.ye(2,:), ...
                this.wmat);
            carriesWeight=d(:,3)>0;
            if nargin==1 
                displayData=d(carriesWeight, 1:2);
                displayWeight=d(carriesWeight, 3);
                if this.truncated
                    e=zeros(length(this.onScale), 1);
                    e(this.onScale)=this.eventBinIdxs;
                else
                    e=this.eventBinIdxs;
                end
                [binIdxsWithEvents, eventIdxsForBin]=unique(e);
                binIdxsWithWeight=find(carriesWeight);
                carriesWeightAndEvents=ismember(binIdxsWithWeight, ...
                    binIdxsWithEvents);
                binIdxsWithWeightAndEvents=binIdxsWithWeight(carriesWeightAndEvents);
                w2e=bsearch(binIdxsWithEvents,binIdxsWithWeightAndEvents);
                display2Event=eventIdxsForBin(w2e);
                displayData=displayData(carriesWeightAndEvents, :);
                displayWeight=displayWeight(carriesWeightAndEvents, :);
            else
                displayData=d(carriesWeight, 1:2);
                displayWeight=d(carriesWeight, 3);
                e=this.eventBinIdxs;
                [binIdxsWithEvents, eventIdxsForBin]=unique(e);
                binIdxsWithWeight=find(carriesWeight);
                mapW2E=bsearch(binIdxsWithEvents,binIdxsWithWeight);
                
                if this.truncated
                    data=data(this.onScale, :);
                end
                assert( issorted( binIdxsWithEvents ) );
                assert( isempty(find(~ismember(...
                    binIdxsWithEvents, binIdxsWithWeight), 1)));
                ranges=this.maxs-this.mins;
                idx2BinIdxsWithWeightAndEvents=find(ismember(...
                    binIdxsWithWeight, binIdxsWithEvents));
                N=length(idx2BinIdxsWithWeightAndEvents);
                
                for i=1:N
                    pick=idx2BinIdxsWithWeightAndEvents(i);
                    %first re-assert the above find/ismember op
                    %to make sure no parity error in RAM, ROM or CPU
                    assert( ismember(binIdxsWithWeight(pick), binIdxsWithEvents));
                    idx=mapW2E(pick);
                    eventIdx=eventIdxsForBin( idx );
                    nearestBinIdx=binIdxsWithEvents(idx);
                    assert( binIdxsWithEvents(idx ) == e(eventIdx) );
                    %equivalent assert without documenting meaning of each
                    %index IS:
                    assert(binIdxsWithEvents(mapW2E(pick) )==...
                        e( eventIdxsForBin( mapW2E(pick))))
                    dataAtEvent=data(eventIdx, :);
                    dataAtBinIdx=displayData(pick, [1:2]);
                    result=abs(dataAtEvent-dataAtBinIdx);
                    %the statement below does not WORK because it produces scalar:
                    % perc=result/ranges*100;
                    perc=result./ranges*100;
                    if max(perc)>.5
                        px=String.encodeRounded(perc(1),3);
                        py=String.encodeRounded(perc(2),3);
                        xBin=String.encodeRounded(dataAtBinIdx(1),4);
                        yBin=String.encodeRounded(dataAtBinIdx(2),4);
                        xEvent=String.encodeRounded(dataAtEvent(1),4);
                        yEvent=String.encodeRounded(dataAtEvent(2),4);
                        [xx,yy]=Density.ToMatrixIdxs(nearestBinIdx, this.M);
                        xLinSpace=String.encodeRounded(this.ye(1, xx),4);
                        yLinSpace=String.encodeRounded(this.ye(2, yy),5);
                        fprintf(['High gap of x=%s%% && y=%s%% between X, Y at bin #%d: [%s, %s]\n'...
                            '   && X/Y at event associated with THIS same bin: [%s, %s]\n'], ...
                            px, py, binIdxsWithEvents(idx ), xBin, yBin, xEvent, yEvent);
                        assert(min(perc)<1, 'Either gap in X or gap in Y MUST be less than .25');
                        
                    end
                end
                
                %some bin idxs have weight without data points because they
                %are close to events that are part of density smoothing
                binIdxsWithWeightButNotEvents=find(~ismember(...
                    binIdxsWithWeight, binIdxsWithEvents));
                N=length(binIdxsWithWeightButNotEvents);
                for i=1:N
                    pick=binIdxsWithWeightButNotEvents(i);
                    %first re-assert the above find/ismember op
                    %to make sure no parity error in RAM, ROM or CPU
                    assert( ~ismember(binIdxsWithWeight(pick), binIdxsWithEvents));
                    idx=mapW2E(pick);
                    eventIdx=eventIdxsForBin( idx );
                    nearestBinIdx=binIdxsWithEvents(idx);
                    assert( nearestBinIdx==e(eventIdx) );
                    dataAtEvent=data(eventIdx, :);
                    dataAtBinIdx=displayData(pick, [1:2]);
                    result=abs(dataAtEvent-dataAtBinIdx);
                    perc=result./ranges*100;
                    if max(perc)>3 && Density.IsDebugging2014
                        actualBinIdx=MatBasics.NthTrueInBoolVector(carriesWeight,pick);
                        px=String.encodeRounded(perc(1),3);
                        py=String.encodeRounded(perc(2),3);
                        xBin=String.encodeRounded(dataAtBinIdx(1),4);
                        yBin=String.encodeRounded(dataAtBinIdx(2),4);
                        xEvent=String.encodeRounded(dataAtEvent(1),4);
                        yEvent=String.encodeRounded(dataAtEvent(2),4);
                        [xx,yy]=Density.ToMatrixIdxs(nearestBinIdx, this.M);
                        xLinSpace=String.encodeRounded(this.ye(1, xx),4);
                        yLinSpace=String.encodeRounded(this.ye(2, yy),4);
                        fprintf(['High gap of x=%s%% && y=%s%% between '...
                            'X, Y at actual bin #%d: [%s, %s] && X/Y at'...
                            ' event associated with nearest '...
                            'bin #%d: [%s, %s]\n'], ...
                            px, py, actualBinIdx, xBin, yBin, ...
                            nearestBinIdx,xEvent, yEvent);
                        if min(perc)>1
                            disp('Either gap in X or gap in Y MUST be less than .25');
                        end
                    end
                end
                if this.truncated
                    e=zeros(length(this.onScale), 1);
                    e(this.onScale)=this.eventBinIdxs;
                else
                    e=this.eventBinIdxs;
                end
                [binIdxsWithEvents, eventIdxsForBin]=unique(e);
                binIdxsWithWeight=find(carriesWeight);
                carriesWeightAndEvents=ismember(binIdxsWithWeight, ...
                    binIdxsWithEvents);
                binIdxsWithWeightAndEvents=binIdxsWithWeight(carriesWeightAndEvents);
                w2e=bsearch(binIdxsWithEvents,binIdxsWithWeightAndEvents);
                display2Event=eventIdxsForBin(w2e);
                displayData=displayData(carriesWeightAndEvents, :);
                displayWeight=displayWeight(carriesWeightAndEvents, :);
            end
        end
        
        function clusterMat=getClusterMat(this)
            clusterMat=reshape(this.pointers, this.M, this.M);
        end
        
        function y=getY(this)
            [ygrid,xgrid]=meshgrid(this.ye(2,:), this.ye(1,:));
            y=[xgrid(:) ygrid(:)]';
        end
        
        function H=contour(this, ax, probabilityContourPercent, ...
                color, lineWidth)
            if ~isempty(this.probabilityDensity)
                H=this.probabilityDensity.drawContours(ax, ...
                    probabilityContourPercent, color, lineWidth);
                return;
            end
            if isempty(this.contourZ)
                this.contourDensity;
            end
            if isempty(this.contourV) || ...
                    this.probabilityContourPercent ~= probabilityContourPercent
                this.probabilityContourPercent=probabilityContourPercent;
                if probabilityContourPercent==2.4
                    contourLevels=8;
                elseif probabilityContourPercent==2.2
                    contourLevels=6;
                else
                    contourLevels=floor(100/probabilityContourPercent);
                end
                MM=this.M^2;
                T=reshape(this.contourZ, 1,MM);
                T=sort(T);
                CT=cumsum(T);
                NT=CT/CT(end);
                this.contourV=zeros(1, contourLevels);
                if probabilityContourPercent==2.4
                    numerators=[ .5/100 1/100 2/100, 4/100, 8/100, 16/100, 32/100, 64/100];
                    for numerator=1:contourLevels
                        idx=bsearch(NT, numerators(numerator));
                        this.contourV(numerator)=T(idx);
                    end
                    elseif probabilityContourPercent==2.2
                    numerators=[ 2/100, 4/100, 8/100, 16/100, 32/100, 64/100];
                    for numerator=1:contourLevels
                        idx=bsearch(NT, numerators(numerator));
                        this.contourV(numerator)=T(idx);
                    end
                else
                    for numerator=1:contourLevels
                        idx=bsearch(NT, numerator/contourLevels);
                        this.contourV(numerator)=T(idx);
                    end
                end
            
            end
            [~,H]=contour(ax, this.contourXm, this.contourYm, ...
                this.contourZ, this.contourV, 'k', 'color', color, ...
                'LineStyle', '-', 'LineWidth', lineWidth);
        end
        
        function contourV=getProbabilityDensities(this, contourLevels)
            if isempty(this.contourZ)
                this.contourDensity;
            end
            MM=this.M^2;
            T=reshape(this.contourZ, 1,MM);
            T=sort(T);
            CT=cumsum(T);
            NT=CT/CT(end);
            contourV=zeros(1, contourLevels);
            for numerator=1:contourLevels
                idx=bsearch(NT, numerator/contourLevels);
                contourV(numerator)=T(idx);
            end
        end
        
        function density2D(this, data, ax, colorRangeStart, colorRangeEnd)
            if ~isempty(this.probabilityDensity)
                if nargin<5
                    this.probabilityDensity.drawJetColors(ax);
                    return;
                else
                    this.probabilityDensity.drawColors(ax, 64, data, ...
                        colorRangeStart, colorRangeEnd);
                    return;
                end
            end
            if isempty(this.fmatVector2)
                this.contourDensity;
            end
            isPseudoColor=nargin<4;
            if size(data, 1)==length(this.onScale)
                data=data(this.onScale, :);
                [x1,~,x3]=unique(this.eventBinIdxs);
            else
                onScale2=MatBasics.FindOnScale(data, this.mins, this.maxs);
                data=data(onScale2, :);
                z=reshape(1:this.M^2,this.M,this.M);
            	eb=interp2(this.xgrid, this.ygrid, z',...
                    data(:,1),data(:,2),'nearest');  %this associates each data point with its nearest grid point
                [x1,~,x3]=unique(eb);
            end
            if isPseudoColor
                colors=jet(128);
                N_=size(colors,1);
                if isempty(this.pseudoZ)
                    try
                        this.pseudoZ=this.getProbabilityDensities(N_);
                    catch ex
                    end
                end
                probabilityRange=this.pseudoZ;
            else
                if nargin<5 
                    colors=colorRangeStart;
                else
                    a1=mean(colorRangeEnd);
                    a2=mean(colorRangeStart);
                    if a1<a2
                        m1=max(colorRangeEnd);
                        if m1>.75
                            f=.75/m1;
                            colorRangeEnd=colorRangeEnd*f;
                        end
                    end
                    N_=32;
                    colors=zeros(N_,3);
                    colors(1,:)=colorRangeStart;
                    colors(N_,:)=colorRangeEnd;
                    gap=zeros(1,3);
                    for i=1:3
                        gap(i)=colorRangeEnd(i)-colorRangeStart(i);
                    end
                    for i=2:N_-1
                        for j=1:3
                            colors(i,j)=colors(1,j)+(i/N_*gap(1,j));
                        end
                    end
                end
                N_=length(colors);
                probabilityRange=this.getProbabilityDensities(N_);
            end
            if size(data,1)<10
                color=colors(1, :);
                plot(ax, data(:,1), data(:,2), 'd',...
                    'markersize', 2, 'MarkerEdgeColor',...
                    color, 'LineStyle', 'none');
                return;
            end
            try
                colormap(colors);
            catch ex
                disp('huh');
            end
            densities=this.fmatVector2(x1);
            lookup=bsearch(probabilityRange,densities);
            eventColors=lookup(x3);
            usedColors=unique(eventColors);
            N2=length(usedColors);
            sz=size(data,1);
            if sz<10000
                marker='d';
                ms=2;
            else
                marker='.';
                ms=2;
            end
            for i=1:N2
                colorIdx=usedColors(i);
                li=eventColors==colorIdx;
                plot(ax, data(li,1), data(li,2), marker,...
                    'markersize', ms, 'MarkerEdgeColor',...
                    colors(colorIdx, :), 'LineStyle', 'none');
            end
            
        end
        
        function [x, y, xIdx, yIdx, mxDns, avgDns]=getPeak(this, clue)
            bi=find(this.pointersWithEvents==clue);
            if isempty(bi)
                x=0;
                y=0;
                xIdx=0;
                yIdx=0;
                mxDns=0;
                avgDns=0;
            else
                [mxDns,mi]=max(this.fmatVector(bi));
                if nargout>5
                    avgDns=median(this.fmatVector(bi));
                end
                maxIdx=bi(mi);
                [xIdx, yIdx]=ind2sub([this.M this.M], maxIdx);
                x=this.mins(1)+((xIdx-1)*this.deltas(1));
                y=this.mins(2)+((yIdx-1)*this.deltas(2));
            end
        end

        function [xIdx, yIdx, mxEig, avgEig]=getPeakEig(this, clue)
            bi=find(this.pointersWithEvents==clue);
            if isempty(bi)
                xIdx=0;
                yIdx=0;
                mxEig=0;
                avgEig=0;
            else
                [mxEig,mi]=max(this.eigVec(bi));
                if nargout>3
                    avgEig=median(this.eigVec(bi));
                end
                maxIdx=bi(mi);
                [xIdx, yIdx]=ind2sub([this.M this.M], maxIdx);
            end
        end
        
        function binScale=toBinScale(this, value, axis)
            binScale=((value-this.mins(axis))./this.deltas(axis))+1;
        end
        
        function [xx, yy]=toData(this, x, y)
            xx=this.mins(1)+((x-1)*this.deltas(1));
            yy=this.mins(2)+((y-1)*this.deltas(2));
        end
        
        function xx=toDataScale(this, value, axis)
            xx=this.mins(axis)+((value-1)*this.deltas(axis));
        end        
        function [xx,yy]=allGridData(this)
            [x,y]=ind2sub([this.M this.M], [1:this.M^2]);
            xx=this.mins(1)+((x-1)*this.deltas(1));
            yy=this.mins(2)+((y-1)*this.deltas(2));
        end
        
        function bins=getNeighborIdxs(this, bin)
            [x,y]=ind2sub([this.M this.M], bin);
            xIdxs=[];
            yIdxs=[];
            for i=-1:1
                nx=x+i;
                if nx>0 && nx<=this.M
                    for j=-1:1
                        if i==0&&j==0
                            continue;
                        end
                        ny=y+j;
                        if ny>0 && ny<=this.M
                            xIdxs(end+1)=nx;
                            yIdxs(end+1)=ny;
                        end
                    end
                end
            end
            bins=sub2ind([this.M this.M], xIdxs, yIdxs);
        end
        
        function f2D=get2ndDerivative(this)
            f=this.fmat;
            M_=this.M;
            fx=zeros(M_, M_);
            fy=zeros(M_, M_);
            D1=this.deltas(1);
            D2=this.deltas(2);
            for x=2:M_-1
                for y=2:M_-1
                    fx(x,y)=f(x-1,y)+f(x+1,y);
                    fy(x,y)=f(x,y-1)+f(x,y+1);
                end
            end
            x=1;
            %left side minus corners
            for y=2:M_-1
                fx(x,y)=f(x+1,y-1)+f(x+1,y);
                fy(x,y)=f(x,y-1)+f(x,y+1);
            end
            
            %right side minus corners
            x=M_;
            for y=2:M_-1
                fx(x,y)=f(x-1,y)+f(x-1,y-1);
                fy(x,y)=f(x,y-1)+f(x,y+1);                
            end

            %bottom side minus corners
            y=M_;
            for x=2:M_-1
                fx(x,y)=f(x-1,y)+f(x+1,y);
                fy(x,y)=f(x,y-1)+f(x-1,y-1);                
            end
            
            %top side minus corners
            y=1;
            for x=2:M_-1
                fx(x,y)=f(x-1,y)+f(x+1,y);
                fy(x,y)=f(x+1,y+1)+f(x,y+1);                
            end

            %top left
            fx(1,1)=f(2,1)+f(2,2);
            fy(1,1)=f(2,2)+f(1,2);

            %bottom left
            fx(1,M_)=f(2,M_)+f(2,M_-1);
            fy(1,M_)=f(1,M_-1)+f(2, M_-1);

            %top right
            fx(M_,1)=f(M_-1,1)+f(M_-1,2);
            fy(M_,1)=f(M_,2)+f(M_-1, 2);
            %bottom right
            fx(M_,M_)=f(M_-1,M_)+f(M_-1,M_-1);
            fy(M_,M_)=f(M_,M_-1)+f(M_-1, M_-1);
            
            f2D=(fx-2*f)/D1 + (fy-2*f)/D2;

        end

        function numClusts=setLabels(this, labels, ...
                clusterNames, clusterColors, clusterMdns, clusterLabels)
            %this.onScale and this.eventBinIdxs defined by Constructor
            %must create pointers and pointersWithEvents
            assert(length(labels)==length(this.onScale));
            %labels=labels(this.onScale);
            u=unique(labels);
            hasZero=any(labels==0);
            if hasZero
                numClusts=length(u)-1;
            else
                numClusts=length(u);
            end
            this.labels=labels;
            this.clusterNames=clusterNames;
            this.clusterColors=clusterColors;
            this.clusterMdns=clusterMdns;
            this.clusterLabels=clusterLabels;
        end
    end
    methods(Static )
        function [density, weight, idxs, bCoords]=Get3D(data, nBins, bandWidth)
            if nargin<2
                nBins=64;
            end
            if nargin<3
                bandWidth=floor(nBins/5);
                if mod(bandWidth, 2)==0
                    bandWidth=bandWidth+1;
                end
            end
            x=data(:,1);
            y=data(:,2);
            z=data(:,3);
            N=numel(x);
            xBins=linspace(min(x),max(x), nBins);
            yBins=linspace(min(y),max(y), nBins);
            zBins=linspace(min(z),max(z), nBins);
            weight=zeros(nBins, nBins, nBins);
            bCoords=[xBins' yBins' zBins'];
            if nargout>1
                idxs=zeros(N,3);
                for i=1:N
                    xi=find((x(i)>=xBins), 1, 'last');
                    yi=find((y(i)>=yBins), 1, 'last');
                    zi=find((z(i)>=zBins), 1, 'last');
                    weight(xi, yi, zi)=weight(xi, yi, zi)+1;
                    idxs(i,:)=[xi, yi, zi];
                end
            else
                for i=1:N
                    xi=find((x(i)>=xBins), 1, 'last');
                    yi=find((y(i)>=yBins), 1, 'last');
                    zi=find((z(i)>=zBins), 1, 'last');
                    weight(xi, yi, zi)=weight(xi, yi, zi)+1;
                end
            end
            density=smooth3(weight, 'gaussian', bandWidth, .8);
        end
        
        function [ok, isMatLabVersion]=HasDbScan(askForDownloadIfNeeded)
            if nargin<1
                askForDownloadIfNeeded=true;
            end
            clusters=[];
            data=randi(100, 50, 3);
            try
                clusters=dbscan(data, .5, 15);
                isMatLabVersion=true;
            catch ex
                try
                    isMatLabVersion=false;
                    % see if DBSCAN is available
                        clusters=DBSCAN(data, .5, 15);
                    %YES
                catch ex
                    disp('No dbscan ... before r2019a');
                    if askForDownloadIfNeeded
                        if ~isdeployed
                            Density.DownloadDbScan;
                        end
                    end
                end
            end
            ok=~isempty(clusters);
        end
        
        function DownloadDbScan(where)
            if nargin<1
                where='south';
            end
            questDlg(struct('where', where, 'msg', Html.Wrap(...
                    ['Neither implementation of dbscan were found<ol>'...
                    '<li>MatLab''s dbscan (fast) requires r2019a or greater'...
                    '<li>DBSCAN (slower) from the MathWorks File Exchange'...
                    '</ol><br>' Html.WrapSmallTags([...
                    '(Click <u>Download</u> button  to get DBSCAN ' ...
                    'from MathWorks File Exchange)']) '<hr>' ]), ...
                    'checkFnc', @(jd, answer)download(jd, answer), ...
                    'modal', false, 'pauseSecs', 0), 'Get DBSCAN?', ...
                    'Download', 'Ok', 'Download');
            function ok=download(jd, answer)
                if isequal('Download', answer)
                    web(['https://www.mathworks.com/matlabcentral/'...
                        'fileexchange/52905-dbscan-clustering-algorithm']);
                    msg(Html.WrapHr(['After downloading & installing be'...
                        '<br>sure to call addpath to the DB_SCAN folder']));
                    msg(Html.WrapHr(['You must restart MatLab <br>'...
                        'after installing for this to take effect!!']), 0);
                end
                ok=true;
                jd.dispose;
            end
        end
    
        function [epsilon, minpts]=GetDbscanParameters(detail, epsilon, minpts)
            if ~strcmpi(detail, 'dbscan arguments')
                if strcmpi(detail, 'very low')
                    epsilon=3;
                    minpts=15;
                elseif strcmpi(detail, 'low')
                    epsilon=2;
                    minpts=15;
                elseif strcmpi(detail, 'medium')
                    epsilon=1.5;
                    minpts=15;
                elseif strcmpi(detail, 'high')
                    epsilon=1;
                    minpts=15;
                elseif strcmpi(detail, 'very high')
                    epsilon=.66;
                    minpts=5;
                elseif strcmpi(detail, 'most high')
                    epsilon=.5;
                    minpts=4;
                end
            end    
        end
        function [numClusters, clusterIds, density]=FindClusters(...
                data, detail, method2D, pu, epsilon, minpts, ...
                distance, mins, maxs)
            if nargin<9
                maxs=[];
                mins=[];
                if nargin<7
                    distance='euclidean';
                    if nargin<6
                        epsilon=1;
                        if nargin<5
                            minpts=5;
                            if nargin<5
                                pu=[];
                                if nargin<3
                                    method2D='dbm';
                                    if nargin<2
                                        detail='medium';
                                    end
                                end
                            end
                        end
                    end
                end
            end
            wantsDbscan=size(data,2)>2 || strcmpi(method2D, 'dbscan');
            app=BasicMap.Global;
            canNotDo=app.noDbscan || ((verLessThan('matlab', '9.6') ... 
                && isdeployed)); %can't add DBSCAN to compiled runtime;
            if wantsDbscan && ~canNotDo
                if isa(pu, 'PopUp')
                    pu2=pu;
                    priorText=pu.label.getText;
                    pu.label.setText(['<html>Finding "' detail ...
                        '" clusters with dbscan<hr></html>']);
                elseif (islogical(pu) || isnumeric(pu)) && ~isempty(pu) && pu
                    pu2=PopUp(['Finding "' detail ...
                        '" clusters with dbscan'], 'north', ...
                        'Clustering...', false);
                else
                    pu2=[];
                end
                try
                    if nargout>2
                        if isempty(mins)
                            mins=min(data);  %vector of min of each column of data
                            maxs=max(data);   %vector of max of each column of data
                        end
                        density=Density.New(data(:,1:2), detail, ...
                            mins(1:2), maxs(1:2));
                    end
                    [epsilon, minpts]=Density.GetDbscanParameters(...
                        detail, epsilon, minpts);
                    if isempty(distance)
                        clusterIds=dbscan(data, epsilon, minpts);
                    else
                        clusterIds=dbscan(data, epsilon, minpts, ...
                            'Distance', distance);
                    end
                catch ex
                    try
                        if ~isempty(pu2)
                            if ~isa(pu, 'PopUp')
                                pu2.setText(Html.WrapHr(['MATLAB''s dbscan not '...
                                    'available ... <br>using DBSCAN from'...
                                    'MathWorks File Exchange.<br><br><i>Note that'...
                                    ' DBSCAN is quite slow!</i>']));
                            else
                                pu.label.setText('Using 3rd party DBSCAN ...');
                            end
                        end
                        clusterIds=DBSCAN(data, epsilon, minpts);
                    catch ex
                        disp(ex);
                        clusterIds=[];
                        if ~isdeployed && ~isempty(pu) && ~isempty(pu2)
                            msg(Html.WrapHr(['Density:FindClusters '...
                                'needs dbscan <br>in MatLab r2019a or'...
                                ' later<br>or a <b>download</b> of DBSCAN '...
                                '<br>from MathWorks FileExchange<br>'...
                                '<br>2D clustering will be done<br>'...
                                'instead with loss of accuracy...']), ...
                            8, 'north=+');
                            Density.DownloadDbScan;
                        end
                        app.noDbscan=true; 
                        canNotDo=true;
                    end
                end
                if ~isempty(pu2)
                    if ~isa(pu, 'PopUp')
                        pu2.close;
                    else
                        pu.label.setText(priorText);
                    end
                end
                numClusters=sum(unique(clusterIds)>0);                
            end
            if wantsDbscan && canNotDo
                wantsDbscan=false;
                warning('WARNING Shaving data to 2D....no dbscan or DBSCAN found');
                data=data(:,1:2);
                if ~isempty(mins)
                    mins=mins(1:2);
                end
                if ~isempty(maxs)
                    maxs=maxs(1:2);
                end
            end
            if ~wantsDbscan 
                if strcmpi(detail, 'low')
                    [numClusters,clusterIds, density]=...
                        Density.ClusterLow(data, mins, maxs);
                elseif strcmpi(detail, 'most high')
                    [numClusters, clusterIds, density]=...
                        Density.ClusterMostHigh(data, mins, maxs);
                elseif strcmpi(detail, 'very high')
                    [numClusters, clusterIds, density]=...
                        Density.ClusterVeryHigh(data, mins, maxs);
                elseif strcmpi(detail, 'high')
                    [numClusters,clusterIds, density]=...
                        Density.ClusterHigh(data, mins, maxs);
                elseif strcmpi(detail, 'medium')
                    [numClusters,clusterIds, density]=...
                        Density.ClusterMedium(data, mins, maxs);
                elseif strcmpi(detail, 'very low')
                    [numClusters,clusterIds, density]=...
                        Density.ClusterVeryLow(data, mins, maxs);
                elseif strcmpi(detail, 'adaptive')
                    [numClusters,clusterIds, density]=...
                        Density.ClusterAdaptive(data, mins, maxs);
                elseif strcmpi(detail, 'nearest neighbor')
                    [numClusters,clusterIds, density]=...
                        Density.ClusterNearestN(data, mins, maxs);
                else
                    warning([detail ' is not a DBM method, using medium']);
                    [numClusters,clusterIds, density]=...
                        Density.ClusterMedium(data, mins, maxs);
                end
            end
        end
        
        function [numClusts, clusterIds, density]=ClusterMostHigh(...
                xy, mins, maxs)
            if nargin<3
                maxs=[];
                mins=[];
            end
            [numClusts, clusterIds, density]=Density.Cluster(...
                xy, 0.95, 256, 2, .1, false, mins, maxs);
        end

        
        function [numClusts, clusterIds, density]=ClusterVeryHigh(...
                xy, mins, maxs)
            if nargin<3
                maxs=[];
                mins=[];
            end
            [numClusts, clusterIds, density]=Density.Cluster(...
                xy, 1.5, 256, 2, .1, false, mins, maxs);
        end

        function [numClusts, clusterIds, density]=ClusterHigh(...
                xy, mins, maxs)
            if nargin<3
                maxs=[];
                mins=[];
            end
            [numClusts, clusterIds, density]=Density.Cluster(...
                xy, 1.6, 256, 1, 4.3, false, mins, maxs);
        end

        function [numClusts, clusterIds, density]=ClusterMedium(...
                xy, mins, maxs)
            if nargin<3
                maxs=[];
                mins=[];
            end
            [numClusts, clusterIds, density]=Density.Cluster(...
                xy, 2.3, 256, 1, 4.3, false, mins, maxs);
        end

        function [numClusts, clusterIds, density]=ClusterLow(...
                xy, mins, maxs)
            if nargin<3
                maxs=[];
                mins=[];
            end
            [numClusts, clusterIds, density]=Density.Cluster(...
                xy, 3.2, 256, 1, 4, false, mins, maxs);
        end

        function [numClusts, clusterIds, density]=ClusterVeryLow(...
                xy, mins, maxs)
            if nargin<3
                maxs=[];
                mins=[];
            end
            [numClusts, clusterIds, density]=Density.Cluster(...
                xy, 4, 256, 2, 1, true, mins, maxs);
        end
        
        %This clusters EXACTLY as described in 2009 publication
        %http://cgworkspace.cytogenie.org/GetDown2/demo/dbm.pdf
        function [numClusts, clusterIds, density]=ClusterAdaptive(...
                xy, mins, maxs)
            if nargin<3
                maxs=[];
                mins=[];
            end
            [numClusts, clusterIds, density]=Density.Cluster(...
                xy, 0, 256, 4, 4.3, false, mins, maxs);
        end
        
        function [numClusts, clusterIds, density]=ClusterNearestN(...
                xy, mins, maxs)
            if nargin<3
                maxs=[];
                mins=[];
            end
            [numClusts, clusterIds, density]=Density.Cluster(...
                xy, -1, 128, 4, 4.3, false, mins, maxs);
        end
        
        function [numClusts, clusterIds, density]=Cluster(xy, ...
                bandWidth, M, backgroundType, backgroundFactor, ...
                slopeSignificance, mins, maxs)
            if nargin<8
                maxs=[];
                if nargin<7
                    mins=[];
                    if nargin<6
                        slopeSignificance=false;
                        if nargin<5
                            backgroundFactor=4.3;
                            if nargin<4
                                backgroundType=1;
                                if nargin<3
                                    M=256;
                                    if nargin<2
                                        bandWidth=2.3;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            options.hWait=0;
            options.mergeSmallNeighboringClusters=1;
            options.DbmMsncPerc= 0;
            options.DbmMsncDist=2;
            options.DbmBandwidth=bandWidth;
            options.DbmM=M;
            options.DbmBackgroundFactor=backgroundFactor;
            options.DbmSlopeIsSignificant=slopeSignificance;
            options.DbmBackgroundType=backgroundType;
            options.DbmNmin=5000;
            if isempty(mins) || isempty(maxs)
                density=Density.Create(xy,options);
            else
                density=Density.CreateWithScale(xy, options,mins,maxs);
            end
            [numClusts, clusterIds]=density.clusterAnalyze(options);
        end
    end
    
    methods
        function [numClusts, eventClusterIds]=clusterAnalyze(this, options)
            this.opts=options;
            debugNewJavaMerging=false;
            if debugNewJavaMerging>0
                startTime=tic;
            end
            if isfield(options, Density.DbmBackgroundFactor)
               bgFactor=options.DbmBackgroundFactor;
            else
                bgFactor=4.3;  %background threshold parameter is NOT 4.3 like paper!
            end
            if this.isBiased
                this.backgroundType=-1;
                this.isSlopeSignificant=0;
                [numClusts, eventClusterIds]=DensityBiased.ClusterAnalyze(this, bgFactor);
                return;
            end
            if isfield(options, Density.DbmSlopeIsSignificant)
                slopeIsSignificant=options.DbmSlopeIsSignificant;
            else
                slopeIsSignificant=0;  %pointer-assignment threshold parameter
            end
            M_=this.M;
            MM=M_^2;

            if ~isfield(options,'hWait')
                title=sprintf( 'Finding clusters for %s cells', ...
                    String.encodeNumber(this.N));
                options.hWait=Gui.ProgressBar(title);
                options.percentDone=0;
                options.totalPercent=1;
            elseif options.hWait==0
                options=rmfield(options, 'hWait');
            end
            
            %PUBREF=STEP 1
            sigmat = 1/(this.N*(this.N-1))*...
                convn(this.wmat,this.phimat.^2,'same') - ...
                1/(this.N-1)*this.fmat.^2;  %d-dim matrix of standard error of estimated densities
            f=this.fmatVector;  %f in single vector
            sig=reshape(sigmat,1,MM);  %sigmat in single vector
            if isfield(options,'hWait')
                waitbar2a(options.percentDone+(.08*options.totalPercent), options.hWait, 'Densities & stderr ');
            end
            
            dfmat=cell(1,Density.D);
            Delta=this.deltas;
            for i=1:Density.D
                PhimatD = -Delta(i)/this.h(i)^2*this.L{i}'.*this.phimat;
                %PUBREF=STEP 2.C
                dfmat{i}=1/this.N*convn(this.wmat,PhimatD,'same');
            end
            
            w=this.wmatVector;
            if ~isfield(options, Density.DbmBackgroundType)
                bgType = 2;
                bgFactor=1;
            else
                bgType=options.DbmBackgroundType;
            end
            %alpha=-1;
            stdErr=sqrt(sig);
            if ~all(this.maxs<=1) && bgType==1 
                bgFactor=4.3^2;
                bgType=4;
            end
            if bgType==4
                back = f <= bgFactor*stdErr; %points designated background
                nonBack = f > bgFactor*stdErr; %points not in the background
            elseif bgType==3  % no scatter
                %secondDerivative=0-bg_thresh;
                %silvermanSigApprox=f/(this.N*this.h(1)*this.h(2)*4*pi);
                secondDerivative=reshape(this.get2ndDerivative,1, MM);
                rightSide=bgFactor*stdErr+(.5*this.h(1)*this.h(2)*secondDerivative);
                back=f<=rightSide;
                nonBack=f>rightSide;
                
                if Density.IsDebugging
                    sigma43=bgFactor*stdErr;
                    halfH1H2D=.5*this.h(1)*this.h(2);
                    
                    vIdx=sub2ind([this.M this.M], 69, 15);
                    fprintf(['@69,15 f=%f AND %f*sigma^=%f AND '...
                        '0.5 h1 h2 (..)=%f! (2nd derivative=%f,  '...
                        'whole rightSide=%f)\n'], f(vIdx), bgFactor, sigma43(vIdx), ...
                        halfH1H2D, secondDerivative(vIdx), rightSide(vIdx));
                    vIdx=sub2ind([this.M this.M], 37, 22);
                    fprintf(['@69,15 f=%f AND %f*sigma^=%f AND '...
                        '0.5 h1 h2 (..)=%f! (2nd derivative=%f,  '...
                        'whole rightSide=%f)\n'], f(vIdx), bgFactor, sigma43(vIdx), ...
                        halfH1H2D, secondDerivative(vIdx), rightSide(vIdx));
                end
            elseif bgType==1  % no scatter
                g=1/this.N*convn(this.wmat,this.phimat.^2,'same');% (normal density)^2 divided by plot frequency
                g=reshape(g,1,MM);
                K2=(bgFactor^2);%backgroundFactor default is 4.3
                leftSide=(this.N*(f.^3))+((K2-1)*(f.^2));
                rightSide=K2*g;
                back=leftSide<=rightSide;
                nonBack=leftSide>rightSide;                
            elseif bgType==2
                if bgFactor==0
                    nonBack=f>0;
                    back=[];
                else
                    [~, inds] = sort(f,'descend');
                    prop = cumsum(w(inds))/sum(w);
                    perc=bgFactor/100;
                    lastPoint=find(prop >= 1 - perc, 1);
                    backInds=inds(1:lastPoint);
                    nonBack=zeros(1,length(f)); 
                    nonBack(backInds)=true(1, length(backInds)); 
                    nonBack = logical(nonBack);
                    back = ~nonBack;
                end
            end
            this.backgroundFactor=bgFactor;
            this.backgroundType=bgType;
            
            %PUBREF = STEP 2.E
            kappa = sum(nonBack)*sum(w(nonBack))/(2*pi*this.N*prod(this.h)*sum(f(nonBack)));
            %norminv is the paper's q(x) critical value.  q(.5) would be 0
            %because half the area of the curve is to the left of 0
            
            Critval=norminv([0 0.95^(1/kappa)],0,1); %norminv(x,mu,sigma)mu=mean,sigma=standard deviation
            critval=Critval(2);  %q(0.095^(1/kappa))
            %% compute Pointers STEP 2 ... temporary version
            %The F matrix is STRICTLY to find the most dense of 8 neighbors
            %in the 2D neighborhood by having 81 MxM lattices that are
            %rotated so that the statement max(F,[],3) sees ONLY
            %each lattice point's 8 neighbors plus self which is 0
            %On XY plot 1-9 indices mean:  
            %               1=northEast 2=north 3=northWest
            %               4=east      5=self  6=west
            %               7=southEast 8=south 9=southWest 
            F=zeros([M_ M_ 9]);
            P1=reshape(1:MM,[M_ M_]);
            P=zeros([M_ M_ 9]);
            Ind=1;
            for i=-1:1
                for j=-1:1
                    F(:,:,Ind)=circshift(this.fmat,[j i]);  %matching up densities of neighbors
                    P(:,:,Ind)=circshift(P1,[j i]);  %matching up corresponding indices of neighbors
                    Ind=Ind+1;
                end
            end
            
            F(end,:,[1 4 7])=nan; %don't let density at one edge switch to opposite edge
            F(:,end,[1 2 3])=nan;
            F(1,:,[3 6 9])=nan;
            F(:,1,[7 8 9])=nan;
            
            %% (07/2011 RJ) Normalize by length of distance vector to actually select steepest gradient (not in paper)
            Fcenter = cat(3, F(:,:,5), F(:,:,5), F(:,:,5), F(:,:,5));
            F(:,:,[1 3 7 9]) = Fcenter + (F(:,:,[1 3 7 9]) - Fcenter)/ sqrt(2);
            %% 
            %PUBREF=STEP 2.A part 1
            [~,maxNeigh]=max(F,[],3); %find 1-9 neighbor index of where the max density is
            %vvecs are directions using e the Euclidean norm.  Positive 
            %direction means moving east on X axis or north on Y.
            %[northEast, north, northWest,
            %   east,self,west,...
            %   southEast,south,southWest]
            Dnorm=sqrt(sum(Delta.^2));
            vvecs{1}=[Delta(1)/Dnorm 0 -Delta(1)/Dnorm ...
                1 0 -1 ...
                Delta(1)/Dnorm 0 -Delta(1)/Dnorm]; %unit vectors in rectangular grid
            vvecs{2}=[Delta(2)/Dnorm 1 Delta(2)/Dnorm ...
                0 0 0 ...
                -Delta(2)/Dnorm -1 -Delta(2)/Dnorm];
            
            xu=vvecs{1}(maxNeigh);
            yu=vvecs{2}(maxNeigh);
            
            %PUBREF = STEP 2.B
            %positive dfmat increasing in the positive (AKA southEast) direction
            %negative dfmat decreasing in the positive (AKA southEast) direction
            Slope=xu.*dfmat{1} + yu.*dfmat{2}; %df in direction of greatest increasing density
            if slopeIsSignificant ~= 0
                d=2;                
                Amat=cell(d,d);
                for i=1:d
                    for j=1:d
                        PhimatA=Delta(i)/this.h(i)^2*Delta(j)/this.h(j)^2*this.L{i}'.*this.L{j}'.*this.phimat.^2;
                        %PUBREF=STEP 2.G
                        Amat{i,j}=1/this.N*convn(this.wmat,PhimatA,'same');
                    end
                end
                %Connor this is not the euclidian norm, but instead a
                %rectangular grid?
                uvecs{1}=[1/sqrt(2) 0 -1/sqrt(2) 1 0 -1 1/sqrt(2) 0 -1/sqrt(2)]; %unit vectors in square grid
                uvecs{2}=[1/sqrt(2) 1 1/sqrt(2) 0 0 0 -1/sqrt(2) -1 -1/sqrt(2)];
                
                %compute Sigma and lambda in direction of greatest increasing density
                Sigma=zeros(this.M, this.M);
                for i=1:2
                    for j=1:2
                        %PUBREF=STEP 2.F (Note that 1/(n - 1) is left out and then replaced in STEP 2.D)
                        Sigma=Sigma+uvecs{i}(maxNeigh).*uvecs{j}(maxNeigh).*(Amat{i,j} - dfmat{i}.*dfmat{j});
                    end
                end
                %PUBREF=STEP 2.D
                lambda=critval*sqrt(1/(this.N-1)*Sigma);
                if this.debug
                    this.gradients=lambda;
                end
                %lambda=critval*sqrt(1./(this.N*this.fmat-1).*Sigma);
                %slopeIsSignificant=1.8;
                bigSlopes=Slope > slopeIsSignificant*lambda;  %which gridpoints have sufficient increase to create a pointer
            else
                bigSlopes=Slope>0;
            end
            this.isSlopeSignificant=slopeIsSignificant;
            Pointmat=zeros(M_,M_);
            for i=1:9
                if i~=5
                    %PUBREF=STEP 2.A part 2
                    a=bigSlopes & maxNeigh==i;  %finds max neighbors in ith direction that have 
                    thisP=P(:,:,i); %sufficiently large slope grid points indices in ith direction
                    Pointmat(a)=thisP(a);  %makes association pointers
                    %disp( sprintf('neighbor# %d, %d, %d!', i, sum(a(:)), sum(thisP(:))));
                end
            end
            
            Pointers=reshape(Pointmat,[1 MM]);
            Pointers(back)=-1;  %sets pointers of background gridpoints to -1
            
            %% compute Pointers STEP 3
            if isfield(options,'hWait')
                waitbar2a(options.percentDone+(.16*options.totalPercent), options.hWait, 'Significant densities ');
            end
            
            %PUBREF=STEP 3
            unassigned=find(Pointers==0); %the indices of gridpoints with no association pointers
            pointed_to=Pointers(Pointers>0); %the indices of gridpoints with association pointers into them
            pathends=intersect(unassigned,pointed_to);  %the indices of pathends: they have pointers into them but no pointers out of them
            sigdens=f(pathends) >= critval*stdErr(pathends);  %these pathends' densities are significant
            
            Pointers(pathends(sigdens))= -1 - pathends(sigdens);  %assign these pathends pointers to dummy state representing a cluster
            insig_pathends=pathends(~sigdens);  %pathends that have insignificant densities            
            insig_paths=false(1,MM);
            insig_paths(insig_pathends)=true;  %keeps track of gridpoints that are on paths to insig_pathends
            
            if isfield(options,'hWait')
                waitbar2a(options.percentDone+(.19*options.totalPercent), options.hWait, 'Starting gradient climbs ');
            end
            
            while ~isempty(insig_pathends)
                to_insigs=ismember(Pointers,insig_pathends);  %finds gridpoints that point into insig_pathends
                insig_pathends=find(to_insigs);  %updates insig_pathends to be the gridpoints that point into them
                insig_paths(to_insigs)=true;  %adds gridpoints that are on paths to insig_pathends
            end
            
            Pointers(insig_paths)=-1;  %sends to background all gridpoints on path to an insignificant pathend
            peaks=sum(sigdens);
            if this.debug
                this.notes=sprintf(...
                    ['thresh=%d, #unassigned=%d, #pointed_to=%d, \n'...
                    '#peaks=%d, #signifPeaks=%d, kappa=%0.6f, critval=%0.6f,\n'...
                    '    #bigSlopes=%d of %d, slopeSum=%0.5f,'...
                    'neighSum=%0.0f, ptrSum1=%0.0f, ptrSum2=%0.0f\n'...
                    '    h=%0.3f/%0.3f Z=%d/%d\n'], ...
                    slopeIsSignificant, length(unassigned), length(pointed_to), ...
                    length(pathends), peaks, kappa, critval, ...
                    sum(bigSlopes(:)), sum(Slope(:)>0), sum(Slope(:)), ...
                    sum(maxNeigh(:)), sum(P(:)), sum(Pointers), ...
                    this.h(1), this.h(2), this.Z(1), this.Z(2));
                disp(this.notes);
                this.pathEnds=pathends;
                this.sigPathEnds=pathends(sigdens);
                this.slopes=Slope;
                this.isSlopeBig=bigSlopes;
                this.maxNbrs=maxNeigh;
                this.Kappa=kappa;
                this.criticalValue=critval;
                this.backGround=reshape(back, this.M, this.M);
            end
            outerLoop=0;
            neighborhood=cell(1, MM);
            possibleClusterTears=[];
            pu=[];
            changes=1;
            if debugNewJavaMerging>0
                toc(startTime);
                fprintf('MatLab density calculations completed\n\n');
                dlmwrite('density.txt', f, 'delimiter', ',', 'precision', 16);
                dlmwrite('pointers.txt', Pointers,'delimiter',  ',','precision', 16);
                dlmwrite('stdErr.txt', stdErr, 'delimiter', ',','precision', 16);
                javaMerging=tic;
            end
            noJava=false;
            try
                javaDbm=edu.stanford.facs.swing.Dbm(M_);
                javaDbm.pointers=Pointers;
                javaDbm.density=f;
                javaDbm.stdErr=stdErr;
                if debugNewJavaMerging>0
                    javaDbm.debugging=1;
                else
                    javaDbm.reportChangeCount=false;
                end
                javaDbm.merge;
            catch ex
                javaDbm=[];
                noJava=true;
                debugNewJavaMerging=false;
                javaMerging=tic;
                startTime=tic;
            end
            if noJava || debugNewJavaMerging>0
                if debugNewJavaMerging>0
                    save('densityBasedMerge', 'MM', 'Pointers', 'f', 'stdErr')
                end 
                toc(javaMerging);
                fprintf('JAVA merging completed\n\n');
                matLabMerging=tic;
                while outerLoop==0 || changes>0 %PUBREF=STEP 5 handle necessary repetitions
                    outerLoop=outerLoop+1;
                    %PUBREF=STEP 4
                    pu=reportClusterMerging(outerLoop, changes, peaks, pu, startTime);
                    newPointers=Pointers;
                    Dummies = find(Pointers < -1); %indices of gridpoints with pointers to dummy states representing clusters
                    fDummies = f(Dummies);  %densities at gridpoints that have pointers to dummy states
                    numDummies = length(Dummies);  %number of dummy states
                    
                    [~,ix]=sort(fDummies,'descend');
                    newDummies=Dummies(ix);  %indices of gridpoints with pointers to dummy states in order of decreasing density
                    innerLoop=0;
                    for i = 1:numDummies
                        innerLoop=innerLoop+1;
                        %make A
                        A = newDummies(i);
                        test=(f(newDummies(i)) - stdErr(newDummies(i)));
                        if Pointers(A)>0
                            % Rachel originally commented:
                            %I can't remember why this is here because
                            %everything in newDummies should have pointers
                            %to dummy states, but I am too lazy to
                            %remove it, assuming I put it here for a%
                            %reason initially.
                            
                            % Answer to Rachel's comment is that the new
                            % merging she did on 6_22_10 WILL cause this when
                            % if ~isempty(starts) && any(f_starts>maxQ)
                            
                            continue
                        end
                        sizeOfA = 0;
                        for k=1:MM
                            newSizeOfA=length(A);
                            if newSizeOfA > sizeOfA  %if not all of A has had its neighbors checked yet
                                if k==newSizeOfA %if checking the neighbors of the current last member of A, update sizeofA
                                    sizeOfA=k;
                                end
                                neighbors=neighborhood{A(k)};
                                if isempty(neighbors)
                                    neighbors=this.getNeighborIdxs(A(k));
                                    neighborhood{A(k)}=neighbors;
                                end
                                A=[A neighbors(Pointers(neighbors)==0 & f(neighbors)+stdErr(neighbors)>test)];
                                [~, whereInA]=unique(A, 'first');  %list of indices of unique values in AA
                                A=A(sort(whereInA));  %unique values of A in order in which they originally appeared, so we don't check the same things twice
                            else
                                break  %if all of A was checked and nothing got added to AA during the last iteration, move on
                            end
                        end
                        
                        %make B
                        neighborsOfAandA=unique([cell2mat(neighborhood(A)) A]);
                        B=neighborsOfAandA(Pointers(neighborsOfAandA)<-1);
                        [fQ,yQ] = max(f(B));
                        newPeak=Pointers(B(yQ));
                        tearAble=Pointers(A(Pointers(A)<-1));
                        useOldMerging=true;
                        neighborsOfMaxB=neighborhood{B(yQ)};
                        starts=Pointers(neighborsOfMaxB)>0;
                        if ~isempty(starts)
                            C=neighborsOfMaxB(starts);
                            f_starts=f(C);
                            if any(f_starts>fQ)
                                [~,maxf]=max(f_starts);
                                newPeak=Pointers(C(maxf));
                                useOldMerging=false;
                            end
                        end
                        
                        if debugNewJavaMerging>0
                            %A(1)=[];  %remove newDummies(i) from A
                            assert(~isempty(find(B==A(1), 1)));
                            if useOldMerging
                                %B(yQ)=[];  %remove q from B
                                debug=B(yQ);
                                assert(Pointers(debug)==newPeak);
                            end
                        end
                        
                        Pointers(A)=newPeak; %create pointers from all points in AA to q: changed to Pointers(q) from q
                        Pointers(B)=newPeak; %create pointers from all points in B to q:  changed to Pointers(q) from q
                        
                        if debugNewJavaMerging>1
                            checkSum=sum(Pointers);
                            fprintf(['loop #%d.%d, newPeak=%d, '...
                                'A=%d, B=%d, check sum=%d\n'], ...
                                outerLoop, innerLoop, newPeak,length(A), ...
                                length(B), checkSum);
                        end
                        
                        if any(tearAble ~= newPeak)
                            possibleClusterTears=[possibleClusterTears tearAble];
                        end
                    end
                    if debugNewJavaMerging==1
                        checkSum=sum(Pointers);
                        fprintf('Loop #%d.%d, check sum=%d\n', ...
                            outerLoop, innerLoop, checkSum);                        
                    end                    
                    changes=sum(Pointers~=newPointers);
                end
                possibleClusterTears=unique(possibleClusterTears);
                if ~noJava
                    try
                        if isempty(javaDbm.possibleClusterTears)
                        elseif any(possibleClusterTears~=javaDbm.possibleClusterTears')
                            msgBox('New JAVA merging POINTER problem');
                        end
                        if ~isequal(Pointers, javaDbm.pointers')
                            msgBox('JAVA merging POINTER problem');
                        end
                    catch
                    end
                    %% put in final clusters
                    javaDbm.fixClusterTear;
                end
                if ~isempty(possibleClusterTears)
                    possibleClusterTears=unique(possibleClusterTears);
                    if Density.IsDebugging
                        fprintf('*Possible* cluster tears ---> %s\n', ...
                            MatBasics.toString(possibleClusterTears));
                    end
                    [Pointers, actualClusterTears]=fixClusterTear(...
                        this.M, Pointers,  neighborhood, possibleClusterTears);
                    if isempty(actualClusterTears)
                        disp('No actual cluster tears');
                    else
                        nTears=sum(actualClusterTears(:,2));
                        fprintf('*ACTUAL* cluster tears --> %d\n',nTears);
                    end
                    if ~noJava
                        if any(javaDbm.tears ~= actualClusterTears)
                            disp('New java merging cluster tear problem');
                        end
                        if any(Pointers~=javaDbm.pointers')
                            disp('New JAVA merging POINTER problem after knitting torn clusters');
                        end
                    end
                end
                this.pointers=Pointers; %this is to save Pointers for making vector plot later
                this.rawPointers=Pointers;
                
                fprintf('MatLab merging completed\n\n');
                toc(matLabMerging);
            else
                javaDbm.fixClusterTear;
                Pointers=javaDbm.pointers';
                this.pointers=Pointers; %this is to save Pointers for making vector plot later
                this.rawPointers=Pointers;
            end
            p=this.pointers>0;
            while any(p)
                this.pointers(p)=this.pointers(this.pointers(p));  %follows the path of all positive pointers until they all go to dummy states
                p=this.pointers>0;
            end
            %PUBREF = STEP 6
            this.pointers(this.pointers==-1)=0;  %send background gridpoints to label 0
            
            eventClusterIds=this.pointers(this.eventBinIdxs); %assigns each original data point the dummy cluster number of its closest gridpoint
            
            clusts=unique(eventClusterIds);
            clusts=clusts(clusts~=0);  %all non-background clusters            
            numClusts=length(clusts);
            j=1;
            for i=clusts
                eventClusterIds(eventClusterIds==i)=j;  %relabel clusters from 1 to NumClusts
                this.pointers(this.pointers==i)=j;  %relabel pointers
                j=j+1;
            end
            if this.truncated
                ca=zeros(1, length(this.onScale));
                ca(this.onScale)=eventClusterIds;
                eventClusterIds=ca;
            end
            clustersWithoutEventBinIdxs=this.pointers<-1;
            if sum(clustersWithoutEventBinIdxs)>0
                if Density.IsDebugging
                    fprintf('There are %d grid points with clusters but no events!', sum(clustersWithoutEventBinIdxs));
                end
                this.pointers(clustersWithoutEventBinIdxs)=0;  %send background gridpoints that had clusters but no events
            end
            this.computePointersWithEvents(eventClusterIds);
            if options.mergeSmallNeighboringClusters && options.DbmMsncPerc>0
                [numClusts, eventClusterIds]=...
                    mergeSmallNeighboringClusters(this, numClusts, stdErr, ...
                    eventClusterIds, ...
                    options.DbmMsncPerc,...
                    options.DbmMsncDist);
            end
            if ~isempty(pu)
                pu.close;
            end
        end
        
        function setPointersWithEvents(this, pointersWithEvents)
            this.pointersWithEvents=pointersWithEvents;
        end
        
        function p=computePointersWithEvents(this, eventClusterIds)
            if isempty(eventClusterIds)
                this.pointersWithEvents=[];
            else
                if this.truncated
                    eventClusterIds=eventClusterIds(this.onScale);
                end
                this.pointersWithEvents=zeros(1, this.M^2);
                this.pointersWithEvents(this.eventBinIdxs)=eventClusterIds;
            end
            if nargout>0
                p=this.pointers;
            end
        end
        
        function [bi,xi,yi]=toBinIdxs(this, y,x)
            xi=floor((x-this.mins(1))./this.deltas(1))+1;
            yi=floor((y-this.mins(2))./this.deltas(2));
            bi=yi*this.M+xi;
            yi=yi+1;
        end
        
        function [result, resultStr]=computeDensity(this, gridX, gridY)
            phi = @(x) 1/sqrt(2*pi)*exp(-x.^2./2);
            xN=this.Z(1)*2+1;
            yN=this.Z(2)*2+1;
            xZ=-this.Z(1):this.Z(1);
            yZ=-this.Z(2):this.Z(2);
            phiXs=phi(xZ*this.deltas(1)/this.h(1))/this.h(1);
            phiYs=phi(yZ*this.deltas(2)/this.h(2))/this.h(2);
            result=0;
            for i=1:xN
                x=gridX+xZ(i);
                phiX=phiXs(i);
                if x>0 && x<=this.M
                    for j=1:yN
                        y=gridY+yZ(j);
                        if y>0 && y<=this.M
                            Wm=this.wmat(x,y);
                            num=Wm*(phiX*phiYs(j));
                            result=result+num;
                        end
                    end
                end
            end
            result=1/this.N*result;
            if nargout>1
                resultStr=String.encodeNumber(result,4);
            end
        end
        
        function decimals=getFmatDecimals(this)
            if isempty(this.decimals)
                this.decimals=abs(floor(log10(median(this.fmat(:)))));
                if isinf(this.decimals)
                    this.decimals=0;
                end
            end
            decimals=this.decimals;
        end
        
        function setOptions(this, options)
            this.opts=options;
        end
        
        function yes=hasSameOptions(this, options)
             yes=false;
             if isempty(this.opts) || isempty(options)
                 return;
             end
             if options.mergeSmallNeighboringClusters==this.opts.mergeSmallNeighboringClusters
                if options.DbmMsncPerc==this.opts.DbmMsncPerc
                    if options.DbmMsncDist==this.opts.DbmMsncDist
                        if options.DbmBandwidth==this.opts.DbmBandwidth
                            if options.DbmM==this.opts.DbmM
                                if options.DbmBackgroundFactor==this.opts.DbmBackgroundFactor
                                    if options.DbmSlopeIsSignificant==this.opts.DbmSlopeIsSignificant
                                        if options.DbmBackgroundType==this.opts.DbmBackgroundType
                                            if options.DbmNmin==this.opts.DbmNmin
                                                yes=true;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
             end
        end
    end
    
    
    methods(Access=private)
        function this=Density(events, options, mins_, maxs_)
            this.mins=mins_;
            this.maxs=maxs_;
            this.truncated=false;
            this.opts=options;
            M_=options.DbmM;
            if isfield(options, Density.DbmBandwidth)
                this.bandWidthInput=options.DbmBandwidth;
                this.bandwidth=this.bandWidthInput/100;
            else
                this.bandwidth=0;
            end
            this.onScale=MatBasics.FindOnScale(events, mins_, maxs_);
            cntEdge=size(events, 1)-sum(this.onScale);
            if cntEdge>0
                this.truncated=true;
                events=events(this.onScale, :);
            end
            this.M=M_;
            N_=size(events,1);  %number of data points
            this.N=N_;
            d=size(events,2);  %number of dimensions
            assert(d==Density.D, 'DBM is currently designed for 2 dimensions');
            MM = M_^d; %number of total gridpoints
            deltas_ = 1/(M_-1)*(maxs_-mins_);  %vector of distances between neighbor grid points in each dimension
            this.deltas=deltas_;
            % ym=gridpoints(d,1:M); %list of every possible combo of 1:M
            
            this.ye=zeros(d,M_);
            pointLL=zeros(N_,d);  %this will be the "lower left" gridpoint to each data point
            for i = 1:d
                this.ye(i,:) = linspace(mins_(i),maxs_(i),M_);
                pointLL(:,i)=floor((events(:,i)-mins_(i))./deltas_(i)) + 1;
            end
            pointLL(pointLL==M_)=M_-1;  %this avoids going over grid boundary
            %% assign each data point to its closest grid point
            [this.xgrid, this.ygrid]=meshgrid(this.ye(1,:),this.ye(2,:));
            z=reshape(1:MM,M_,M_);
            this.eventBinIdxs=interp2(this.xgrid, this.ygrid,z',...
                events(:,1),events(:,2),'nearest');  %this associates each data point with its nearest grid point
            %% compute w
            Deltmat=repmat(deltas_,N_,1);
            shape=M_*ones(1,d);
            
            wmat_=zeros(M_,M_);
            for i=0:1  %number of neighboring gridpoints in d dimensions
                for j=0:1
                    pointm=pointLL+repmat([j i],N_,1);  %indices of ith neighboring gridpoints
                    pointy=zeros(N_,d);
                    for k=1:d
                        pointy(:,k)=this.ye(k,pointm(:,k));  %y-values of ith neighboring gridpoints
                    end
                    %PUBREF=REF 2.1
                    W=prod(1-(abs(events-pointy)./Deltmat),2);  %contribution to w from ith neighboring gridpoint from each datapoint
                    wmat_=wmat_+accumarray(pointm,W,shape);  %sums contributions for ith gridpoint over data points and adds to wmat
                end
            end
            this.wmat=wmat_;
            
            %% compute f, sig, df and A
            this.h=zeros(1,d);
            this.contourH=zeros(1,d);
            this.bandwidth=DensityBiased.ConfirmBandWidth(this);
            if this.bandwidth>0
                if isfield(options, Density.DbmNmin)
                    this.N_min=options.DbmNmin;
                end
                for i =1:d
                    this.h(i)=this.bandwidth*(maxs_(i)-mins_(i));
                    this.contourH(i)=Density.CONTOUR_BANDWIDTH *(maxs_(i)-mins_(i));
                    if N_<this.N_min
                        if this.N_min<=9000
                            this.h(i)=(N_/this.N_min)^(-1/6)*this.h(i);
                        else
                            this.h(i)=this.h(i)/sqrt(N_/this.N_min);
                        end
                    end
                end
            else
                n6 = N_^(-1/6);
                for i =1:d
                    this.contourH(i)=Density.CONTOUR_BANDWIDTH*(maxs_(i)-mins_(i));
                    this.h(i) = std(events(:,i))*n6;
                end
            end
            Z=zeros(1,d);
            Zin=cell(1,d);
            for i =1:d
                Z(i) = min(floor(4*this.h(i)/deltas_(i)), M_-1);
                Zin{i}=-Z(i):Z(i);
            end
            this.Z=Z;
            %fprintf('h1=%f, h2=%f and Z1=%f Z2-%f\n', this.h(1), this.h(2), Z(1), Z(2))
            [bbwMat, nearN, nnVec, eigVec]=DensityBiased.Create(this, events);
            if ~isempty(bbwMat)
                this.nearN=nearN;
                this.isBiased=true;
                this.nnVec=nnVec;
                this.nnMat=reshape(nnVec, [M_ M_]);
                this.eigVec=eigVec;
                this.eigMat=reshape(eigVec, [M_ M_]);
                this.fmat=bbwMat;
                this.fmatVector=reshape(this.fmat,[1, MM]);
                this.fmatVector2=this.fmatVector;
                this.contourZ=bbwMat';
                [this.contourXm, this.contourYm]=meshgrid(this.ye(1,:),this.ye(2,:));
                return;
            end
            
            % gaussian phi(x) produces the value representing the height of
            % the bell curve at point x.  So phi(0) is the peak of the bell
            % curve which is .3989.  x from -5 to 5 gets a reasonable
            % result.  Usually what you put it into a phi function is
            % the result of some calculation that gets the domain between 5
            % and 5.  In DBM Z is defined such that Z*delta/h is within -4
            % to 4
            
            
            phi = @(x) 1/sqrt(2*pi)*exp(-x.^2./2);
            
            [this.L{1},this.L{2}]=meshgrid(Zin{1},Zin{2});
            
            Phix=phi(this.L{1}*deltas_(1)./this.h(1))./this.h(1);
            Phiy=phi(this.L{2}*deltas_(2)./this.h(2))./this.h(2);
            
            this.phimat = (Phix.*Phiy)';   %matrix of Phi for inputting into convn
            %PUBREF=REF 2.2
            this.fmat = 1/N_*convn(wmat_,this.phimat,'same');  %d-dim matrix of estimated densities
            this.fmatVector=reshape(this.fmat,[1, MM]);
            this.wmatVector=reshape(this.wmat, [1, MM]); 
            
            %fprintf(['Density sum=%4.5f mean=%4.7f stdDev=%4.7f;'...
            %    ' Data cnt=%d mean %1.4f/%1.4f \n'], ...
            %    sum(this.fmatVector), mean(this.fmatVector), ...
            %    std(this.fmatVector), size(events,1), ...
            %    mean(events(:,1)), mean(events(:,2)));
            
        end
        
        function contourDensity(this)
            d=2;
            h_=zeros(1,d);
            Z_=zeros(1,d);
            Zin=cell(1,d);
            for i =1:d
                h_(i)=this.contourH(i);
                if this.N<this.N_min
                    h_(i)=(this.N/this.N_min)^(-1/6)*1.7*h_(i);
                end
                Z_(i)=min(floor(4*h_(i)/this.deltas(i)), this.M-1);
                Zin{i}=-Z_(i):Z_(i);
            end
            phi = @(x) 1/sqrt(2*pi)*exp(-x.^2./2);
            [L_{1},L_{2}]=meshgrid(Zin{1},Zin{2});
            Phix=phi(L_{1}*this.deltas(1)./h_(1))./h_(1);
            Phiy=phi(L_{2}*this.deltas(2)./h_(2))./h_(2);
            Phimat = (Phix.*Phiy)';   
            fMat = 1/this.N*conv2(this.wmat,Phimat,'same');
            this.fmatVector2=reshape(fMat,[1, this.M^2]);
            this.contourZ=fMat';
            [this.contourXm, this.contourYm]=meshgrid(this.ye(1,:),this.ye(2,:));
        end         
    end
    
    methods(Static)
        function options=VeryHighOptions
            options.hWait=0;
            options.mergeSmallNeighboringClusters=1;
            options.DbmMsncPerc=0;
            options.DbmMsncDist=2;
            %bandWidth, M, backgroundType, backgroundFactor, slopeSignificance
            %1.5, 256, 2, .1, false
            options.DbmBandwidth=1.5;
            options.DbmM=256;
            options.DbmBackgroundFactor=.1;
            options.DbmSlopeIsSignificant=false;
            options.DbmBackgroundType=2;
            options.DbmNmin=5000;            
        end
        
        function DebugClusterPointers(fg, scrolling, mode)
            if nargin<3
                mode=1;
            end
            fg.density.debug=true;
            fg.dbm;
            fg.density.debug=false;
            html=['<html>' Qc.Html(fg, scrolling, mode) '</html>'];
            %File.SaveTextFile('qc.html', html);
            %Html.BrowseFile('qc.html');
            Html.BrowseString(html);
        end
        function n=LARGE_EVENTS
            n=250000;
        end
        function dbmVersion=DbmVersion
            %dbmVersion='v4';
            %if switching to v4 to support cluster boundaries 
            %below the W of the variable (See FcsGater.lowEnd)
            %then you MUST also change the static JAVA variable
            %edu.stanford.facs.swing.GatingTreeProps.DbmVersion='v4'
            %and RE-COMPILE (sigh)
            dbmVersion='v3';
        end
        
        function density=New(xy, detail, mins, maxs)
            if nargin<4
                maxs=[];
                if nargin<3
                    mins=[];
                    if nargin<2
                        detail='medium';
                    end
                end
            end
            if strcmpi(detail, 'low')
                    options=Density.Options(3.2, 256, 1, 4, false);
                elseif strcmpi(detail, 'most high')
                    options=Density.Options(0.95, 256, 2, .1, false);
                elseif strcmpi(detail, 'dbscan arguments')
                    options=Density.Options(1.5, 256, 2, .1, false);
                elseif strcmpi(detail, 'very high')
                    options=Density.Options(1.5, 256, 2, .1, false);
                elseif strcmpi(detail, 'high')
                    options=Density.Options(1.6, 256, 1, 4.3, false);
                elseif strcmpi(detail, 'medium')
                    options=Density.Options(2.3, 256, 1, 4.3, false);
                elseif strcmpi(detail, 'very low')
                    options=Density.Options(4, 256, 2, 1, true);
                elseif strcmpi(detail, 'adaptive')
                    options=Density.Options(0, 256, 4, 4.3, false);                    
                elseif strcmpi(detail, 'nearest neighbor')
                    options=Density.Options(-1, 128, 4, 4.3, false);
                else
                    warning([detail ' is not a DBM method, using medium']);
                    options=Density.Options(2.3, 256, 1, 4.3, false);
                end

            if isempty(mins) || isempty(maxs)
                density=Density.Create(xy,options);
            else
                density=Density.CreateWithScale(xy, options,mins,maxs);
            end
        end
        
        function options=Options(bandWidth, M, backgroundType,...
                backgroundFactor, slopeSignificance)
            if nargin<5
                slopeSignificance=false;
                if nargin<4
                    backgroundFactor=4.3;
                    if nargin<3
                        backgroundType=1;
                        if nargin<2
                            M=256;
                            if nargin<1
                                bandWidth=2.3;
                            end
                        end
                    end
                end
            end
            options.hWait=0;
            options.mergeSmallNeighboringClusters=1;
            options.DbmMsncPerc= 0;
            options.DbmMsncDist=2;
            options.DbmBandwidth=bandWidth;
            options.DbmM=M;
            options.DbmBackgroundFactor=backgroundFactor;
            options.DbmSlopeIsSignificant=slopeSignificance;
            options.DbmBackgroundType=backgroundType;
            options.DbmNmin=5000;
        end
        
        function this=Create(data, options)
            mins_=min(data);  %vector of min of each column of data
            maxs_=max(data);   %vector of max of each column of data
            this=Density(data, options, mins_, maxs_);
        end
        
        function this=CreateWithScale(data, options, mins, maxs)
            this=Density(data, options, mins, maxs);
        end
        
        function m=ToData(x, y, w)
            assert(length(x)==length(y) && length(x)==size(w,1) ...
                && length(x)==size(w,2));
            M_=length(x);
            M1=M_-1;
            m=zeros(M_^2,3);
            ww=w;
            xx=x';
            yy=y';
            for i=0:M_-1
                j=(i*M_)+1;
                ii=i+1;
                m(j:j+M1,:)=[xx, zeros(M_,1)+y(ii), ww(:,ii)];
            end
        end
        %It seems x and y are flipped from M matrix to MM vector but not so
        % Essentially in a MxM matrix the y index value is the same with each 
        % cell in the SAME row but it differs with each cell in the SAME
        % column.  Since y is usually height/column and x is usually
        % width/row this arrangement makes sense.
        function [x,y]=ToMatrixIdxs(vectorIdx, M)
            y=floor((vectorIdx-1)/M)+1;
            x=mod(vectorIdx-1,M)+1;
        end
        
        function vectorIdx=ToVectorIdx(x, y, M)
            y_=(y-1)*M;
            vectorIdx=y_+x;            
        end
        
        function ok=IsDebugging
            ok=false;
        end
        
        function ok=IsDebugging2014
            ok=false;
        end
        
        function Debug(M)
            if nargin<1
                M=32;
            end
            xSpace=linspace(8000,9000, 32);
            ySpace=linspace(2000,3000, 32);
            wMatrix=reshape(1:32^2, 32,32);
            wVector=reshape(wMatrix, 1, M^2);
            d=Density.ToData(xSpace,ySpace,wMatrix);
            matrixIndices=[M,1;M,3;3,M;1,M;M,2;2,M];
            for i=1:size(matrixIndices,1)
                x=matrixIndices(i,1);
                y=matrixIndices(i,2);
                v=Density.ToVectorIdx(x,y, M);
                [x_, y_]=Density.ToMatrixIdxs(v, M);
                assert(x_==x && y_==y);
                assert(wMatrix(x,y)==d(v,3) && wVector(v)==wMatrix(x,y));                
            end
            
            %test rectangular bin test for pixel resolution of data
            MX=10;
            MY=100;
            x=[2:9];
            y=[26:33];            
            pxMatrix=randi(500, MX, MY);
            pxVector=reshape(pxMatrix, 1, MX*MY);
            zz=Density.ToVectorIdx(x,y, MX);
            [xx, yy]=Density.ToMatrixIdxs(zz, MX);
            assert( isempty(find(x~=xx, 1)) && isempty(find(yy~=y, 1)));
            for i=1:length(xx)
                value=pxMatrix(xx(i), yy(i));
                vIdx=Density.ToVectorIdx(x(i), y(i), MX);
                assert(value==pxVector(vIdx));
                
            end
        end
        
        function wmat_=Weight(events,  M_, mins_, maxs_)
            deltas_ = 1/(M_-1)*(maxs_-mins_);  
            N_=size(events,1);  %number of data points
            d=size(events,2);  %number of dimensions
            Deltmat=repmat(deltas_,N_,1);
            shape=M_*ones(1,d);
            wmat_=zeros(M_,M_);
            data4Idx=zeros(d,M_);
            lowIdx4Data=zeros(N_,d);  %this will be the "lower left" gridpoint to each data point
            for i = 1:d
                data4Idx(i,:)=linspace(mins_(i),maxs_(i),M_);
                lowIdx4Data(:,i)=floor((events(:,i)-mins_(i))./deltas_(i)) + 1;
            end
            lowIdx4Data(lowIdx4Data==M_)=M_-1;  %this avoids going over grid boundary

            for i=0:1  %number of neighboring gridpoints in d dimensions
                for j=0:1
                    curIdxForData=lowIdx4Data+repmat([j i],N_,1);  %indices of ith neighboring gridpoints
                    dataAtCurIdx=zeros(N_,d);
                    for k=1:d
                        dataAtCurIdx(:,k)=data4Idx(k,curIdxForData(:,k));  %y-values of ith neighboring gridpoints
                    end
                    %PUBREF=REF 2.1
                    W=prod(1-(abs(events-dataAtCurIdx)./Deltmat),2);  %contribution to w from ith neighboring gridpoint from each datapoint
                    wmat_=wmat_+accumarray(curIdxForData,W,shape);  %sums contributions for ith gridpoint over data points and adds to wmat
                end
            end
        end
        
    end 
end
