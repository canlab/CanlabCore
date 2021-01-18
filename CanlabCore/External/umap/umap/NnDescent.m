%
%   AUTHORSHIP
%   Primary Developer:  Stephen Meehan <swmeehan@stanford.edu> 
%   Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%                           Jonathan Ebrahimian <jebrahimian@mail.smu.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
classdef NnDescent<handle
    
    methods(Static)
        function TestNnDescent(X, Y, varargin)
            if nargin<3
                varargin={};
                if nargin<2
                    Y=nan;
                    if nargin<1
                        X='cytofExample1';
                    end
                end
            end
            varargin{end+1}='example_finder_callback';
            varargin{end+1}=@UmapUtil.RelocateExamples;
            KnnFind.Test(X, Y, varargin{:});
        end
        
        function TestSelfOther(other, self, dropEndLabelColumn, ...
                n_neighbors, metric, transform_queue_size, dist_args,percentRows)
            if nargin<8
                percentRows=1;
                if nargin<7
                    dist_args=[];
                    if nargin<6
                        transform_queue_size=UmapUtil.DFLT_TRANSFORM_Q_SZ;
                        if nargin<5
                            metric='euclidean';
                            if nargin <4
                                n_neighbors=15;
                                if nargin<3
                                    dropEndLabelColumn=true;
                                    if nargin<2 || isempty(self)
                                        if nargin==0
                                            other=0;
                                        elseif ischar(other)
                                            other=str2double(other);
                                        else
                                            other=other;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if isempty(metric)
                metric='euclidean';
            end
            if strcmp(metric, 'euclidean')
                sfx='';
            else
                sfx=['_' metric];
            end
            if isempty(n_neighbors)
                n_neighbors=15;
            end
            if isempty(dropEndLabelColumn)
                dropEndLabelColumn=true;
            end
            if isnumeric(other)
                code=other;
                if isnan(code)
                    code=0;
                end
                if code==29
                    self='s2_samusikImported_29D';
                    other='s1_samusikImported_29D';
                elseif code==27
                    self='s9_genentech_27D';
                    other='s1_genentech_27D';
                elseif code==10
                    other='sampleBalbcLabeled55k';
                    self='sampleRagLabeled148k';
                else
                    self='testSet';
                    dropEndLabelColumn=false;
                    other='trainingSet';
                end

            end
            matFile1=[other sfx '.mat'];
            if ~String.EndsWithI(other, '.csv') 
                other=[other '.csv'];
            end
            if ~String.EndsWithI(self, '.csv') 
                self=[self '.csv'];
            end
            resolvedFiles=UmapUtil.RelocateExamples({other, self});
            p=fileparts(resolvedFiles{1});
            matFile1=fullfile(p, matFile1);
            [otherData, labels]=read(other);
            [R1, C1]=size(otherData);
            otherSize=[String.encodeInteger(R1) ' X ' ...
                String.encodeInteger(C1) ];
            if ~exist(matFile1, 'file')
                fprintf(['FIRST doing MatLab knnsearch for EXACT search '...
                    'graph for training set (%s)\n'], otherSize);
                [knnIndices, knnDists]=KnnFind.Approximate(otherData, [],...
                    'K', n_neighbors, 'metric', metric, ...
                    'dist_args', dist_args);
                search_graph=NnDescent.SearchGraph(...
                    otherData, n_neighbors, knnIndices, knnDists, labels);
                save(matFile1, 'search_graph');
            else
                load(matFile1, 'search_graph');
            end

            selfData=read(self);
            [R2, C2]=size(selfData);
            selfSize=[String.encodeInteger(R2) ' X ' ...
                String.encodeInteger(C2) ];
            allSize=[selfSize ' in ' otherSize ', n_neighbors=' ....
                num2str(n_neighbors) ' "' metric '"'];
            
            timeCpp=tic;
            fprintf('Starting C++ knnsearch of %s\n', ...
                allSize);
            cb=[];
            %cb=@test_callback;
            
            knnIndices2=KnnFind.Approximate(otherData, selfData, ...
                'K', n_neighbors, 'metric', metric, ...
                'dist_args', dist_args, 'progress_callback', cb, ...
                'X_search_graph', search_graph, ...
                'nn_descent_transform_queue_size', transform_queue_size);
            fprintf('Done ... C++ took %s searching %s \n', ...
                String.MinutesSeconds(toc(timeCpp)), allSize);

            fprintf(['Starting MatLab knnsearch of %s\n'], allSize);
            timeMatLab=tic;
            knnIndices1=KnnFind.Determine(otherData, ...
                selfData, struct('K', n_neighbors, 'Distance', metric, ...
                'dist_args', dist_args));
            fprintf('Done ... MatLab took %s searching %s \n', ...
                String.MinutesSeconds(toc(timeMatLab)), allSize);

       
            p = KnnFind.AssessApproximation(knnIndices1, knnIndices2);
            disp(['        C++ accuracy is ' num2str(100*p) ' percent!']);
            disp('DONE C++ knnsearch of self in other');
            
            function [data, labels]=read(file)
                f=UmapUtil.RelocateExamples({file});
                file=f{1};
                data=File.ReadCsv(file);
                if dropEndLabelColumn
                    data=data(:,1:end-1);
                    labels=data(:,end);
                else
                    labels=[];
                end
                if percentRows <1 && percentRows>0
                    [R,C]=size(data);
                    rows=floor(percentRows*R);
                    data=data(1:rows,:);
                end
            end
            
            function ok=test_callback(msg)
                disp(msg);
                ok=true;
            end

        end

        function search_graph=SearchGraph(...
            data, n_neighbors, knnIndices, knnDists, labels, target_weight)
            if nargin<6
                target_weight=.5;
            end
            graph = fuzzy_simplicial_set(data, n_neighbors, 'euclidean',...
                    'knn_indices', knnIndices, 'knn_dists', knnDists, ...
                    'set_op_mix_ratio', 1,'local_connectivity', 1);
            if nargin>4 && ~isempty(labels)
                far_dist = 2.5 * (1 / (1 - target_weight));
                graph=categorical_simplicial_set_intersection(graph, labels, 1, far_dist);
            end
            search_graph = (graph > 0) + speye(size(graph,1));
        end
        
        function [indptr,indices,search_graph]=SearchGraphIndptrIndices(...
            data, n_neighbors, knnIndices, knnDists, labels, target_weight)
            if nargin<6
                target_weight=.5;
            end
            search_graph=NnDescent.SearchGraph(data, n_neighbors, ...
                knnIndices, knnDists, labels, target_weight);
            [indptr, indices]=KnnFind.IndPtrAndIndices(search_graph);
        end
        
        function MakeTestExample(other, self, dropEndLabelColumn, n_neighbors)
            if nargin <4
                n_neighbors=15;
                if nargin<3
                    dropEndLabelColumn=true;
                    if nargin<2
                        self='testSet';
                        dropEndLabelColumn=false;
                        if nargin<1
                            other='trainingSet'; 
                        end
                    end
                end
            end
            [otherData, other, labels]=read(other);
            %MUST create unlabeled file for C++ testing
            File.SaveMatrix([other '.csv' ], otherData, true); 
            selfTime=tic;
            [knnIndices, knnDists]=nearest_neighbors(otherData, n_neighbors, 'euclidean');
            disp('Done finding knn of self in self');
            toc(selfTime)
            File.SaveMatrix([other '.knnIndices.csv' ], int32(knnIndices-1), true);
            File.SaveMatrix([other '.knnDists.csv' ], knnDists, true);
            [indptr, indices]=NnDescent.SearchGraphIndptrIndices(...
                otherData, n_neighbors, knnIndices, knnDists, labels);
            File.SaveMatrix([other '.indptr.csv' ], int32(indptr), true);
            File.SaveMatrix([other '.indices.csv' ], int32(indices), true);
            
            [selfData, self]=read(self);
            %MUST create unlabeled file for C++ testing
            File.SaveMatrix([self '.csv' ], selfData, true); 
            selfOtherTime=tic;
            [knnIndices, knnDists] = knnsearch(otherData, selfData,'K', n_neighbors,'Distance', 'euclidean');
            disp('Done finding knn of self in other');
            toc(selfOtherTime);
            File.SaveMatrix([self '.knnIndices.csv' ], int32(knnIndices-1), true);
            File.SaveMatrix([self '.knnDists.csv' ], knnDists, true);
            
            function [data, file, labels]=read(file)
                if ~String.EndsWithI(file, '.csv') 
                    file=[file '.csv'];
                end
                f=UmapUtil.RelocateExamples({file});
                file=f{1};
                data=File.ReadCsv(file);
                if nargout>2
                    labels=data(:,end);
                end
                if dropEndLabelColumn
                    data=data(:,1:end-1);
                    labels=data(:,end);
                else
                    labels=[];
                end
                file=[file(1:end-4) '_ul'];
            end
        end

        function MakeBasicTestExample()
            files=UmapUtil.RelocateExamples({'trainingSet.csv', ...
                'trainingSet.labels.csv', 'testSet.csv'});
            trainingSetData=File.ReadCsv(files{1});
            trainingSetFile=files{1}(1:end-4);
            testSetData=File.ReadCsv( files{3});
            [knnIndices, knnDists]=nearest_neighbors(trainingSetData, 15, 'euclidean');
            File.SaveMatrix([trainingSetFile '.knnIndices.csv' ], knnIndices-1, true);
            File.SaveMatrix([trainingSetFile '.knnDists.csv' ], knnDists, true);
            labels=File.ReadCsv(files{2});
            [indptr2, indices2]=NnDescent.SearchGraphIndptrIndices(...
                trainingSetData, 15, knnIndices, knnDists, labels);
            
            testSetFile=files{3}(1:end-4);
            [knnIndices, knnDists] = knnsearch(trainingSetData, testSetData,'K', 15,'Distance', 'euclidean');
            File.SaveMatrix([testSetFile '.knnIndices.csv' ], knnIndices-1, true);
            File.SaveMatrix([testSetFile '.knnDists.csv' ], knnDists, true);
            trainingSetData=[trainingSetData labels];
            [M,U]=run_umap(trainingSetData, 'label_column', 'end', ...
                'label_file', 'balbcLabels.properties');
            [indptr, indices]=KnnFind.IndPtrAndIndices(U.search_graph);
            File.SaveMatrix([trainingSetFile '.indptr.csv' ], indptr, true);
            File.SaveMatrix([trainingSetFile '.indices.csv' ], indices, true);
            disp('CSV files to run C++ nn descent on training and test set');
            disp(['    in folder ' fileparts(files{1}) ' are trainingSet.csv, ']);
            disp( '    testSet.csv, trainingSet.indptr.csv AND trainingSet.indices.csv');
            disp('  ');
            disp('To run same thing in python first CREATE the UST template ')
            disp(['   doUmap.py ' trainingSetFile ...
                '.csv --n_neighbors=15 --metric=euclidean --min_dist 0.3 --firstRow 0 --verbose --saveTemplate ' ...
                trainingSetFile '.template --output_dimensions 2 --labels ' trainingSetFile '.labels.csv'])
            output2=files{3}(1:end-4);
            disp('..and then APPLY the template with  ')
            disp(['   doUmap.py ' output2 '.csv --n_neighbors 15 --metric=euclidean  --min_dist 0.3 --firstRow 0 --verbose  --output_dimensions 2 --useTemplate ' trainingSetFile '.template.umap'])
        end
        

        function Test(varargin)
            if isempty(varargin)
                disp('Usage: NnDescent.Test[vector of column sizes], [list of row sizes], [list of neighbors]');
                disp('   will default to ');
                disp('      Test_NN [10 29],  [1e3 1e4 5e4], [15 30]');
                disp(' here we go ...');
                varargin={ [10, 29]};
            end
            if length(varargin)<2
                varargin{2}=[1e3 1e4 5e4];
            end
            if length(varargin)<3
                varargin{3}=[15 30];
            end
            UmapUtil.Initialize();
            eliverData=[];
            samusikData=[];
            
            Ds=varargin{1};
            for D=1:length(Ds)
                doTest(Ds(D), varargin{2}, varargin{3});
            end
            
            
            function doTest(C, rowTests, neighborTests)
                if C<=10
                    if isempty(eliverData)
                        disp('Reading data from Eliver Ghosn');
                        eliverData=File.ReadCsv(...
                            UmapUtil.RelocateExamples('sample130k.csv'));
                    end
                    data=eliverData(:, 1:C);
                else
                    if isempty(samusikData)
                        disp('Reading data from Nikolay Samusik');
                        samusikData=File.ReadCsv(...
                            UmapUtil.RelocateExamples(...
                                's1_samusikImported_29D_ul.csv'));
                    end
                    data=samusikData(:, 1:C);
                end
                if isempty(neighborTests)
                    neighborTests=[15 30];
                    if nargin<2
                        rowTests=[1e3 1e4 5e4];
                    end
                end
                R=size(data, 1);
                rowTests(end+1)=R;
                dist_args=1.35;
                metric='minkowski';
                for i = rowTests
                    for j = neighborTests
                        scope=sprintf('data=%d-by-%d, neighbors=%d ...',...
                            i,C,j);
                        fprintf(['Running knnsearch (via '...
                            'nearest_neighbors.m) on\n\t%s'], scope);
                        tic
                        [indices1, dists1] = nearest_neighbors(...
                            data(1:i,:), j, metric, dist_args);
                        t=toc;
                        fprintf(' took %s seconds !!\n', String.encodeRounded(t, 2));
                        tic
                        fprintf('    >>> now C++ nn_descent on \n\t%s', scope);
                        indices2=KnnFind.Approximate(data(1:i,:), ...
                            [],'K', n_neighbors, 'metric', metric, ...
                            'dist_args', dist_args, 'progress_callback',...
                            @test_callback);
                        t=toc;
                        fprintf(' took %s seconds ...\n', String.encodeRounded(t, 2));
                        p = KnnFind.AssessApproximation(indices1, indices2);
                        disp(['        C++ accuracy is ' num2str(100*p) ' percent!']);
                    end
                end
            end
            function ok=test_callback(msg)
                disp(msg);
                ok=true;
            end
        end
    end
end
