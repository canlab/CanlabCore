classdef UstTest
%
%   AUTHORSHIP
%   Primary Developer:  Stephen Meehan <swmeehan@stanford.edu>
%   Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause

    properties(Constant)
        DFLT_SAMPLE_SET='samusikImported_29D';
    end
    
    methods(Static)
        function [args, pu]=Initialize(typeOfTesting, ...
                dfltMatchTypes, dfltMatchScenarios, varargin)
            [args, argued]=UmapUtil.Initialize(varargin{:});
            if isempty(args)
                pu=[];
                return;
            end
            pu=PopUp(['Determining scope of "sample ' ...
                typeOfTesting '" testing...'], ...      
                'north west+', 'UstTest: starting...', false, false, ...
                'genieSearch.png');
            args.parent_popUp=pu;
            args=UmapUtil.NewArgDefault(args, argued, ...
                'match_table_fig', false);
            args=UmapUtil.NewArgDefault(args, argued, ...
                'match_histogram_fig', false);
            args=UmapUtil.NewArgDefault(args, argued, ...
                'match_supervisors', dfltMatchTypes);
            args=UmapUtil.NewArgDefault(args, argued, ...
                'match_scenarios',  dfltMatchScenarios);
            args=UmapUtil.NewArgDefault(args, argued, ...
                'verbose', 'text');
            args=UmapUtil.NewArgDefault(args, argued, ...
                'see_training', true);
            args=UmapUtil.NewArgDefault(args, argued, ...
                'cluster_output', 'ignore');
            args=UmapUtil.NewArgDefault(args, argued, ...
                'label_column', 'end');
            args=UmapUtil.NewArgDefault(args, argued, ...
                'ust_test_components', [3 2]);
            args=UmapUtil.NewArgDefault(args, argued, ...
                'n_neighbors', 15);
            args=UmapUtil.NewArgDefault(args, argued, ...
                'ust_test_basic_reduction', true);
            args=UmapUtil.NewArgDefault(args, argued, ...
                'color_file', 'colorsByName.properties');
        end
        
        function [ustFldr, testId, runs]=GetResultsFldr(...
                typeOfTesting, fldr, fn)
            ustFldr=fullfile(fldr, ['ustTest.' typeOfTesting], fn);            
            File.mkDir(ustFldr);
            runs=Map;
            runs.load(fullfile(ustFldr, 'runs.mat'));
            testId=runs.newId;
        end
        
        function Split(varargin)
            [args, pu]=UstTest.Initialize('split', [1 2 3 4],...
                [ 1 2 3 4], varargin{:});
            extras=[];
            tick=tic;
            if isempty(args)
                warning('No arguments provided');
                return;
            end
            csvSampleFile=args.csv_file_or_data;
            assert( ~isempty( csvSampleFile) && ischar(csvSampleFile),...
                'First argument must be csv file of data');
            [files, existence]=UmapUtil.RelocateExamples({csvSampleFile,...
                File.SwitchExtension(csvSampleFile, '.properties')});
            if ~all(existence)
                pu.close;
                return;
            end
            csvSampleFile=files{1};
            labelFile=files{2};
            [fldr, fn]=fileparts(csvSampleFile);
            [ustFldr, testId]=UstTest.GetResultsFldr('split', fldr, fn);
            try
                [dataSet, cn]=File.ReadCsv(csvSampleFile);
            catch ex
                msg(Html.Wrap(['Could not open your csv file<br>&nbsp;&nbsp;&nbsp;'...
                    BasicMap.Global.smallStart csvSampleFile ...
                    '<br><br><b>BECAUSE</b>:<br>&nbsp;&nbsp;&nbsp;&nbsp;' ex.message ...
                    BasicMap.Global.smallEnd '<hr>']), 0, ...
                    'center', 'Error...', 'error.png');
                pu.close
                return;
            end
            if args.ust_test_synthesize>0
                [syntheticData, syntheticLabels]=generate_synthetic_set(...
                    dataSet(:,1:end-1), dataSet(:,end), ...
                    args.ust_test_synthesize*size(dataSet,1));
                dataSet=[syntheticData syntheticLabels];
                strSyn=['_syn' ...
                        num2str(floor(args.ust_test_synthesize*100))];
            else
                strSyn='';
            end
            args.result_folder=ustFldr;
            args.match_file=['match_' num2str(testId) ];
            falsePosNegFile=[UMAP_extra_results.FileName(args, 'falsePosNeg') '.txt'];
            falsePosNegFile=File.AppendSuffix(falsePosNegFile, strSyn);
            if exist(falsePosNegFile, 'file') && ...
                    (ismember(4, args.match_scenarios) ||...
                    ismember(3, args.match_scenarios))
                [answer, cancelled]=Gui.Ask(...
                    'Prior false +/- results exist...',...
                    {'Redo UST classification ...', ...
                    'View +/- of prior classification'}, ...
                    'prevFalsePosNeg', 'Redo?', 2);
                if cancelled
                    pu.close;
                    return;
                elseif answer==2
                    pu.close;
                    FalsePositiveNegative.Plot([0 1], falsePosNegFile)
                    return;
                end
            end
            close all force;
            try
                Gui.DisposeAllJavaWindows;
            catch
            end
            
            csvFile=[UMAP_extra_results.FileName(args) '.csv'];
            csvFile=File.AppendSuffix(csvFile, strSyn);
            args.parameter_names=cn(1:end-1);
            args.see_training=true;
            [R, hiD]=size(dataSet);
            idxsFile=fullfile(ustFldr, ['sampleIdxs' strSyn '.mat']);
            if exist(idxsFile, 'file')
                load(idxsFile, 'sampleIdxs');
            else
                sampleIdxs=1:R;
                sampleIdxs=sampleIdxs(randperm(length(sampleIdxs)));
                save(idxsFile, 'sampleIdxs');
            end
            third=floor(R/3);
            testSetIdxs=sampleIdxs(third+1:end);
            testSetData=dataSet(testSetIdxs, :);
            dataBackUp=testSetData;
            htmlBody='';
            htmlHead1='';
            htmlHead2='';
            htmlHead3='';
            csvBody='';
            csvHead='';
            tabBody='';
            tabHead='';
            
            good=true;
            firstCase=true;
            nTestCases=length(args.ust_test_cases);
            curTest=1;
            if args.ust_test_both_thirds
                nSplits=2;
            else
                nSplits=1;
            end
            [nLoDs, nRunUmaps, dataSetsTxt, dataSetTxt]=...
                UmapUtil.InitProgress(args, pu, nTestCases, nSplits, R, hiD-1);
            for loD=1:nLoDs
                args.n_components=args.ust_test_components(loD);
                reductionSummary=[dataSetTxt '-' ...
                    num2str(args.n_components) 'D'];
                pu.setText(['Reducing ' reductionSummary]);
                if args.ust_test_synthesize>0
                    tmpltFile=fullfile(ustFldr, ['template' strSyn '_' ...
                        num2str(hiD) 'D_' num2str(args.n_neighbors) 'nn_'...
                        num2str(args.n_components) 'D.mat']);
                else
                    tmpltFile=fullfile(ustFldr, ['template_' ...
                        num2str(hiD) 'D_' num2str(args.n_neighbors) 'nn_'...
                        num2str(args.n_components) 'D.mat']);
                end
                args.template_file=tmpltFile;
                if ~exist(tmpltFile, 'file')
                    trainingSetIdxs=sampleIdxs(1:third);
                    reduced=run_umap(dataSet(trainingSetIdxs,:),...
                        'label_column', 'end',...
                        'parameter_names', cn, ...
                        'save_template_file', tmpltFile, ...
                        'n_components', args.n_components, ...
                        'color_file', args.color_file,...
                        'n_neighbors', args.n_neighbors,...
                        'label_file', labelFile);
                    if isempty(reduced)
                        pu.close;
                        return;
                    end
                end
                for i=1:nTestCases
                    testCase=args.ust_test_cases(i);
                    if testCase==2
                        isRand=true;
                        perturbation=0;
                        freqMean=args.ust_test_freq_mean;
                    elseif testCase==3
                        isRand=true;
                        perturbation=args.ust_test_perturbation;
                        freqMean=1;
                    elseif testCase==4
                        isRand=true;
                        perturbation=args.ust_test_perturbation;
                        freqMean=args.ust_test_freq_mean;
                    elseif testCase==12
                        isRand=true;
                        perturbation=0;
                        freqMean=args.ust_test_freq_mean*2;
                    elseif testCase==13
                        isRand=true;
                        perturbation=args.ust_test_perturbation*2;
                        freqMean=1;
                    elseif testCase==14
                        isRand=true;
                        perturbation=args.ust_test_perturbation*2;
                        freqMean=args.ust_test_freq_mean*2;
                    else
                        isRand=false;
                        perturbation=0;
                        freqMean=1;
                    end
                    if perturbation>0
                        testSetData(:, 1:end-1)=perturb_classes(...
                            testSetData(:, 1:end-1), testSetData(:,end));
                    end
                    if freqMean<1
                        elimIdxs=ln_freq_var(testSetData(:,end), freqMean);
                        testSetData(elimIdxs,:)=[];
                    end
                    
                    if ~isRand
                        subset1=testSetData(1:third,:);
                        subset2=testSetData(third+1:end,:);
                    else
                        R2=size(testSetData,1);
                        third=floor(R2/3);
                        a=1:third*2;
                        a=a(randperm(length(a)));
                        split=floor(length(a)/2);
                        idxs1=a(1:split);
                        idxs2=a(split+1:end);
                        subset1=testSetData(idxs1,:);
                        subset2=testSetData(idxs2,:);
                    end
                    
                    if runTest(subset1, isRand, 1)
                        if args.ust_test_both_thirds
                            if ~runTest(subset2, isRand, 2)
                                good=false;
                                break;
                            end
                        end
                    else
                        good=false;
                        break;
                    end
                    testSetData=dataBackUp;
                end
            end
            if ~good
                return;
            end
            tick2=String.MinutesSeconds(toc(tick));
            if ~isempty(csvFile)
                File.SaveTextFile(csvFile, [csvHead csvBody])
            end
            if ~isempty(tabBody) && ~isempty(falsePosNegFile)
                File.SaveTextFile(falsePosNegFile, [tabHead tabBody]);
                FalsePositiveNegative.Plot([0 1], falsePosNegFile)
            end
                
            UmapUtil.SeeTableHtml(htmlHead1, htmlHead2, htmlHead3, ...
                htmlBody, 1, extras.htmlFile, ...
                ['<h1>UST results: ' dataSetsTxt '</h1>'], ...
                ['<h2>Runtime: ' tick2 '</h2>'], ...
                ['<h3>' csvSampleFile '</h3>'] );
            pu.close;
            
            
            function ok=runTest(subset, isRand, part)                
                ok=true;
                if part==1
                    args.cascade_x=0;
                else
                    args.cascade_x=470;
                end
                if isRand
                    sRand=';rand';
                else
                    sRand='';
                end
                args.sample_set=['split ' testId];
                args.test_set=sprintf(...
                    '%d-%d: %s;%s', testCase,...
                    part, String.encodePercent(perturbation, 1, 1), ...
                    String.encodePercent(freqMean, 1, 1));
                args.description=[args.test_set sRand, reductionSummary];
                args.context.testCase=testCase;
                args.context.testPart=part;
                args.context.isRand=isRand;
                args.context.perturbation=perturbation;
                args.context.freqMean=freqMean;
                args.csv_file_or_data=subset;
                
                args.parent_context=sprintf('%d/%d', curTest, ...
                    nLoDs * nRunUmaps);
                [~,~,~,extras]=run_umap(args);
                curTest=curTest+1;
                if ~isempty(extras.matchHtmlHead1)
                    if firstCase
                        htmlHead1=extras.matchHtmlHead1;
                        htmlHead2=extras.matchHtmlHead2;
                        htmlHead3=extras.matchHtmlHead3;
                    end
                    html=[extras.matchHtmlBody];
                    csvHead=extras.matchCsvHead;
                    csvBody=[csvBody extras.matchCsvBody];
                    tabHead=extras.falsePosNegHead;
                    tabBody=[tabBody extras.falsePosNegBody];
                else
                    ok=false;
                    return;
                end
                if args.ust_test_basic_reduction
                    tmplt=args.template_file;
                    args.template_file='';
                    msprv=args.match_supervisors;
                    args.match_supervisors=1;
                    msc=args.match_scenarios;
                    args.match_scenarios=[3 4];
                    args.description='';
                    args.label_file=labelFile;
                    args.parent_context=sprintf('%d/%d', curTest, ...
                        nLoDs*nRunUmaps);
                    [~,~,~,extras]=run_umap(args);
                    curTest=curTest+1;
                    args.template_file=tmplt;
                    args.match_supervisors=msprv;
                    args.match_scenarios=msc;
                    if ~isempty(extras.matchHtmlHead1)
                        if firstCase
                            htmlHead1=[htmlHead1 extras.matchHtmlHead1];
                            htmlHead2=[htmlHead2 extras.matchHtmlHead2];
                            htmlHead3=[htmlHead3 extras.matchHtmlHead3];
                        end
                        html=[html extras.matchHtmlBody];
                        csvHead=extras.matchCsvHead;
                        csvBody=[csvBody extras.matchCsvBody];
                        tabHead=extras.falsePosNegHead;
                        tabBody=[tabBody extras.falsePosNegBody];
                    else
                        ok=false;
                    end
                end
                htmlBody=[htmlBody '<tr>' html '</tr>'];
                firstCase=false;
            end            
        end
        function fileName=Merge(sampleClues, synthesize, sampleSetClue)
            if nargin<3
                sampleSetClue=29;
                if nargin<2
                    synthesize=0;
                    if nargin<1
                        sampleClues=[1 2];
                    end
                end
            end
            fileName=UstTest.MergeSamples(sampleSetClue, true,...
                synthesize, [], sampleClues);
        end

        function mergedFile=MergeSamples(sampleSetClue, ...
                askToRegenerate, synthesize, pu, varargin)
            [~,samples, labels, ~, mergedFile, ~, strSyn]...
                =UstTest.SolveClues(varargin{:},varargin{:}, ...
                sampleSetClue, [], [], synthesize);
            N=length(samples);
            if N<2 || mod(N, 2) ~= 0
                disp('You must provide pairs of samples');
                return;
            end
            [files, existence]=UmapUtil.RelocateExamples(...
                [samples, labels]);
            if ~all(existence)
                return;
            end
            if synthesize>0
                mergedFile=File.AppendSuffix(mergedFile, strSyn);
            end
            fldr=fileparts(files{1});
            mergedFile=fullfile(fldr, mergedFile);
            if exist(mergedFile, 'file')
                if ~askToRegenerate || ...
                        ~askYesOrNo(Html.WrapHr(...
                        'Merger exists ... regenerate anyway?'), ....
                        'Ooops-y ....')
                    return;
                end
            end
            if ~isempty(pu)
                pu.setText2(['Merging ' String.Pluralize2('sample', N)]);
            end
            mergedData=[];
            firstLabelFile=files{N+1};
            for i=1:2:N
                [data, columns]=...
                    LabelBasics.Merge2Samples([], files{i}, ...
                    files{i+N}, files{i+1}, files{N+i+1}, firstLabelFile);
                mergedData=[mergedData; data];
            end
            if synthesize>0
                [syntheticData, syntheticLabels]=...
                    generate_synthetic_set(...
                    mergedData(:,1:end-1), mergedData(:,end), ...
                    synthesize*size(mergedData,1));
                mergedData=[syntheticData syntheticLabels];
                mergedFile=File.AppendSuffix(mergedFile, strSyn);
            end
            fldr=fileparts(files{1});
            t=array2table(mergedData);
            t.Properties.VariableNames=columns;
            writetable(t, mergedFile);
            copyfile(firstLabelFile, File.SwitchExtension(mergedFile, '.properties')) 
        end
        
        function sampleSet=SolveSampleSetClue(sampleSetClue)
            if isempty(sampleSetClue)
                sampleSet=UstTest.DFLT_SAMPLE_SET;
            elseif isnumeric(sampleSetClue)
                switch sampleSetClue
                    case 19
                        sampleSet='samusikImported_19D';
                    case 29
                        sampleSet=UstTest.DFLT_SAMPLE_SET;
                    case 29.1
                        sampleSet=UstTest.DFLT_SAMPLE_SET;
                    case 27
                        sampleSet='genentech';
                    case 29.2
                        sampleSet='samusikManual_29D';
                    otherwise
                        sampleSet=UstTest.DFLT_SAMPLE_SET;
                end
            else
                assert(ischar(sampleSetClue));
                sampleSet=sampleSetClue;
            end
        end
        
        function [tmplt, ss, tSmpl]=Supervise(sampleSetClue, ...
                templateClue, onlyCreate, loDs, nns)
            if nargin<5
                nns=15;
                if nargin4
                    loDs=3;
                    if nargin<3
                        onlyCreate=false;
                        if nargin<2
                            templateClue=1;
                            if nargin<1
                                sampleSetClue=29.2;
                            end
                        end
                    end
                end
            end
            ss=UstTest.SolveSampleSetClue(sampleSetClue);
            tSmpl=['s' num2str(templateClue) '_' ss ];
            nDs=length(loDs);
            nNns=length(nns);
            for h=1:nNns
                for i=1:nDs
                    D=loDs(i);
                    nn=nns(h);
                    tmplt=['ust_' tSmpl '_' num2str(nn) 'nn_' num2str(D) 'D.mat'];
                    tmplt=UmapUtil.RelocateExamples(tmplt);
                    if ~onlyCreate || ~exist(tmplt, 'file') % gotta build template first
                        r=run_umap([tSmpl '.csv'], ...
                            'label_file', [tSmpl '.properties'], ...
                            'label_column', 'end', ...
                            'save_template_file', tmplt ,...
                            'n_components', D, ...
                            'n_neighbors', nn,...
                            'color_file', 'colorsByName.properties');
                        if isempty(r)
                            return;
                        end
                    end
                end
            end
        end
        
        function [sampleSet, samples, labels, templates, templateSample, ...
                templateLabel, strSyn]=SolveClues(sampleClues, ...
                templateClue, sampleSetClue, loDs, nn, synthesize)
            if nargin<6
                synthesize=0;
                if nargin<5
                    nn=15;
                    if nargin<4
                        loDs=[3 2];
                        if nargin<3
                            sampleSetClue=[];
                            if nargin<2
                                templateClue=1;
                                if nargin<1
                                    sampleClues=2;
                                end
                            end
                        end
                    end
                end
            end
            if synthesize>0
                strSyn=['_syn' ...
                    num2str(floor(synthesize*100))];
            else
                strSyn='';
            end            
            sampleSet=UstTest.SolveSampleSetClue(sampleSetClue);
            nD=length(loDs);
            templates=cell(1,nD);
            N=length(sampleClues);
            labels=cell(1, N);
            if ~isnumeric(sampleClues)
                samples=sampleClues;
                assert(StringArray.IsValid(samples), ...
                    'sampleClues must be sample names');
                for i=1:N
                    labels{i}=File.SwitchExtension(samples{i}, '.properties');
                end
                if ischar(templateClue)
                    templateSample=templateClue;
                    [~,tSample]=fileparts(templateSample);
                    for i=1:nD
                        D=loDs(i);
                        templates{i}=['ust_' tSample strSyn ...
                            '_' num2str(nn) 'nn_' num2str(D) 'D.mat'];
                    end 
                    templateLabel=File.SwitchExtension(templateSample, ...
                        '.properties');
                end
            else
                samples=cell(1, N);
                for i=1:N
                    samples{i}=['s' num2str(sampleClues(i)) '_' ...
                        sampleSet '.csv'];
                    labels{i}=File.SwitchExtension(samples{i}, '.properties');
                end
                if ~isempty(templateClue)
                    if isnumeric(templateClue)
                        priorFldr='';
                        prefix='';
                        for i=1:length(templateClue)
                            prefix=[prefix 's' num2str(templateClue(i))];
                        end
                        tSample=[prefix '_' sampleSet];
                    else
                        if ~endsWith(templateClue, '.csv')
                            tSample=templateClue;
                        else
                            tSample=templateClue(1:end-4);
                        end
                        [priorFldr, tSample]=fileparts(tSample);
                        if ~isempty(strSyn)
                            if endsWith(tSample, strSyn)
                                strSyn='';
                            end
                        end
                    end
                    for i=1:nD
                        D=loDs(i);
                        templates{i}=['ust_' tSample strSyn '_' num2str(nn) ...
                            'nn_' num2str(D) 'D.mat'];
                    end
                    if isempty(priorFldr)
                        templateSample=[tSample '.csv'];
                    else
                        templateSample=fullfile(priorFldr, [tSample '.csv']);
                    end
                    templateLabel=File.SwitchExtension(templateSample, '.properties');
                end
            end
            hasTemplates=~isempty(templates) && ischar(templates{1});
            if hasTemplates
                [files, existence]=UmapUtil.RelocateExamples(...
                    [samples, labels templates templateSample ...
                    templateLabel], true, templates);
            else
                if size(samples) ~= size(labels)
                    samples=samples';
                end
                if isnumeric(templateClue) && length(templateClue)>1
                    [files, existence]=UmapUtil.RelocateExamples(...
                        [samples, labels templateSample ], true, ...
                        {templateSample});
                else
                    [files, existence]=UmapUtil.RelocateExamples(...
                        [samples, labels templateSample ]);
                end
            end
            if ~all(existence)
                existence2=existence;
                for i=1:nD
                    existence2((N*2)+i)=2;
                end
                if ~all(existence2)
                    sampleSet='';
                    samples={};
                    return;
                end
            end
            for j=1:N
                samples{j}=files{j};
                labels{j}=files{j+N};
            end
            if hasTemplates
                for i=1:nD
                    templates{i}=files{(N*2)+i};
                end
                templateSample=files{end-1};
                templateLabel=files{end};
            else
                templates={};
            end
        end        
        
        function Pair(templateClue, sampleClue, sampleSetClue, varargin)
            varargin=[{'fake.csv'} varargin];
            close all force;
            try
                Gui.DisposeAllJavaWindows;
            catch
            end
            [args, pu]=UstTest.Initialize('pair', [1 3 4], [1 3 4], varargin{:});
            if nargin<3
                sampleSetClue=[];
                if nargin<2
                    sampleClue=[2 3];
                    if nargin<1
                        templateClue=1;
                    end
                end
            end
            if isnumeric(templateClue) && length(templateClue)>1
                %get the UN-synthesized sample merger
                templateClue=UstTest.MergeSamples(sampleSetClue, ...
                    false, 0, pu, templateClue);
            end
            basicToo=args.ust_test_basic_reduction;
            nn=args.n_neighbors;
            loDs=args.ust_test_components;
            [sampleSet, samples, labels, templates, templateSample,...
                templateLabel, strSyn]=UstTest.SolveClues(sampleClue, ...
                templateClue, sampleSetClue, loDs, nn, ...
                args.ust_test_synthesize);
            if isempty(sampleSet)
                pu.close;
                return;
            end
            nD=length(loDs);
            N=length(samples);
            fldr=fileparts(templates{1});
            [ustFldr, testId]=UstTest.GetResultsFldr('pair', fldr, sampleSet);
            results=struct();
            tick=tic;
            strLoDs=String.Num2Str(loDs, '-');
            sfx=[strSyn '_' num2str(nn) 'nn_' strLoDs 'D.html'];
            if ischar(templateClue)
                matchFile=['matches_' num2str(testId) sfx];
            else
                matchFile=['matches_s' String.Num2Str(sampleClue,'-') ...
                    '_t' num2str(templateClue) sfx];
            end
            matchFile=fullfile(ustFldr, matchFile);
            fpnFile=File.SwitchExtension(matchFile, '.txt');
            if exist(fpnFile, 'file') && ...
                    (ismember(4, args.match_scenarios) ||...
                    ismember(3, args.match_scenarios))
                [answer, cancelled]=Gui.Ask(...
                    'Prior false +/- results exist...',...
                    {'Redo UST classification ...', ...
                    'View +/- of prior classification'}, ...
                    'prevFalsePosNeg', 'Redo?', 2);
                if cancelled
                    pu.close;
                    return;
                elseif answer==2
                    pu.close;
                    FalsePositiveNegative.Plot([0 1], fpnFile)
                    return;
                end
            end
            
            matchTypes=args.match_supervisors;
            scenarios=args.match_scenarios;            
            sampleSetAbbrv=sampleSet;
            
            for i=1:nD
                template=templates{i};
                [~,templateFile]=fileparts(templateSample);
                if ~exist(template, 'file')
                    if args.ust_test_synthesize>0
                        pu.setText2(['Generating ' ...
                            String.encodePercent(args.ust_test_synthesize) ...
                            ' synthetic data']);
                        try
                            [dataSet, cn]=File.ReadCsv(templateSample);
                        catch ex
                            msg(Html.Wrap(['Could not open your '...
                                'csv file<br>&nbsp;&nbsp;&nbsp;'...
                                BasicMap.Global.smallStart templateSample ...
                                '<br><br><b>BECAUSE</b>:<br>&nbsp;'...
                                '&nbsp;&nbsp;&nbsp;'...
                                ex.message BasicMap.Global.smallEnd ...
                                '<hr>']), 0, 'center', 'Error...', ...
                                'error.png');
                            pu.close
                            return;
                        end
                        synRows=floor(args.ust_test_synthesize*size(dataSet,1));
                        [syntheticData, syntheticLabels]=...
                            generate_synthetic_set(...
                            dataSet(:,1:end-1), dataSet(:,end), ...
                            synRows);
                        run_umap([syntheticData syntheticLabels], ...
                            'parameter_names', cn, ...
                            'label_file', templateLabel, ...
                            'label_column', 'end', ...
                            'n_neighbors', nn,...
                            'n_components', loDs(i),...
                            'color_file', args.color_file,...
                            'save_template_file', template);

                    else
                        run_umap(templateSample, ...
                            'label_file', templateLabel, ...
                            'label_column', 'end', ...
                            'n_neighbors', nn,...
                            'n_components', loDs(i),...
                            'color_file', args.color_file,...
                            'save_template_file', template);
                    end
                end
                for j=1:N
                    sample=samples{j};
                    args.label_file=labels{j};
                    [~,sampleFile]=fileparts(sample);
                    try
                        [dataSet, cn]=File.ReadCsv(sample);
                    catch ex
                        msg(Html.Wrap(['Could not open csv file<br>&nbsp;&nbsp;&nbsp;'...
                            BasicMap.Global.smallStart sample ...
                            '<br><br><b>BECAUSE</b>:<br>&nbsp;&nbsp;&nbsp;&nbsp;' ex.message ...
                            BasicMap.Global.smallEnd '<hr>']), 0, ...
                            'center', 'Error...', 'error.png');
                        ok=false;
                        break;
                    end
                    [R, hiD]=size(dataSet);
                    args.csv_file_or_data=dataSet;
                    if i==1 && j==1
                        UmapUtil.InitProgress(args, pu, N, 1, R, hiD-1);
                        if ~isempty(sampleSetAbbrv)
                            suffix=['_' num2str(hiD-1) 'D']; 
                            loi=String.LastIndexOf(sampleSetAbbrv, suffix);
                            if loi>0
                                sampleSetAbbrv=sampleSetAbbrv(1:loi-1);
                            end
                        end
                    end
                    reductionSummary=[String.encodeInteger(R)...
                        'x' num2str(hiD-1) 'D-' num2str(loDs(i)) 'D'] ;
                    if ~isempty(sampleSetAbbrv) && ...
                            contains(sampleFile, sampleSetAbbrv)
                        loi=String.IndexOf(sampleFile, '_');
                        if loi>0
                            sampleFileAbbrv=sampleFile(1:loi-1);
                        else
                            sampleFileAbbrv=sampleFile;
                        end
                        loi=String.IndexOf(templateFile, '_');
                        if loi>0
                            templateFileAbbrv=templateFile(1:loi-1);
                        else
                            templateFileAbbrv=templateFile;
                        end
                        args.sample_set=sampleSetAbbrv;
                        args.test_set=sampleFileAbbrv;
                        args.training_set=templateFileAbbrv;
                        args.description=sprintf(...
                            'Pair ust ss=%s;t=%s;nn=%d;s=%s %s', ...
                            sampleSetAbbrv, templateFileAbbrv, ...
                            nn, sampleFileAbbrv, reductionSummary);
                    else
                        args.description=sprintf(...
                            'Pair ust t=%s;nn=%d;s=%s %s', ...
                            templateFile, nn, sampleFile, reductionSummary);
                    end
                    args.template_file=template;
                    args.match_supervisors=matchTypes;
                    args.match_scenarios=scenarios;
                    args.parameter_names=cn(1:end-1);
                    [~, ~, ~, extras]=run_umap(args);
                    [results, ok]=UmapUtil.CollectResults(extras, results);
                    if ok && basicToo
                        args.description=sprintf(...
                            'Pair basic s=%s;nn=%d;%s', sampleFile, ...
                            nn, reductionSummary);
                        args.template_file='';
                        args.match_supervisors=1;
                        args.match_scenarios=scenarios(scenarios>2);
                        args.parameter_names={};
                        [~, ~, ~, extras]=run_umap(args);
                        if i==1 && j==1
                            results=UmapUtil.ExtendResultsHtmlHead(extras, results);
                        end
                        [results, ok]=UmapUtil.CollectResults(extras, results);
                    end
                    results=UmapUtil.SetResultsHtmlRow(results);
                    if ~ok
                        break;
                    end
                end
            end
            if ok
                UmapUtil.SeeHtml(results, matchFile, 1, tick, ...
                    [sampleSet ' ' matchFile]);
                UmapUtil.SaveFiles(results, matchFile, true);
            end
            pu.close;
        end
      
       
    end
end
