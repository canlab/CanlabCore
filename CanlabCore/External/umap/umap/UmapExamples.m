%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

classdef UmapExamples
    methods(Static)
        function DisplayPub(which)
            disp('PLEASE NOTE:The sample data for this example is made public with the publication');
            fprintf('%s\n\n', UmapExamples.SampleURL(which));
        end
        
        function url=SampleURL(which)
            switch which
                case 1
                    url='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4896314/';
                case 2
                    url='https://www.frontiersin.org/articles/10.3389/fimmu.2019.01194/full';
                otherwise
                    url='https://www.pnas.org/content/107/6/2568';
                        
                    
            end
        end
        function sample=SolveClue(sampleClue)
            if isnumeric(sampleClue)
                switch sampleClue
                    case 27
                        UmapExamples.DisplayPub(2);
                        sample='genentech_d1_27D.csv';
                    case 29
                        UmapExamples.DisplayPub(1);
                        sample='s1s2_samusikImported_29D.csv';
                    case 11
                        UmapExamples.DisplayPub(3);
                        
                        sample='balbc_Bcells_Macrophages.csv';
                    case 11.12
                        UmapExamples.DisplayPub(3);
                        sample='sampleBalbcLabeled12k.csv';
                    
                    otherwise
                        error('Do not understand sample clue %d',...
                            sampleClue );
                end
            else
                sample=sampleClue;
            end
        end
        
        function SplitSlow(sampleClue)
           if nargin<1
               sampleClue=29;
           end
            UstTest.Split(...
                UmapExamples.SolveClue(sampleClue),...
                'ust_test_components', [3 2], ...
                'ust_test_basic_reduction', true, ...
                'cluster_detail', {'most high', 'low'}, ...
                'match_html', 1, 'ust_test_cases', [ 1 4]);
        end
        
        function SplitFast(sampleClue)
            if nargin<1
               sampleClue=29;
           end
            UstTest.Split(...
                UmapExamples.SolveClue(sampleClue),...
               'ust_test_cases', 1, 'ust_test_basic_reduction', true, ...
                'cluster_detail', {'most high', 'low'}, 'match_html', 1);
        end

        function SplitFastest(sampleClue)
            if nargin<1
                sampleClue=29;
            end
            UstTest.Split(...
                UmapExamples.SolveClue(sampleClue),...
                'cluster_detail', 'most high',...
                'ust_test_components', 3,...
                'ust_test_cases', 1, 'ust_test_basic_reduction', false, ...
                'match_html', 1, 'match_supervisors', 3 , ...
                'match_scenarios', 4);
        end

        %The sample data for this example is made public with the publication
        %https://www.frontiersin.org/articles/10.3389/fimmu.2019.01194/full
        function Genentech(clusterDetail, scenarios, matchTypes)
            if nargin<3
                matchTypes=1;
                if nargin<2
                    scenarios=[];
                    if nargin<1
                        clusterDetail={};
                    end
                end
            end
            if isempty(scenarios)
                scenarios=[1 4];
            end
            if isempty(clusterDetail)
                clusterDetail={'most high', 'low'};
            end
            if ischar(clusterDetail)
                clusterDetail={clusterDetail};
            end
            UmapExamples.DisplayPub(2);
            run_umap('genentech_d1_27D_70k.csv', ...
                'label_file', 'genentech_d1_27D.properties', ...
                'label_column', 'end', ...
                'template_file', 'genentech_d1_27D_3D.mat',...
                'see_training', true, ...
                'cluster_detail', clusterDetail, ...
                'color_file', 'colorsByName.properties',...
                'match_supervisors', matchTypes, ...
                'match_histogram_fig', false,...
                'match_html', -1, 'match_scenarios', scenarios);
        end

        
        function Samusik(clusterDetail, scenarios, matchTypes, ...
                sampleSetClue, sampleClue, templateClue, nn, D)
            if nargin<8
                D=3;
                if nargin<7
                    nn=15;
                    if nargin<6
                        templateClue=1;
                        if nargin<5
                            sampleClue=2;
                            if nargin<4
                                sampleSetClue=29;
                                if nargin<3
                                    matchTypes=3;
                                    if nargin<2
                                        scenarios=[];
                                        if nargin<1
                                            clusterDetail={};
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if isempty(scenarios)
                scenarios=[1 4];
            end
            if isempty(clusterDetail)
                clusterDetail='most high';
            end
            if ischar(clusterDetail)
                clusterDetail={clusterDetail};
            end
            [tmplt, ss]=UstTest.Supervise(sampleSetClue, templateClue, ...
                true, D, nn);
            smpl=['s' num2str(sampleClue) '_' ss ];
            rf=fullfile(UmapUtil.LocalSamplesFolder, 'ustTest.pair');
            mf=sprintf('%s_s%d_t%d__mt%s__ms%s_%dnn_%dD', ss, ...
                sampleClue, templateClue, String.Num2Str(matchTypes), ...
                String.Num2Str(scenarios),  nn, D);
            run_umap([smpl '.csv'], ...
                'label_file', [smpl '.properties'], ...
                'label_column', 'end', ...
                'template_file', tmplt ,...
                'see_training', true, ...
                'cluster_detail', clusterDetail, ...
                'color_file', 'colorsByName.properties',...
                'match_supervisors', matchTypes, ...
                'match_histogram_fig', false,...
                'match_html', 1, ...
                'result_folder', rf,...
                'match_file', mf, ...
                'match_scenarios', scenarios);
        end

        %The sample data for this example is made public with the publication
        %https://www.pnas.org/content/107/6/2568
        function Balbc12(clusterDetail, scenarios, matchTypes)
            if nargin<3
                matchTypes=1;
                if nargin<2
                    scenarios=[];
                    if nargin<1
                        clusterDetail={};
                    end
                end
            end
            if isempty(scenarios)
                scenarios=4;
            end
            if isempty(clusterDetail)
                clusterDetail={'most high'};
            end
            if ischar(clusterDetail)
                clusterDetail={clusterDetail};
            end
            rf=fullfile(UmapUtil.LocalSamplesFolder, 'ustTest.pair');
            [~,~,~,extras]=run_umap('sampleBalbcLabeled12k.csv', ...
                'label_file', 'balbcLabels.properties', ...
                'label_column', 'end', ...
                'template_file', 'ustBalbc3D.mat',...
                'see_training', true, ...
                'cluster_detail', clusterDetail, ...
                'color_file', 'colorsByName.properties',...
                'match_supervisors', matchTypes, ...
                'match_histogram_fig', false,...
                'result_folder', rf, ...
                'match_html', -1, ...
                'match_scenarios', scenarios);
            records=extras.getFalsePosNeg;
            if ~isempty(records)
                extras.getFalsePosNeg(true);
                extras.saveMatchFiles;
                FalsePositiveNegative.Plot([0 1], extras.falsePosNegFile);
            end
        end

        %The sample data for this example is made public with the publication
        %https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4896314/
        function PairSamusikManual(D, nn)
            if nargin<2
                D=3;
                if nargin<1
                    nn=15;
                end
            end
            UmapExamples.DisplayPub(1);
            UstTest.Pair(1, [2 10], 29.2, 'verbose', 'text', 'match_supervisors', 1:4, 'ust_test_basic_reduction', false, 'match_table_fig', true, 'match_scenarios', [1 2 3 4], 'ust_test_components', D, 'n_neighbors', nn);
        end

        %The sample data for this example is made public with the publication
        %https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4896314/
        function PairSamusikImport(nn, D, testSets, synthesize, templateClue)
            if nargin<5
                templateClue=1;
                if nargin<4
                    synthesize=0;
                    if nargin<3
                        testSets=[2 10];
                        if nargin<2
                            D=3;
                            if nargin<1
                                nn=30;
                            end
                        end
                    end
                end
            end
            UmapExamples.DisplayPub(1);
            UstTest.Pair(templateClue, testSets, 29, ...
                'ust_test_synthesize', synthesize, ...
                'verbose', 'text', ...
                'match_supervisors', 3, ...
                'ust_test_basic_reduction', false, ...
                'match_scenarios', 4, ...
                'ust_test_components', D, 'n_neighbors', nn);
        end

        %The sample data for this example is made public with the publication
        %https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4896314/
        function PairSamusikImport19(nn, D)
            if nargin<2
                D=3;
                if nargin<1
                    nn=15;
                end
            end
            UmapExamples.DisplayPub(1);
            UstTest.Pair(1, [2 10], 19, 'verbose', 'text', 'match_supervisors', 3, 'ust_test_basic_reduction', false, 'match_table_fig', true, 'match_scenarios', 4, 'ust_test_components', D, 'n_neighbors', nn);
        end

        %The sample data for this example is made public with the publication
        %https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4896314/
        function PairSamusikNn(nn)
            if nargin<1
                nn=20;
            end
            UmapExamples.DisplayPub(1);
            UstTest.Pair(1, 2, 29, 'verbose', 'graphic', 'match_supervisors', 1:4, 'ust_test_basic_reduction', true, 'match_table_fig', true, 'match_scenarios', 4, 'ust_test_components', 3, 'n_neighbors', nn);
        end

        %The sample data for this example is made public with the publication
        %https://www.pnas.org/content/107/6/2568
        function SplitBalbc55
            UmapExamples.DisplayPub(3);
            UstTest.Split('sampleBalbcLabeled55k.csv',...
                'ust_test_cases', 1, 'ust_test_basic_reduction', false, ...
                'cluster_detail', {'most high', 'low'}, 'match_html', 1, ...
                'match_scenarios', [ 1 3 4 ], 'match_supervisors', [1 3]);
        end
        
        %The sample data for this example is made public with the publication
        %https://www.pnas.org/content/107/6/2568
        function SplitBalbc55withBasic
            UmapExamples.DisplayPub(3);
            UstTest.Split('sampleBalbcLabeled55k.csv',...
                'ust_test_cases', 1, 'ust_test_basic_reduction', false, ...
                'cluster_detail', {'most high', 'low'}, 'match_html', 1, ...
                'match_scenarios', [ 1 3 4 ], 'ust_test_basic_reduction', ...
                true);
        end

        %The sample data for this example is made public with the publication
        %https://www.pnas.org/content/107/6/2568
        function Split5
            UmapExamples.DisplayPub(3);
            UstTest.Split('sampleBalbcLabeled55k.csv',...
                'ust_test_cases', 1, 'ust_test_basic_reduction', true, ...
                'cluster_detail', {'most high', 'low'}, 'match_html', 1, ...
                'match_scenarios', [ 1 2 3 4 ], 'ust_test_components', [3 2]);
        end
        
       %The sample data for this example is made public with the publication
        %https://www.pnas.org/content/107/6/2568
        function EliverSplitFast
            UmapExamples.DisplayPub(3);
            UstTest.Split('sampleBalbcLabeled12k.csv',...
                'ust_test_cases', 4, 'ust_test_basic_reduction', false, ...
                'cluster_detail', 'most high', 'match_html', 1, ...
                'match_supervisors', [3 4], ...
                'match_scenarios', 4, 'ust_test_basic_reduction', ...
                false, 'ust_test_components', 2);
        end

        function Split6
            UmapExamples.DisplayPub(3);
            UstTest.Split('sampleBalbcLabeled55k.csv',...
                'ust_test_cases', 1, 'ust_test_basic_reduction', false, ...
                'cluster_detail', {'most high', 'low'}, 'match_html', 1, ...
                'match_scenarios', [ 1 4 ], 'ust_test_basic_reduction', ...
                true, 'ust_test_components', [3 2],'ust_test_cases',...
                3);
        end
        function SplitGenenTech(synthesize, matchTypes)
            if nargin<2
                matchTypes=[3 4];
                if nargin<1
                    synthesize=0;
                end
            end
            UmapExamples.DisplayPub(2);
            UstTest.Split('genentech_d1_27D.csv', ...
                'ust_test_synthesize', synthesize, ...
                'ust_test_components', 2, ...
                'ust_test_basic_reduction', false, ...
                'cluster_detail', 'most high', ...
                'match_html', 1, ...
                'n_neighbors', 30, ...
                'match_supervisors', matchTypes, ...
                'match_scenarios', 4, ...
                'ust_test_cases', [ 1 2 4]);
        end

        
        function RunTest1
            UstTest.RunBig(10);
        end

        
    end
end