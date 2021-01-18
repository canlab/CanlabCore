%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

classdef UMAP_extra_results < handle
    properties
        supervisorMatchedLabels; 
                        %vector of labels matched from a supervised template 
                        % 1 value per row of input data matrix
        fig;            % main figure with UMAP output plot
        qft;            % instance of qf tree object for input data. 
                        % The fig property contains the view of the tree.
        qftSupervisors; % instance of qf tree object for supervisor data. 
                        % The fig property
                        % contains the view of the tree.
        qfd={};         % array qf QFTable instances. 
                        % The function getAverages summarizes the result.
                        % Relevant figure properties are fig for table view, 
                        % qHistFig for view of dismilarity score histogram
                        % fHistFig for view of f measure histogram
       matchHtmlHead1;
       matchHtmlHead2;
       matchHtmlHead3;
       matchHtmlBody;
       matchCsvHead;
       matchCsvBody;
       falsePosNegFile;
       falsePosNegHead;
       falsePosNegBody='';
       app;
       htmlFile;
       csvFile;
       args;
    end
    
    methods
        
        function closeMatchFigs(this)
            N=length(this.qfd);
            for i=1:N
                this.qfd{i}.closeFigs;
            end
        end
        
        function closeTreeFigs(this)
            if ~isempty(this.qft)
                Gui.CloseFig(this.qft.fig);
            end
            if ~isempty(this.qftSupervisors)
                Gui.CloseFig(this.qftSupervisors.fig);
            end
        end
        
        function avgs=getMatchAverages(this)
            N=length(this.qfd);
            avgs=[];
            for i=1:N
                qf=this.qfd{i};
                [medianDissimilarity, meanDissimilarity, ...
                    medianOverlap, meanOverlap]=qf.getAverages;
                avgs(end+1,:)=[medianDissimilarity meanDissimilarity...
                    medianOverlap meanOverlap];
            end
        end
        
        function fig=seeFalsePosNeg(this)
            if ~isempty(this.falsePosNegFile) && ~isempty(this.falsePosNegHead)
                File.SaveTextFile(this.falsePosNegFile, ...
                    [this.falsePosNegHead this.falsePosNegBody])
                [fig, that]=FalsePositiveNegative.Plot([0 1], ...
                    this.falsePosNegFile);
                that.fcnMoreHtml=@()getFalsePosNegMatrixHtml(this);
            else
                fig=[];
            end
        end
        
        function saveMatchFiles(this, h1)
            if nargin<2
                h1='';
            end
            if ~isempty(this.csvFile) && ~isempty(this.matchCsvHead)
                File.SaveTextFile(this.csvFile, [this.matchCsvHead ...
                    this.matchCsvBody])
            end
            if ~isempty(this.falsePosNegFile) && ~isempty(this.falsePosNegHead)
                File.SaveTextFile(this.falsePosNegFile, ...
                    [this.falsePosNegHead this.falsePosNegBody])
            end
            if ~isempty(this.htmlFile) && ~isempty(this.matchHtmlHead1)
                File.SaveTextFile(this.htmlFile, this.getMatchHtml(h1));
            end
        end
        
        function seeMatches(this, how, h1)
            if nargin<3
                h1='';
            end
            if how==-1
                msg(this.getMatchHtml(h1), 0, 'east++', ['Matches: '...
                    this.describeRedution ' reduction'], 'none')
            elseif how==1
                Html.BrowseString(this.getMatchHtml(h1));
            elseif how==2
                Html.BrowseFile(this.htmlFile);
            end
        end
        
        function this=UMAP_extra_results()
            this.app=BasicMap.Global;
        end
        
        function doMatchOutput(this, D)
            if ~isempty(this.qfd)
                this.doMatchHtmlHead(D);
                this.doMatchHtmlBody;
                this.doMatchCsv;
                this.doFalsePosNeg;
            end
        end
        
        function htmls=getFalsePosNegMatrixHtml(this)
            args_=this.args;
            qfds=this.qfd;
            N=length(qfds);
            if N<1 || isempty(find(args_.match_scenarios==4, 1))
                return;
            end
            reduction=args_.reduction;
            if ~isempty(args_.training_set)
                trainingSet=args_.training_set;
            elseif ischar(args_.template_file)
                [~,fn]=fileparts(args_.template_file);
                trainingSet=fn;
            else
                trainingSet='';
            end
            if ~isempty(args_.test_set)
                testSet=args_.test_set;
            elseif ischar(args_.csv_file_or_data)
                [sampleFldr,fn]=fileparts(args_.csv_file_or_data);
                testSet=fn;
            else
                testSet='';
            end
            if ~isempty(args_.sample_set)
                sampleSet=args_.sample_set;
            elseif exist('sampleFldr', 'var')
                [~, sampleSet]=fileparts(sampleFldr);
            else
                sampleSet='';
            end
            htmls=cell(1,N);
            for i=1:N
                qt=qfds{i};
                N2=length(qt.qf.falsePosNegs);
                if N2>0
                    mt=qt.context.matchType;
                    cd=qt.context.clusterDetail;
                    cd=StringArray.IndexOf(Density.DETAILS, cd);
                    htmls{i}=FalsePositiveNegative.MatrixHtml(qt.qf);
                end
            end
        end
        
        function doFalsePosNeg(this, doPlotNow)
            args_=this.args;
            qfds=this.qfd;
            N=length(qfds);
            if N<1 || isempty(find(args_.match_scenarios==4, 1))
                return;
            end
            fullBody=this.falsePosNegBody;
            this.falsePosNegFile=[UMAP_extra_results.FileName(...
                args_, 'falsePosNeg') '.txt'];            
            reduction=args_.reduction;
            if ~isempty(args_.training_set)
                trainingSet=args_.training_set;
            elseif ischar(args_.template_file)
                [~,fn]=fileparts(args_.template_file);
                trainingSet=fn;
            else
                trainingSet='';
            end
            if ~isempty(args_.test_set)
                testSet=args_.test_set;
            elseif ischar(args_.csv_file_or_data)
                [sampleFldr,fn]=fileparts(args_.csv_file_or_data);
                testSet=fn;
            else
                testSet='';
            end
            if ~isempty(args_.sample_set)
                sampleSet=args_.sample_set;
            elseif exist('sampleFldr', 'var')
                [~, sampleSet]=fileparts(sampleFldr);
            else
                sampleSet='';
            end
            for i=1:N
                qt=qfds{i};
                N2=length(qt.qf.falsePosNegs);
                if N2>0
                    mt=qt.context.matchType;
                    cd=qt.context.clusterDetail;
                    cd=StringArray.IndexOf(Density.DETAILS, cd);
                    [body, notFound]=FalsePositiveNegative.TabRows(...
                        qt.qf.falsePosNegs, reduction, sampleSet, ...
                        trainingSet, testSet, args_.n_neighbors, ...
                        args_.hiD, args_.n_components, mt, cd);
                    if ~isempty(notFound)
                        body=[body newline notFound];
                    end
                    fullBody=[fullBody body newline];
                end
            end
            this.falsePosNegHead=FalsePositiveNegative.TabHead;
            this.falsePosNegBody=fullBody;
            if nargin>1 && doPlotNow
                this.seeFalsePosNeg;
            end
        end
        
        
        function html=getMatchHtml(this, h1)
            if nargin<2
                h1='';
            end
            html=['<html>' h1 '<table border="1"><thead>' ...
                '<tr>' this.matchHtmlHead1 '</tr>' ...
                '<tr>' this.matchHtmlHead2 '</tr>' ...
                '<tr>' this.matchHtmlHead3 '</tr></thead>'...
                '<tr>' this.matchHtmlBody '</tr>'...
                '</table></html>'];
        end
        
        function s=describeRedution(this)
            s=UmapUtil.GetReductionLongText(this.args.reductionType);
        end
        
        function records=getFalsePosNeg(this, verbose)
            N=length(this.qfd);
            records={};
            for i=1:N
                qf=this.qfd{i};
                N2=length(qf.qf.falsePosNegs);
                if N2>0
                    if nargin>1 && verbose
                        for j=1:N2
                            record=qf.qf.falsePosNegs{j};
                            fprintf(['#%d "%s"; false pos=%s '...
                                'false neg=%s; match 1 x %d\n'], ...
                                j, record.trainingClass, ...
                                String.encodePercent(record.falsePos, ...
                                record.testSize, 2), ...
                                String.encodePercent(record.falseNeg, ...
                                record.trainingSize, 2), ...
                                length(record.testIds));
                        end
                    end
                    records=[records qf.qf.falsePosNegs];
                end
            end
        end
    end
    
    methods(Access=private)    
        function doMatchHtmlHead(this, D)
            args_=this.args;
            nMatchTypes=length(args_.match_supervisors);
            nMatchScenarios=length(args_.match_scenarios);
            if any(args_.match_scenarios==1)
                nUstMatches=nMatchScenarios-1;
            else
                nUstMatches=nMatchScenarios;
            end
            this.htmlFile=[UMAP_extra_results.FileName(args_) '.html'] ;
            colspan=['colspan="' num2str(nUstMatches*3) '"'];
            if isempty(args_.description)
                htmlDsc='<th></th>';
            else
                htmlDsc='<th>Context</th>';
            end
            sm1=this.app.smallStart;
            sm2=this.app.smallEnd;
            html1='';
            html2='';
            html3='';
            has1=false;
            cluDtls=args_.cluster_detail;
            nCluDtls=length(cluDtls);
            if length(args_.ust_test_components)>1
                loD=[];
            else
                loD=args_.n_components;
            end
            for c=1:nCluDtls
                for i=1:nMatchTypes
                    matchType=args_.match_supervisors(i);
                    if c>1
                        if matchType>=3 % nearest neighbor no clustering
                            continue;
                        end
                    end
                    mt=UmapUtil.GetMatchTypeLongText(matchType, ...
                        args_.reductionType, loD, D);
                    if matchType<3
                        mt=[mt ' ' cluDtls{c}];
                    end
                    html1=[html1 '<th ' colspan '>' mt '</th>'];
                    for j=1:nMatchScenarios
                        scenario=args_.match_scenarios(j);
                        sc=UmapUtil.GetMatchScenarioText(scenario, args_.reductionType) ;
                        if scenario==1
                            if ~has1
                                has1=true;
                                html1_1=['<th colspan="3">Prior<br>'...
                                    sm1 'classifications' sm2 '</th>' ];
                                html2_1=['<th colspan="3">' sc '<br>' ...
                                    sm1 'Dissimilarity' sm2 '</th>'];
                                html3_1='<th>Unmatched</th><th>Median</th><th>Mean</th>';
                            end
                        else
                            if scenario==4
                                htmlStat=['<br>' sm1 'Overlap' sm2];
                            else
                                htmlStat=['<br>' sm1 'Dissimilarity' sm2];
                            end
                            html2=[html2 '<th colspan="3">' sc htmlStat '</th>' ];
                            html3=[html3 '<th>Unmatched</th><th>Median</th><th>Mean</th>'];
                        end
                    end
                end
            end
            if has1
                if nUstMatches>0
                    this.matchHtmlHead1=[htmlDsc html1_1  html1];
                    this.matchHtmlHead2=['<th></th>' html2_1 html2];
                    this.matchHtmlHead3=['<th></th>' html3_1 html3];
                else
                    this.matchHtmlHead1=[htmlDsc html1_1];
                    this.matchHtmlHead2=['<th></th>' html2_1];
                    this.matchHtmlHead3=['<th></th>' html3_1];
                end
            else
                this.matchHtmlHead1=[htmlDsc html1];
                this.matchHtmlHead2=['<th></th>'  html2];
                this.matchHtmlHead3=['<th></th>' html3]; 
            end
        end
        

        function doMatchHtmlBody(this)
            args_=this.args;
            qfds=this.qfd;
            dsc=args_.description;
            td='<td align="right">';
            td_='</td>';
            N_=length(qfds);
            htmlDsc=['<td>' dsc td_];
            html='';
            mdn=zeros(1, N_);
            mn=zeros(1,N_);
            adjMdn=zeros(1, N_);
            adjMn=zeros(1,N_);
            matchRate=zeros(1,N_);
            matches=cell(1, N_);
            scenarios=zeros(1,N_);
            isDis=true(1, N_);
            for i=1:N_
                qf=qfds{i};
                matchStrategy=qf.context.matchStrategy;
                [tUnmatched, tN, sUnmatched, sN]=qf.getMatchCounts;
                matches{i}=sprintf('%d/%d, %d/%d', tUnmatched, tN,...
                    sUnmatched, sN);
                if matchStrategy==1
                    scenarios(i)=qf.context.matchScenario;
                    [~,mdn(i), mn(i)]=qf.getData(true);
                else
                    scenarios(i)=4;
                    [~,mdn(i), mn(i)]=qf.getData(false);
                    isDis(i)=false;
                end
                tM=tN-tUnmatched;
                sM=sN-sUnmatched;
                matchedSubsetCnt=tM+sM;
                subsetCnt=tN+sN;
                
                if tUnmatched>0 || sUnmatched>0
                    if matchStrategy==1
                        addDs=(tUnmatched+sUnmatched)*100;
                        adjMdn(i)=(matchedSubsetCnt*mdn(i)+addDs)/subsetCnt;
                        adjMn(i)=(matchedSubsetCnt*mn(i)+addDs)/subsetCnt;
                    else
                        adjMdn(i)=matchedSubsetCnt*mdn(i)/subsetCnt;
                        adjMn(i)=matchedSubsetCnt*mn(i)/subsetCnt;
                    end
                else
                    adjMdn(i)=mdn(i);
                    adjMn(i)=mn(i);
                end
                matchRate(i)=matchedSubsetCnt/subsetCnt;
            end
            usc=unique(scenarios);
            nSc=length(usc);
            bestMedian=zeros(1, 4);
            bestMean=zeros(1, 4);
            bestMatches=zeros(1,4);
            worstMatches=zeros(1,4);
            for i=1:nSc
                u=usc(i);
                if u==4
                    bestMedian(u)=max(adjMdn(scenarios==u));
                    bestMean(u)=max(adjMn(scenarios==u));
                else
                    bestMedian(u)=min(adjMdn(scenarios==u));
                    bestMean(u)=min(adjMn(scenarios==u));
                end
                bestMatches(u)=max(matchRate(scenarios==u));
                worstMatches(u)=min(matchRate(scenarios==u));
            end
            has1=false;
            for i=1:N_
                scenario=scenarios(i);
                sMdn=String.encodeRounded(mdn(i), 1, true);
                sMn=String.encodeRounded(mn(i), 1, true);
                match=matches{i};
                ms1='';
                ms2='';

                if scenario>1 
                    if matchRate(i)==bestMatches(scenario)
                        ms1='<b><font color="#CC00FF">';
                        ms2='</font></b>';
                    elseif matchRate(i)==worstMatches(scenario)
                        ms1='<i><font color="red">';
                        ms2='</font></i>';
                    end
                end
                if scenario>1 && adjMdn(i)==bestMedian(scenario)
                    s1='<b><font color="blue">';
                    s2='</font></b>';
                else
                    s1='';
                    s2='';
                end
                if scenario==1
                    html1=[td ms1 match ms2 td_ td s1 sMdn  s2 td_];
                else
                    html=[html td ms1 match ms2 td_ td s1 sMdn  s2 td_];
                end
                if scenario>1 && adjMn(i)==bestMean(scenario)
                    s1='<b><font color="#00FFFF">';
                    s2='</font></b>';
                else
                    s1='';
                    s2='';
                end
                if scenario==1
                    has1=true;
                    html_1=[html1 td s1 sMn s2 td_ ];
                else
                    html=[html td s1 sMn s2 td_ ];
                end
            end
            if has1
                html=[htmlDsc html_1 html];
            else
                html=[htmlDsc html];
            end
            this.matchHtmlBody=html;
        end
        
        
        function doMatchCsv(this)
            args_=this.args;
            qfds=this.qfd;
            N=length(qfds);
            if N>0
                this.csvFile=[UMAP_extra_results.FileName(args_) '.csv'] ;
                lf=newline;
                this.matchCsvHead=[ StringArray.toString(fieldnames(...
                    qfds{1}.context), ',') 'trainingUnmatched,'...
                    'trainingSubsets,testUnMatched,testSubsets,'...
                    'median,mean' lf];
                csv='';
                for i=1:N
                    qf=qfds{i};
                    [tUnmatched, tN, sUnmatched, sN]=qf.getMatchCounts;
                
                    if qf.context.matchStrategy==1
                        [~,mdn, mn]=qf.getData(true);
                    else
                        [~,mdn, mn]=qf.getData(false);
                    end
                    csv=[csv String.toString(qf.context) ...
                        num2str(tUnmatched) ',' num2str(tN) ',' ...
                        num2str(sUnmatched) ',' num2str(sN) ',' ...
                        num2str(mdn) ',' num2str(mn) ',' lf];
                end
                this.matchCsvBody=csv;
            end
        end
    end
    
    methods(Static)
        function file=FileName(args, prefix)
            if nargin<2
                prefix='match';
            end
            if isfield(args, [prefix '_file'])
                fn=getfield(args, [prefix '_file']);
            else
                fn='';
            end
            if isempty(fn)
                strMatchTypes=strrep(strrep(num2str(...
                    args.match_supervisors), ' ', '_'), '__', '_');
                strMatchScenarios=strrep(strrep(num2str(...
                    args.match_scenarios), ' ', '_'), '__', '_');
                fn=[prefix '_' strMatchTypes...
                    '_for_' strMatchScenarios ];
            end
            file=fullfile(args.result_folder, fn);
            File.mkDir(args.result_folder);
        end
    end
end