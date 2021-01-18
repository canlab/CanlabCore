%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef LabelBasics
    methods(Static)
        function [labels, loss]=RemoveOverlap(labels)
            nCols=size(labels,2);
            if nCols>1
                y=labels(:,1);
                nCols=size(labels,2);
                
                for col=2:nCols
                    ll=y==0;
                    more=labels(ll, col);
                    y(ll)=more;
                end
                if nargout>1
                    lbls=unique(labels);
                    lbls=lbls(lbls>0);
                    h1=histc(sort(labels(:)), lbls);
                    h2=histc(sort(y), lbls);
                    loss=h2-h1;
                end
                labels=y;
            end
        end
        
        % assumes AutoGate modules are on path
        function map=GetLabelMap(gt, labels)
            ids=unique(labels);
            [names, N, strIds]=GatingTree.GetUniqueNames(gt.tp, ids);
            map=java.util.Properties;
            hl=gt.highlighter;
            for i=1:N
                dfltClr=Gui.HslColor(i,N);
                if ids(i)~=0
                    id=strIds{i};
                    name=names{i};
                    if isempty(name)
                        name=['label=' id];
                    end
                    map.put(id, name);
                    clr=num2str(floor(hl.getColor(id, dfltClr)*256));
                    map.put([id '.color'], clr);
                end
            end
        end
        
    
        function [key, keyColor, keyTraining, keyTest]=Keys(lbl)
            key=num2str(lbl); % prevent java.lang.Character
            keyColor=[key '.color'];
            keyTraining=[key '.trainingFrequency'];
            keyTest=[key '.testFrequency'];
            key=java.lang.String(key);
        end
        
        function Frequency(lbls, lblMap, training, match, how)
            if nargin<5
                how=-1; %contains, 0=eq, 1 = startsWith
                if nargin<4
                    match='';
                    if nargin<3
                        training=true;
                    end
                end
            end
            [trainFreq, testFreq]=LabelBasics.HasFrequencies(lblMap);
            total=length(lbls);
            tab=sprintf('\t');
            u=unique(lbls);
            nLbls=length(u);
            frequencies=histc(sort(lbls), u);
            [sortFrequencies,II]=sort(frequencies, 'descend');
            if nargin>3
                if ~isempty(training) && training
                    disp([num2str(nLbls) ' training set labels']);
                else
                    disp([num2str(nLbls) ' test set labels']);
                end
            end
            for i=1:nLbls
                lbl=u(II(i));
                freq=sortFrequencies(i);
                pFreq=String.encodePercent(freq, total, 1);
                [key, ~, keyTraining, keyTest]=LabelBasics.Keys(lbl);
                if ~isempty(training)
                    if training
                        lblMap.put(keyTraining, pFreq);
                    else
                        lblMap.put(keyTest, pFreq);
                    end
                end
                if nargin>3
                    name=lblMap.get(key);
                    if trainFreq && ~training
                        start=['From ' lblMap.get(keyTraining) ' to '];
                    else
                        start='';
                    end
                    if ~isempty(match)
                        start=[tab start];
                        if how==-1
                            if String.Contains(name, match)
                                start=['* ' start];
                            end
                        elseif how==1
                            if String.StartsWith(name, match)
                                start=['* ' start];
                            end
                        else
                            if strcmp(name, match)
                                start=['* ' start];
                            end
                        end
                    end
                    sFreq=String.encodeInteger(freq);
                    fprintf(['%s%s=%s events for #%d="%s"\n'], start, ...
                        pFreq, sFreq, lbl, name);
                end
            end
            if nargin>3
                if ~isempty(training) && training
                    disp([num2str(nLbls) ' training set labels']);
                else
                    disp([num2str(nLbls) ' test set labels']);
                end
            end
            if ~isempty(training)
                if training
                    lblMap.put('hasTrainingFrequencies', 'yes');
                else
                    lblMap.put('hasTestFrequencies', 'yes')
                end
            end
        end
        
        function [training, test]=HasFrequencies(lblMap)
            if isempty(lblMap)
                training=false;
                test=false;
            else
                training=strcmpi('yes', lblMap.get('hasTrainingFrequencies'));
                test=strcmpi('yes', lblMap.get('hasTestFrequencies'));
            end
        end
        
        function [addTrainingHtml, sup1, sup2, trStart, trEnd]=...
                AddTrainingHtml(lblMap, needHtml)
            addTrainingHtml=false;
            if needHtml
                sup1=BasicMap.Global.supStart;
                sup2=BasicMap.Global.supEnd;
                if ~isempty(lblMap)
                    trStart=[sup1 '<font color="#11DDAA"><b> training '];
                    trEnd=['</b></font>' sup2];
                    [trainFreq, testFreq]=LabelBasics.HasFrequencies(lblMap);
                    addTrainingHtml=trainFreq&&testFreq;
                else
                    trStart=[]; trEnd=[];
                end
            else
                sup1=[]; sup2=[]; trStart=[]; trEnd=[];
            end

        end
        
        function[data, columnNames1]=Merge2Samples(...
                fileName, sample1, label1, sample2, label2, label3)
            [data1, columnNames1]=File.ReadCsv(sample1);
            [data2, columnNames2]=File.ReadCsv(sample2);
            if ~isequal(columnNames1(1:end-1), columnNames2(1:end-1))
                msgWarning('Training & test set labels do not match');
                return;
            end
            lblMap2=loadLabels(label2);
            lblMap1=loadLabels(label1);
            if isequal(label3, label1)
                lbls=data2(:,end);
                reLbls=LabelBasics.RelabelIfNeedBe(lbls, lblMap1, lblMap2);
                if size(reLbls,2)~=size(lbls,2)
                    data=[];
                    return;
                end
                data2(:,end)=reLbls;
            else
                lblMap3=loadLabels(label3);
                lbls=data1(:,end);
                reLbls=LabelBasics.RelabelIfNeedBe(lbls, lblMap3, lblMap1);
                if size(reLbls,2)~=size(lbls,2)
                    data=[];
                    return;
                end
                data1(:,end)=reLbls;
                lbls=data2(:,end);
                reLbls=LabelBasics.RelabelIfNeedBe(lbls, lblMap3, lblMap2);
                if size(reLbls,2)~=size(lbls,2)
                    data=[];
                    return;
                end
                data2(:,end)=reLbls;
            end
            data=[data1;data2];
            if ~isempty(fileName)
                try
                    fu=edu.stanford.facs.wizard.FcsUtil(fileName);
                catch ex
                    msgError('Can not load java edu.stanford.facs.wizard.FcsUtil' );
                end
                problem=fu.createTextFile(fileName, [], data, [], columnNames1);
                copyfile(label1, File.SwitchExtension(fileName, '.properties')) 
                if ~isempty(problem)
                    msgError(problem, 12);
                end
            end
            
            function lblMap=loadLabels(lblFile)
                try
                    lblMap=java.util.Properties;
                    lblMap.load(java.io.FileInputStream(lblFile));
                catch ex
                    ex.getReport
                    lblMap=[];
                end
            end
        end
        
        function lbls=RelabelIfNeedBe(lbls, trainingMap, testMap)
            if isempty(trainingMap) || trainingMap.size()==0
                msgWarning('Training set map is empty');
                lbls=[];
                return;
            end
            if isempty(testMap) || testMap.size()==0
                msgWarning('Test set map is empty');
                lbls=[];
                return;
            end
            u=unique(lbls);
            N=length(u);
            trainingIdByName=LabelBasics.IdByName(trainingMap);
            for i=1:N
                lbl=u(i);
                if lbl ~= 0
                    key=java.lang.String(num2str(lbl));
                    if ~trainingMap.containsKey(key)
                        name=testMap.get(key);
                        if isempty(name)
                           warning(['Test set properties lack'...
                               ' name for label "' char(key) '"']); 
                        else
                            newLbl=trainingIdByName.get(name);
                            if ~isempty(newLbl)
                                newLbl=str2double(newLbl);
                                lbls(lbls==lbl)=newLbl;
                            else
                                warning(['Training set properties lack'...
                                    ' name for label "' char(key) '"'...
                                    'named "' name '"']);
                            end
                        end
                    end
                end
            end            
        end
        
        function map=IdByName(inMap)
            map=java.util.TreeMap;
            it=inMap.keySet.iterator;
            while it.hasNext
                key=char(it.next);
                if ~endsWith(key, '.color')
                    name=inMap.get(key);
                    if map.containsKey(key)
                        warning(['Duplicate use of ' name]);
                    else
                        map.put(name, key);
                    end
                end
            end
        end
        
        
    end
end
