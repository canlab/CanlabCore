%
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

classdef Template < handle

    methods(Static)
        function [umap, badCnt, canLoad, reOrgData, reducedParams]=...
                Get(inData, d1, umapFile, maxSdu) 
            canLoad=true;
            if nargin<4
                maxSdu=1;
                if nargin<3
                    umapFile=[];
                end
            end
            reOrgData=[];
            reducedParams=[];
            umap=[];
            badCnt=0;
            dimNow=size(inData, 2);
            while true
                umap=[];
                if isempty(umapFile)
                    umapFile=getUmapFile;
                end
                if ~isempty(umapFile)
                    sdus=[];
                    dimUmap=0;
                    try
                        try
                            load(umapFile, 'umap');
                        catch
                            canLoad=false;
                            break;
                        end
                        umapFile=[];
                        try
                            if ~isempty(umap.supervisors)
                                umap.supervisors.prepareForTemplate;
                            end                            
                            d2=umap.dimNames;
                            if ~isempty(d2) && ...
                                    ~StringArray.AreSameOrEmpty(d2, d1)
                                [reOrgData, reducedParams]=MatBasics.ReOrg(...
                                    inData, d1, d2, true, false);
                                isOk=~isempty(reOrgData);
                                if isOk && ~isempty(reducedParams)
                                    html=Html.To2Lists(d2,d1,'ol', ...
                                        'Template', 'Current', true, 25);
                                    answ=ask(['<html><b>'...
                                        'Parameters are a subset of each other</b><hr>'...
                                        html '<br><br><center>'...
                                        'Accept reduced parameters?'...
                                        '</center></html>'], 'Problem...',...
                                        'error');
                                    if answ==-1
                                        umap=[];
                                        return;
                                    elseif answ==0
                                        continue;
                                    end
                                end
                                if ~isOk
                                    html=Html.To2Lists(d2,d1,'ol', ...
                                        'Template', 'Current', true, 25);
                                    showMsg(['<html><font color="red"><b>'...
                                        'Parameters differ</b><hr>'...
                                        html '</html>'], 'Problem...', 'error');
                                    badCnt=badCnt+1;
                                    continue;
                                else
                                    inData=reOrgData;
                                end
                            end
                        catch ex
                            disp(ex);
                        end
                        if isprop(umap, 'rawMeans') && ~isempty(umap.rawMeans)
                            sdus=MatBasics.SduDist2(inData, umap.rawMeans, umap.rawStds);
                        else
                            sdus=MatBasics.SduDist(inData, umap.raw_data);
                        end
                        dimUmap=size(umap.raw_data, 2);
                    catch ex
                        ex.getReport
                    end
                    if isempty(sdus)
                        s=sprintf(['Chosen template has %d data dimensions'...
                            '<br>and current # of dimensions is %d'], ...
                            dimUmap, dimNow);
                        showMsg(Html.WrapHr(s),'Incompatible...', 'error');
                    elseif any(sdus>maxSdu)
                        badCnt=badCnt+1;
                        showMsg(Html.WrapHr([num2str(sum(sdus>maxSdu))...
                            ' standard deviation unit(s) ' ...
                            '<br>are greater than ' num2str(maxSdu) ...
                            '<br><b>' MatBasics.toRoundedTable(sdus, 2, ...
                            find(sdus>maxSdu)) '</b>']), ...
                            'Incompatible...', 'error');
                    else
                        if isempty(d2)
                            warning(['This template has no dimension names'...
                                ' to aid data compatibility checking!!']);
                            msg(Html.WrapHr(...
                                ['This template has no '...
                                'dimension names<br>to aid data '...
                                'compatibility checking!!<br><br>'...
                                'Good luck....']), 8, 'north east+', ...
                                'Template is vague');
                        end
                        break;
                    end
                else
                    break;
                end
            end
            
            function umapFile=getUmapFile()
                umapFile=FileBasics.UiGet('*.umap.mat', pwd, ...
                    'Select prior compatible UMAP template');
            end
            
        end
        
        function ok=Save(umap, inputFile)
            umapFile=getNewUmapFile(inputFile);
            ok=~isempty(umapFile);
            if ok
                umap.progress_callback=[];
                umap.graph=[];
                pu=PopUp('Saving template');
                save(umapFile, 'umap');
                pu.close;
            end
            
            
            function umapFile=getNewUmapFile(file)
                [fldr, fl, ~]=fileparts(file);
                [fldr, file]=FileBasics.UiPut(fldr, [fl '.umap.mat'], ...
                    'Save UMAP as guiding template');
                if isempty(fldr)
                    umapFile=[];
                else
                    if ~String.EndsWith(file, '.umap.mat')
                        file=[file(1:end-4) '.umap.mat'];
                        if exist(fullfile(fldr,file), 'file')
                            answer=questdlg(['Template "' ...
                                file '" already '...
                                'exists ... Replace?']);
                            if isempty(answer) || isequal('Cancel', answer)
                                umapFile=[];
                                return;
                            elseif ~yes
                                umapFile=this.getNewUmapFile;
                                return;
                            end
                        end
                    end
                    umapFile=fullfile(fldr, file);
                end
            end
        end
        
        function [percNewSubsets, newSubsetIdxs, newSubsetCnt]=...
                CheckForUntrainedFalsePositives(template, inData, ...
                sduLimit, parameterLimit)
            if isempty(template.supervisors)
                percNewSubsets=0;
                newSubsetIdxs=[];
                newSubsetCnt=0;
                return;
            end
            if nargin<3
                sduLimit=3.66;
                parameterLimit=2;
            end
            [newSubsetCnt, newSubsetIdxs]=detectUnsupervised(template, inData, ...
                sduLimit, parameterLimit);
            R=size(inData,1);
            percNewSubsets=newSubsetCnt/R*100;
        end
        
        function [choice, cancelled]=Ask(perc)
            html=Html.WrapHr([String.encodeRounded(perc, 1) '% of the '...
                'rows/events appear NEW <br> or unseen in the data '...
                'that trained<br> this supervised template<br>']);
            choices={'Use anyway (it''s okay)', ...
                'Re-supervise via SDU distance'};
            [choice, cancelled]=Gui.Ask(html, choices, ...
                'umapNew2Template', 'New subsets?', 1);
        end
       
    end
end