%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

classdef FileBasics<handle
    methods(Static)
        function [ok,errMsg, errId]=MkDir(folder)
            errMsg='';
            errId=0;
            if exist(folder, 'dir')
                ok=true;
            else
                [ok, errMsg,errId]=mkdir(folder);
                if ~ok
                    disp(errMsg);
                end
            end
        end
        function [folder, file]=UiPut(dfltFolder, dfltFile, ttl)
            if nargin<3
                ttl='Save to which folder & file?';
            end
            FileBasics.MkDir(dfltFolder);
            [~,~,ext]=fileparts(dfltFile);
            done=false;
            if ismac
                jd=FileBasics.MsgAtTopScreen(ttl);
            else
                jd=[];
            end
            [file, folder]=uiputfile(['*' ext], ttl, ...
                fullfile(dfltFolder, dfltFile));
            if ~isempty(jd)
                jd.dispose;
            end
            if isempty(folder) || isnumeric(folder)
                folder=[];
                file=[];
            end
        end
        
        function out=UiGet(clue, folder, ttl)
            out=[];
            if ismac
                jd=FileBasics.MsgAtTopScreen(ttl);
            else
                jd=[];
            end
            [file, fldr]=uigetfile(clue,ttl, [folder '/']);
            if ~isempty(jd)
                jd.dispose;
            end
            if ~isnumeric(file) && ~isnumeric(fldr)
                out=fullfile(fldr,file);
            end
            if isempty(out)
                return;
            end
            
        end
        function jd=MsgAtTopScreen(ttl, pauseSecs)
            if nargin<2
                pauseSecs=8;
            end
            pu=showMsg(ttl, 'Note...', 'north east++', false, false, pauseSecs);
            jd=pu.dlg;
            jd.setAlwaysOnTop(true);
            fig=get(0, 'currentFigure');
            if ~isempty(fig)
                [~,~,~,~,pe]=Gui.FindJavaScreen(gcf);
                jd.setLocation(pe.x+(pe.width/2)-100, pe.y);
            end
        end
        
        function out=Canonize1(in, mMileIfContained, cur)
            if nargin<3
                cur=pwd;
                if nargin<2
                    mMileIfContained=false;
                end
            end
            ff=java.io.File(in);
            if ~ff.isAbsolute
                ff=java.io.File(fullfile(cur, in));
                out=char(ff.getCanonicalFile.toString);
            else
                out=char(ff.toString);
            end
            if mMileIfContained && exist(out, 'dir')
                f=FileBasics.FullPathList(out, true, true);
                if ~isempty(f)
                    out=f;
                end
            end
        end
        
        function out=Canonize(in, mFileIfContained, cur)
            if nargin<3
                cur=pwd;
                if nargin<2
                    mFileIfContained=false;
                end
            end
            if iscell(in)
                N=length(in);
                out=cell(1,N);
                for i=1:N
                    out{i}=FileBasics.Canonize1(in{i}, mFileIfContained, cur);
                end
            else
                out=FileBasics.Canonize1(in, mFileIfContained, cur);
            end
        end
        
        function l=FullPathList(filePath, mFilesOnly, firstFile)
            if nargin<3
                firstFile=false;
                if nargin<2
                    mFilesOnly=true;
                end
            end
            if iscell(filePath)
                l={};
                nPaths=length(filePath);
                for i=1:nPaths
                    l=[l FileBasics.FullPathList(filePath{i}, mFilesOnly, firstFile)];
                end
                return
            end
            if isempty(filePath)
                l={};
                return;
            end
            xxx=dir(filePath);
            N=length(xxx);
            l={};
            for i=1:N
                fl=xxx(i);
                [~,~,ext]=fileparts(fl.name);
                if ~mFilesOnly || (~fl.isdir && strcmpi(ext, '.m'))
                    if firstFile
                        l=fullfile(filePath, fl.name);
                        return;
                    end
                    l{end+1}=fullfile(filePath, fl.name);
                end
            end
        end
        
        %MatLab file search order when loading m file
        %   1 Current folder
        %   2 Already loaded in memory 
        %   3 Top of the path list down
        %Matlab search order for relative file refs with which or exist
        %   1 Current folder
        %   2 Top of the path list down
        function [ok, nConflicts, nSystemConflicts, nRemovals]...
                =AddNonConflictingPaths(addThese, forbidThese, ...
                talkToUserIfConflicts, checkAllFiles, ignoreSystemConflicts,...
                talkToUserIfNoConflicts)
            if nargin<6
                talkToUserIfNoConflicts=false;
                if nargin<5
                    ignoreSystemConflicts=true;
                    if nargin<4
                        checkAllFiles=false;
                        if nargout<3
                            talkToUserIfConflicts=true;
                            if nargin<2
                                forbidThese={};
                                if nargin<1
                                    addThese={mfilename('fullpath')};
                                end
                            end
                        end
                    end
                end
            end
            homeFldr=File.Home;
            originalAddThese=addThese;
            if ischar(addThese)
                addThese={addThese};
            end
            if ischar(forbidThese)
                forbidThese={forbidThese};
            end
            ok=true;
            try
                rootFldr=ctfroot;
            catch 
                rootFldr=[];
            end
            mapRemovable=TreeMapOfMany;
            mapNonRemovable=TreeMapOfMany;
            nToAdd=length(addThese);
            for ii=1:nToAdd
                addpath(addThese{ii});
            end
            app=BasicMap.Global;
            if ispc
                paths=strsplit(path, ';');
            else
                paths=strsplit(path, ':');
            end
            nPaths=length(paths);
            if checkAllFiles
                pu=PopUp('Exhaustively checking path conflicts');
                addThese=FileBasics.Canonize(addThese, false);            
                addThese=FileBasics.FullPathList(addThese, true, false);
                forbidThese=FileBasics.FullPathList(forbidThese, true, false);
                pu.initProgress(length(addThese)+length(forbidThese));
            else
                addThese=FileBasics.Canonize(addThese, true);
            end
            purify(addThese);
            purify2(forbidThese);
            nConflicts=mapRemovable.getUniqueValueCount;
            nSystemConflicts=mapNonRemovable.getUniqueValueCount;
            rerunCheckingAll=false;
            removals=[];
            if exist('pu', 'var')
                pu.close;
            end
            cancelled=false;
            if nConflicts>0 
                uniquePaths=mapRemovable.getUniqueValues;
                nUniquePaths=length(uniquePaths);
                if talkToUserIfConflicts
                    conflicted=mapRemovable.getKeyCountStrings(...
                        'conflicts with', 'path');
                    if nConflicts==1
                        h2=Html.H2(['There is ' ...
                        String.Pluralize2('path conflict', nConflicts) '.'], false);
                    else
                        h2=Html.H2(['There are ' ...
                            String.Pluralize2('path conflict', nConflicts) '.'], false);
                    end
                    if nSystemConflicts>0
                        conflicted=[conflicted ...
                            mapNonRemovable.getKeyCountStrings(...
                            'system conflicts with', 'path')];
                        strSysConflicts=[app.smallStart '(<b>and ' ...
                            String.Pluralize2('system path conflict', ...
                            nSystemConflicts) '</b>)<br><br>' app.smallEnd];
                    else
                        strSysConflicts='';
                    end

                    theMsg=Html.Wrap([h2 strSysConflicts ...
                        'This could lead to system instability<br>'...
                        Html.ToLimitedList(conflicted, 10) '<hr>']);
                    if ~checkAllFiles
                        choices={...
                            'Remove conflicting paths',...
                            'Keep conflicting paths',...
                            'Check every file in paths added'};
                        property='pathConflict1';
                    else
                        property='pathConflict2';
                        choices={'Remove conflicting paths',...
                            'Keep conflicting paths'};
                    end
                    [answer, cancelled]=Gui.Ask(theMsg, choices, ...
                        property, 'MatLab path conficts', 1);
                    if answer==1
                        if nSystemConflicts>0
                            disableIdxs=[nConflicts+1:nConflicts+nSystemConflicts];
                            up=[uniquePaths mapNonRemovable.getUniqueValues];
                        else
                            disableIdxs=[];
                            up=uniquePaths;
                        end
                        [removals, cancelled]=mnuMultiDlg(struct(...
                            'msg', 'Select conflicting non-system paths to remove', ...
                            'disableIdxs', disableIdxs), ...
                            'Removing paths', up, 0:nUniquePaths-1,false, true);
                    elseif answer==3
                        rerunCheckingAll=true;
                    else
                        ok=false;
                    end
                end
            elseif talkToUserIfNoConflicts
                if nSystemConflicts>0
                    strSysConflicts=[app.smallStart '(<b>but there IS ' ...
                        String.Pluralize2('system path conflict', ...
                        nSystemConflicts) '</b>)<br><br>' app.smallEnd];
                else
                    strSysConflicts='<br>';
                end
                theMsg=Html.WrapHr(['OK ... a quick scan suggests no path '...
                    'conflicts...<br>'... 
                    strSysConflicts '<b>Scan exhaustively every file in path added?</b>']);
                [rerunCheckingAll, cancelled]=askYesOrNo(...
                    theMsg, 'Adding path..', 'center', true, [], ...
                    'pathConflict4');
            end
            if cancelled
                ok=false;
            end
            if rerunCheckingAll
                [ok, nConflicts, nSystemConflicts, nRemovals]=...
                    FileBasics.AddNonConflictingPaths(...
                    originalAddThese, forbidThese, ...
                    talkToUserIfConflicts, true, false);
            else
                nRemovals=length(removals);    
                for ii=1:nRemovals
                    rIdx=removals(ii);
                    disp(['Removing ' uniquePaths(rIdx)]);
                    rmpath(uniquePaths{rIdx})
                end
            end
            
            function ok=isSystemPath(thePath)
                if ((~isempty(rootFldr) && ...
                        startsWith(thePath, rootFldr)) ...
                        || (~startsWith(thePath, homeFldr) && ...
                        contains(thePath, [filesep 'MATLAB'])))
                    ok=true;
                else
                    ok=false;
                end
            end
            
            function purify2(items)
                N=length(items);
                for i=1:N
                    item=items{i};
                    thePath=which(item);
                    if ~isempty(thePath)
                        if isSystemPath(thePath)
                            mapNonRemovable.put(item, thePath);
                        else
                            mapRemovable.put(item, thePath);
                        end
                    end
                    if exist('pu', 'var')
                        pu.incrementProgress;
                    end
                end
            end
            function purify(items)
                N=length(items);
                for i=1:N
                    item=items{i};
                    [fldr, fn, ext]=fileparts(item);
                    fileName=[fn ext];
                    [~,lastFldr, lastFldrExt]=fileparts(fldr);
                    lastFldr=[lastFldr lastFldrExt];
                    if isempty(fldr)
                        warning([ 'Ignoring no folder in "' item '"']);
                    else
                        for j=1:nPaths
                            thePath=paths{j};
                            file=fullfile(thePath, fileName);
                            if isequal(file, item)
                            elseif ignoreSystemConflicts ...
                                    && isSystemPath(thePath)
                                break;
                            elseif exist(file, 'file')
                                if isSystemPath(thePath)
                                    mapNonRemovable.put(fullfile(lastFldr,fileName), thePath);
                                else
                                    mapRemovable.put(fullfile(lastFldr,fileName), thePath);
                                end
                            end
                        end
                        if exist('pu', 'var')
                            pu.incrementProgress;
                        end
                    end
                end
            end
        end
    end
end