%   Class for probability binning described at
%   https://www.ncbi.nlm.nih.gov/pubmed/11598946
%
%   Bioinformatics inventors
%   Original:  Roederer M1, Moore W, Treister A, Hardy RR, Herzenberg LA.
%   Incorporation into QF matching:  Darya Orlova <dyorlova@gmail.com>
%
%   Software Developers: Connor Meehan <connor.gw.meehan@gmail.com>
%           Stephen Meehan <swmeehan@stanford.edu> 
%
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef File
    properties(SetAccess=private)
        file='';
        absolute='';
        firstFolder='';
    end
    
    properties(Constant)
        DEBUG=0;
        EXPORT_TYPES={'*.txt', 'tab-delimited file';...
            '*.csv', 'comma separated values';...
            '*.xls', 'Excel workbook'};
        PROP_EXPORT='exportType2';
        PROP_XLS_IMG='xlsGateImg';
        PROP_XLS_ADD='lastXls';
        PROP_XLS_OPEN='isOpenSlidesWorkbook';
    end
    
    methods        
        function [this]= File(file)
            this.file=file;
            this.absolute=File.canonicalPath(file);
            if ~isdir(this.absolute)
                [this.firstFolder, ~, ~]=fileparts(this.absolute);
            else
                this.firstFolder=this.absolute;
            end
        end
        
        function [levels]=countAncestors(this, ancestorFolder)
            ancestorFolder=File.canonicalPath(ancestorFolder);
            N=length(ancestorFolder);
            folder=this.firstFolder;
            levels=-1;
            while strncmp(ancestorFolder, folder, N)
                levels=levels+1;
                [folder, ~, ~]=fileparts(folder);
            end
        end
    end
    
    methods(Static)
        function [size, map]=Tally(folder, spec, parent, map, cur)
            if nargin<3
                parent=[];
            end
            if nargin<4
                map=java.util.TreeMap;
            end
            size=0;
            fsp=fullfile(folder,spec);
            xxx=dir(fsp);
            N=length(xxx);
            if nargin<5
                cur=now;
            end
            if isempty(parent)
                scan=true;
            else
                idx=String.LastIndexOf(folder, filesep);
                if idx>0
                    lastPath=String.SubString(folder, idx+1);
                    scan=strcmp(lastPath, parent);
                else
                    scan=false;
                end
            end
            if scan
                for i=1:N
                    fl=xxx(i);
                    fname=fl.name;
                    if ~fl.isdir
                        subSize=fl.bytes;
                        elapsed=cur-fl.datenum;
                        path=fullfile(folder, fname);
                        map.put(elapsed, path);
                        size=size+subSize;
                    end
                end
            end
            if ~isempty(parent)
                path=fullfile(folder, parent);
                if exist(path, 'dir')
                    subSize=File.Tally(path, spec, parent, map, cur);
                    size=size+subSize;
                end
            end
            fsp=folder;
            xxx=dir(fsp);
            N=length(xxx);
            for i=1:N
                fl=xxx(i);
                fname=fl.name;
                if fl.isdir
                    if ~strcmp('.', fname) && ~strcmp('..', fname)
                        if isempty(parent) || ~strcmp(fname, parent)
                            path=fullfile(folder, fname);
                            subSize=File.Tally(path, spec, parent, map, cur);
                            size=size+subSize;
                        end
                    end
                end
            end
        end
        
        function [sz, str]=Size2(to, matches)
            sz=0;
            existType=exist(to);
            if existType>0
                files=dir(to);
                if existType==2 %file
                    sz=files(1).bytes;
                elseif existType==7 %directory
                    N=length(files);
                    for i=1:N
                        if files(i).isdir
                            n=files(i).name;
                            if ~isequal(n, '.') && ~isequal(n, '..')
                                sz=sz+File.Size2(fullfile(to, ...
                                    files(i).name), matches);
                            end
                        else
                            sz=sz+files(i).bytes;
                        end
                    end
                end
                File.SizeMatches(to, matches);
            end
            if nargout>1
                str='';%['<td align="right">' String.encodeBytes(sz) '</td>'];
                it=matches.keySet.iterator;
                while it.hasNext
                    sz2=matches.get(it.next);
                    str=[str '<td align="right">' ...
                        String.encodeBytes(sz2) '</td>'];
                end
            end
        end
        
        function SizeMatches(dirName, matches)
            it=matches.keySet.iterator;
            while it.hasNext
                match=char(it.next);
                sp=fullfile(dirName, match);
                dd=dir(sp);
                N=length(dd);
                sz=0;
                for i=1:N
                    sz=sz+dd(i).bytes;
                end
                matches.put(match, matches.get(match)+sz);
            end
        end
        
        function [sz, str]=Size(to)
            sz=0;
            existType=exist(to);
            if existType>0
                files=dir(to);
                if existType==2 %file
                    sz=files(1).bytes;
                elseif existType==7 %directory
                    N=length(files);
                    for i=1:N
                        if files(i).isdir
                            n=files(i).name;
                            if ~isequal(n, '.') && ~isequal(n, '..')
                                sz=sz+File.Size(fullfile(to,files(i).name));
                            end
                        else
                            sz=sz+files(i).bytes;
                        end
                    end
                end
            end
            if nargout>1
                str=String.encodeBytes(sz);
            end
        end
        
        function CopyBySpec(from, someSpec, to)
            spec=fullfile(from, someSpec);
            File.mkDir(to);
            copyfile(spec, to);
        end
        
        function ok=Copy(from, to, retrying)
            ok=true;
            try
                copyfile(from, to);
            catch ex
                ok=false;
                if ~exist(fileparts(to), 'dir')
                    File.mkDir(fileparts(to))
                    ok=File.Copy(from, to, true);
                    return;
                end
                if nargin>2 && retrying
                    rethrow(ex);
                else
                    ok=File.Copy1by1(from, to);
                end
            end
        end
        function ok=Copy1by1(from, to)
            ok=true;
            a=dir(from);
            N=length(a);
            for i=1:N
                n=a(i).name;
                if isequal(n, '.') || isequal(n, '..')
                else
                    ff=fullfile(a(i).folder, a(i).name);
                    if ~File.Copy(ff, to, true)
                        ok=false;
                        return;
                    end
                end
            end
        end

        function ok=isFile(input)
            ok=false;
            if exist(input, 'file')
                ok=~isdir(input);
            end
        end
        
        function file=AppendSuffix(file, suffix)
            [fldr, ff, ext]=fileparts(file);
            file=fullfile(fldr, [ff suffix ext]);
        end
        
        function fullFile=SwitchExtension(fullFile, type2)
            lSepIdx1=String.LastIndexOf(type2, filesep);
            lSepIdx2=String.LastIndexOf(fullFile, filesep);
            lExtIdx1=String.LastIndexOf(type2, '.');
            lExtIdx2=String.LastIndexOf(fullFile, '.');
            if lExtIdx1>0 && lExtIdx1 > lSepIdx1 && lExtIdx2>0 ...
                    && lExtIdx2>lSepIdx2
                ext1=String.SubString(type2, lExtIdx1);
                ext2=String.SubString(fullFile, lExtIdx2);
                if ~strcmp(ext1,ext2)
                    fullFile=[fullFile(1:lExtIdx2-1) ext1];
                end
            end
        end
        
        function fullFile=SwitchExtension2(fullFile, newExt)
            if endsWith(fullFile, newExt)
                return;
            end
            sepIdx=String.LastIndexOf(fullFile, filesep);
            idx=String.LastIndexOf(fullFile, '.');
            if idx>0 && idx>sepIdx
                fullFile=[fullFile(1:idx-1) newExt];
            end
        end
        
        function floats=ReadFloat32(fileName)
            f=[];
            floats=[];
            try
                f = fopen(fileName,'r');
                sz=fread(f, [1 2], 'int32', 'b');
                floats=fread(f, sz, 'float32', 'b');
                fclose(f);
            catch ex
                ex.getReport
                if ~isempty(f)
                    fclose(f);
                end
            end
        end
        
        function WriteFloat32(fileName, floats)
            f=[];
            try
                f = fopen(fileName,'w');
                if f==-1
                    fileattrib(fileName, '+w');
                    f = fopen(fileName,'w');
                    if f<0
                        return;
                    end
                end
                fwrite(f, size(floats), 'int32', 'b');
                fwrite(f, floats, 'float32', 'b');
                fclose(f);
            catch ex
                ex.getReport
                if ~isempty(f)
                    fclose(f);
                end
            end
        end
        
        function fn=SaveTempHtml(html)
            fn=File.SaveTempTextFile(html, '.html');
        end
        
        function fn=SaveTempTextFile(text, fileExtension)
            if nargin<2
                fileExtension='.txt';
            end
            fn=[tempdir num2str(floor(now)) fileExtension];
            File.SaveTextFile(fn, text)
        end
        
        function ok=SaveMatrix(fileName, matrix, useLabels)
            if nargin<3
                useLabels=true;
            end
            ok=false;
            try
                if useLabels
                    f = fopen(fileName,'w');
                    if f==-1
                        fileattrib(fileName, '+w');
                        f = fopen(fileName,'w');
                        if f<0
                            return;
                        end
                    end
                    [~,C]=size(matrix);
                    for c=1:C
                        if c>1
                            fprintf(f, ', ');
                        end
                        fprintf(f, 'column%d', c);
                    end
                    fprintf(f, '\n');
                    fclose(f);
                    if isinteger(matrix)
                        dlmwrite(fileName, matrix, '-append',  'precision', '%13d');
                    else
                        dlmwrite(fileName, matrix, '-append', 'precision', 13);
                    end
                else
                    csvwrite(fileName, matrix);
                end
                ok=true;
            catch ex
                ok=false;
                ex.getReport
            end
        end

        function [matrix, ok]=ReadMatrix(fileName)
            try
                matrix=csvread(fileName, 1);
                ok=true;
            catch ex
                matrix=[];
                ok=false;
                ex.getReport
            end
        end

        function props=ReadProperties(fileName, getEmpty)
            props=java.util.Properties;
            try
                fis=java.io.FileInputStream(strrep(fileName, '~', File.Home));
                props.load(fis);
                fis.close;
            catch ex
                fldr=fileparts(fileName);
                if isempty(fldr)
                    onMatLabPath=which(fileName);
                    if ~isempty(onMatLabPath)
                        warning(['Found properties on MatLab path?:  "' onMatLabPath '"']);
                        props=File.ReadProperties(onMatLabPath);
                        return;
                    end
                end
                if nargin<2 || ~getEmpty
                    ex.getReport
                    props=[];
                end
            end    
        end
        
        function props=GetProperties(fileName)
            props=File.ReadProperties(fileName, true);
        end
        
        function ok=SaveProperties2(fileName, props, headerString)
            ok=false;
            try
                fos=java.io.FileOutputStream(fileName);
                if nargin<3
                    headerString='';
                end
                props.save(fos, headerString);
                fos.close;
                ok=true;
            catch ex
                ex.getReport
            end
        end
        
        function SaveTextFile(fileName, text)
            f=[];
            try
                f = fopen(fileName,'w');
                if f==-1
                    fileattrib(fileName, '+w');
                    f = fopen(fileName,'w');
                    if f<0
                        return;
                    end
                end
                fwrite(f, text);
                fclose(f);
            catch ex
                ex.getReport
                if ~isempty(f)
                    fclose(f);
                end
            end
        end
        
        function fid=Save(nameOrFileHandle, joStyle, ci, ci2, labels,...
                matrix, close)
            if ischar(nameOrFileHandle)
                fid=fopen(nameOrFileHandle, 'w');
            else
                fid=nameOrFileHandle;
            end
            if joStyle
                eol='\r';
            else
                eol='\n';
            end
            N2=length(ci);
            if ~joStyle
                fprintf(fid, '\t');
            end
            for i = 1:N2
                fprintf(fid, '%s', labels{ci(i)});
                if i<N2
                    fprintf(fid, '\t');
                end
            end
            fprintf(fid, eol);
            for i = 1:N2
                if ~joStyle
                    fprintf(fid, '%s\t', labels{ci(i)});
                end
                for j = 1:N2
                    fprintf(fid, '%1.8f', matrix(ci2(i), ci2(j)));
                    if j<N2
                        fprintf(fid, '\t');
                    end
                end
                fprintf(fid, eol); 
            end
            fprintf(fid, eol);
            if nargin>6 && close
                fclose(fid);
            end
        end
        
        function ConfirmParentFolderExists(file)
            [folder, ~]=fileparts(file);
            File.mkDir(folder);
        end
        
        function [ok,errMsg, errId]=moveFile(from, to, force)
            ok=true;
            errMsg='';
            errId=0;            
            try
                removeFrom=false;
                if ispc
                    if ~contains(from, '*') &&...
                            ~contains(to, '*')
                        if exist(from, 'dir')
                            if ~exist(to, 'dir') && ~exist(to, 'file')
                                fromDir=from;
                                from=fullfile(from, '*.*');
                                File.mkDir(to);
                                removeFrom=true;
                                pause(.25);
                                exist(to)
                                disp('Avoiding windows Copy question yes/skip/cancel');
                            end
                        end
                    end
                end
                if nargin>2 && force
                    [ok,errMsg, errId]=movefile(from, to, 'f');                    
                else
                    [ok,errMsg, errId]=movefile(from, to);
                end
                if removeFrom
                    if exist(to, 'dir')
                        File.rmDir(fromDir);
                    end
                end
                if ~ok
                    disp(errMsg);
                end
            catch ex
                disp(ex);
            end
        end
        
        function [ok,errMsg, errId]=rmDir(to, killIfFile)
            ok=true;
            errMsg='';
            errId=0;
            if exist(to, 'dir')
                try
                    [ok,errMsg, errId]=rmdir(to, 's');
                catch ex
                    if exist([to sprintf('/Icon\r')], 'file')
                        delete([to sprintf('/Icon\r')]);
                        rmdir(to, 's');
                    end
                end
            elseif nargin>1 && killIfFile && exist(to, 'file')
                delete(to);
            elseif exist(to, 'file')
                [ok,errMsg, errId]=rmdir(to, 's');
            end
            if ~ok
                disp(errMsg);
            end
        end

        function [ok,errMsg, errId]=emptyDir(to, killIfFile)
            if nargin<2
                killIfFile=false;
            end
            File.rmDir(to, killIfFile);
            [ok,errMsg, errId]=File.mkDir(to);
        end


        function [ok,errMsg, errId]=mkDir(folder)
            errMsg='';
            errId=0;
            if isempty(folder)
                ok=true;
                return;
            end
            if exist(folder, 'dir')
                ok=true;
            else
                [ok, errMsg,errId]=mkdir(folder);
                if ~ok
                    disp(errMsg);
                end
            end
        end
    
        function path=TrimUserRoot(path)
            ur=File.Home;
            if ~isempty(path) && String.StartsWith(path, ur)
                N=length(ur);
                path=['HOME' String.SubString2(path, N+1, length(path)+1)];
            end
        end
                    
        function ok=isSameOrAncestor(ancestor, descendant)
            file_=File(descendant);
            ok=file_.countAncestors(ancestor)>=0;
        end
        
        function ok=endsWith( path, name )
            path=String.PruneSuffix(path, filesep);
            name=String.PruneSuffix(name, filesep);
            [p, f, ~]=fileparts(path);
            ok=strcmp(f,name);
        end
        
        function folder=parent(path)
            [folder,~,~]=fileparts(File.canonicalPath( path ) );
        end
        
        function  path=canonicalPath( path )
            path=File.absolutePath( path, false);
        end       
        
        function yes=isEmpty(filespec)
            files=dir(filespec);
            yes=false;
            if isempty(files)
                yes=true;
            else
                if length(files)==2
                    if isdir(filespec) || String.EndsWith(filespec,'*.*')
                        if strcmp(files(1).name,'.')
                            if strcmp(files(2).name, '..')
                                yes=true;
                            end
                        end
                    end
                end
            end                 
        end
        
        function ok=IsAbsolute(path)            
            f = java.io.File(path);
            ok=f.isAbsolute;
        end
        
        function [relativePath, parentFolderCount]=...
                GetSubPath(descendent, ancestor)
            relativePath=descendent;
            d=java.io.File(descendent);
            if d.isAbsolute
                p={};
                a=java.io.File(ancestor);
                while ~isempty(d)
                    if d.equals(a)
                        parentFolderCount=length(p);
                        if parentFolderCount==0
                            relativePath='';
                        else
                            for i=parentFolderCount:-1:1
                                if i==parentFolderCount
                                    relativePath=p{i};
                                else
                                    relativePath=[relativePath ...
                                        filesep p{i}];
                                end
                            end
                        end
                        return;
                    end
                    p{end+1}=char(d.getName);
                    d=d.getParentFile;
                end
            end
            parentFolderCount=-1;
        end
        
        function [path, added]=AddParentPath(fromPath, toPath)
            added=~File.IsAbsolute(fromPath);
            if ~added
                path=fromPath;
            else
                path=fullfile(toPath, fromPath);
            end
        end
        
        
        function abs_path=absolutePath( path,  throwErrorIfFileNotExist )            
            % 2nd parameter is optional:
            if nargin < 2
                throwErrorIfFileNotExist = false;
            end
            
            %build absolute path
            file = java.io.File(path);
            if ~file.isAbsolute()
                file=java.io.File(cd, path);
            end
            abs_path = char(file.getCanonicalPath());
            
            %check that file exists
            if throwErrorIfFileNotExist && ~exist(abs_path, 'file')
                throw(MException('absolutePath:fileNotExist', 'The file %s doesn''t exist', abs_path));
            end
        end
    
        function path=addSeparatorIfNecessary(path)
            if ~String.EndsWith(path, filesep)
                path=strcat(path,filesep);
            end
        end
        
        function [ok, parentFolder]=IsLastFolder(...
                path, lastFolder)
            ok=false;            
            f=java.io.File(path);
            parentFolder=char(f.getParent().toString());
            if f.isDirectory()
                if strcmp(f.getName(), lastFolder)
                    ok=true;
                end
            end
        end
        
        function name=Name(path)
            [~, name, ext]=fileparts(path);
            name=[name ext];
        end
        
        function folder=GetRoot(files, homePrefix)
            if nargin==1
                homePrefix=false;
            end
            N=length(files);
            if iscell(files) && N>0
                folder=File.getFolder(files{1});
                for i=2:N
                    if ~isempty(files{i})
                        folder=File.getCommon(folder, files{i});
                    end
                end
            else
                folder=File.getFolder(files);
            end
            folder=File.addSeparatorIfNecessary(folder);
            if homePrefix
                folder=File.TrimUserRoot(folder);
            end
        end
        
        function [folder]=getCommon(f1,f2)
            folder='/';
            while 1
                folder1=File.getFolder(f1);
                folder2=File.getFolder(f2);
                N1=length(folder1);
                N2=length(folder2);
                if N1<=N2
                    if strncmp(folder1, folder2, N1)
                        folder=folder1;
                        break;
                    end
                else
                    if strncmp(folder2, folder1, N2)
                        folder=folder2;
                        break;
                    end
                end
                [f1, ~, ~]=fileparts(folder1);
                [f2, ~, ~]=fileparts(folder2);
            end
        end
        
        function [folder]=getFolder(f)
            if isdir(f)
                folder=f;
            else
                [folder,~,~]=fileparts(f);
            end
            folder=File.canonicalPath(folder);
        end
        
        function [commonParentFolder, list]=Abbreviate(files)
            list={};
            if ~isempty(files)
                commonParentFolder=File.GetRoot(files);
                prefix=length(commonParentFolder)-1;
                N=length(files);
                for i=1:N
                    file=File.canonicalPath(files{i});
                    if ~isempty(file)
                        [path, name, ext]=fileparts(file);
                        N2=length(path);
                        if N2 < prefix
                            list{i}='.';
                        elseif N2 > prefix
                            subPath=path(1, [prefix+2:N2]);
                            list{i}=[name ext ' (in ' subPath ')'];
                        else
                            list{i}=[name ext];
                        end
                    end
                end
            else
                commonParentFolder='';
            end
        end
        
        
        function ok=ExistsOrOk(file, type)
            if ~isempty(file)
                if nargin == 1
                    type='This file';
                end
                if ~exist(file)
                    no='Are you kidding? .. no, no, NO!';
                    yes='Sure .. (why not?)';
                    [folder,fileName, ext]=File.Parts(file);
                    ok=strcmp('Ok', questdlg([ '"' fileName ext ...
                        '" does not exist in "' folder '"'...
                        '.  Do you wish to continue?'], [type ' is not found '],...
                        yes, no, no));
                else
                    ok=true;
                end
            else
                ok=false;
            end
        end
        
         function [ newPath, name, ext ] = Parts( path )
             %UNTITLED Summary of this function goes here
             %   Detailed explanation goes here
             [newPath, name, ext]=fileparts(path);
             if isempty(newPath)
                 newPath = ['.' filesep];
             else
                 newPath = strcat(newPath, filesep);
             end
         end
         
         function [list1, list2] = AddToLists(folder, file, list1, list2)
             item=[folder file ];
             idx=getnameidx(list1, item);
             if (idx==0)
                 list1{end+1}=item;
                 list2{end+1}=[file ' (' folder ')' ];
             end
         end
         
         function [set] = AddToSet(folder, file, set)
             item=[folder file ];
             idx=getnameidx(set, item);
             if (idx==0)
                 set{end+1}=item;
             end
         end
         function files=ToList(pathString, canonical)
             if strcmp('cell', class(pathString))
                 files=pathString;
                 return;
             end
             l=dir(pathString);
             if nargin==1
                 canonical=false;
             end
             if canonical
                 pathString=File.canonicalPath(pathString);
             end
             [p, ~, ~]=File.Parts(pathString);
             n=length(l);
             files=cell(1, n);
             for i=1:n
                 f=[p l(i).name];
                 files{i}=f;
             end
         end

         
        function [home]=Home
            home=char( java.lang.System.getProperty('user.home') );
        end

        function [ok, file]=SwitchRoot(file, sub)
            str=fullfile(File.Home, sub);
            if String.StartsWith(file, str)
                idx=length(str);
                file=fullfile(sub, String.SubString(file, idx+1));
                ok=true;
            else
                ok=false;                
            end
        end
        
        function abbreviation=AbbreviateRoot(file)
            h=File.Home;
            [ok, abbreviation]=File.SwitchRoot(file, 'Desktop');
            if ~ok
                [ok, abbreviation]=File.SwitchRoot(file, 'Documents');
                if ~ok
                    [ok, abbreviation]=File.SwitchRoot(file, 'Downloads');
                    if ~ok
                        home=char( java.lang.System.getProperty('user.name') );
                        [ok, abbreviation]=File.SwitchRoot(file, home);
                        if ~ok
                            abbreviation=file;
                        end
                    end
                end
                
            end
        end
        
        function description=Describe(file)
            [folder,name,ext]=fileparts(file);
            [folder, parentName, parentExt]=fileparts(folder);
            folder=File.AbbreviateRoot(folder);
            description=[parentName parentExt filesep name ext ' (in ' folder ')'];
        end
        
        function files=FindAll(path, fileSpec)
            l=dir(path);
            N=length(l);
            files={};
            for i=1:N
                d=l(i);
                if d.isdir && ~String.StartsWith(d.name, '.')
                    subFiles=File.FindAll(fullfile(path, d.name), fileSpec);
                    if ~isempty(subFiles)
                        if isempty(files)
                            files=subFiles;
                        else
                            files=[files(1,:) subFiles(1,:)];
                        end
                    end
                end
            end
            l=dir(fullfile(path,fileSpec));
            N=length(l);
            subPath=fileparts(fileSpec);
            for i=1:N
                f=fullfile(path, subPath, l(i).name);
                if ~l(i).isdir
                    files{end+1}=f;
                end
            end
        end
        function txt=ReadTextFile(file)
            txt=char(edu.stanford.facs.swing.CpuInfo.readTextFile(file));
        end
        function yes=AreEqual(fileName1, fileName2)
            file_1 = javaObject('java.io.File', fileName1);
            file_2 = javaObject('java.io.File', fileName2);
            yes=javaMethod('contentEquals','org.apache.commons.io.FileUtils',...
                file_1, file_2);
        end
        function yes=TextFilesAreEqual(file1, file2)
            yes=false;
            fid1 = fopen(file1, 'r');
            fid2 = fopen(file2, 'r');
            if fid1 ~= -1 && fid2 ~= -1
                lines1 = textscan(fid1,'%s','delimiter','\n');
                lines2 = textscan(fid2,'%s','delimiter','\n');
                lines1 = lines1{1};
                lines2=lines2{1};
                yes = isequal(lines1,lines2);
            end
            
            if fid1~=-1
                fclose(fid1);
            end
            if fid2~= -1
                fclose(fid2);
            end
        end
        
        function ok=WriteTextFile(file, lines)
            if isnumeric(lines)
                lines=MatBasics.ToStrs(lines);
            end
            try
                fid=fopen(file, 'wt');
                N=length(lines);
                if isa(lines, 'java.lang.Object[]')
                    for i=1:N
                        fprintf(fid, '%s\n', lines(i));
                    end
                elseif iscell(lines) 
                    for i=1:N
                        fprintf(fid, '%s\n', lines{i});
                    end
                else
                    fprintf(fid, '%s\n', char(lines));
                end
                fclose(fid);
                ok=true;
            catch ex
                disp(ex.message);
                ok=false;
            end
        end
        
        function file=FindFirst(path, fileSpec)
            file=[];
            l=dir(fullfile(path,fileSpec));
            N=length(l);
            for i=1:N
                f=fullfile(path, fileparts(fileSpec), l(i).name);
                if ~l(i).isdir
                    file=f;
                    return;
                end
            end
            l=dir(path);
            isub = [l(:).isdir]; %# returns logical vector
            nameFolds = {l(isub).name}';
            N=length(nameFolds);
            for i=1:N
                d=nameFolds{i};
                if d(1) ~= '.'
                    subFile=File.FindFirst(fullfile(path, d), fileSpec);
                    if ~isempty(subFile)
                        file=subFile;
                        return;
                    end
                end
            end
        end
        
        function SaveProperties(p,f, sortKeys)
            l={};
            if sortKeys
                ts=java.util.TreeSet;
                ts.addAll(p.keySet);
                ks=ts.iterator;
            else
                ks=p.keySet.iterator;
            end
            while ks.hasNext
                k=char(ks.next);
                v=char(p.get(k));
                l{end+1}=[k '=' v];
            end
            File.WriteTextFile(f, l);
        end
        
        function out=Html(path)
            out='';
            if isempty(path)
                return;
            end
            cnt=0;
            home=File.Home;
            prefix='';
            while true
                [path, name, ext]=fileparts(path);
                cnt=cnt+1;
                if cnt==1
                    out=[name ext out];
                elseif cnt==2
                    out=[name ext filesep out];
                end
                if isempty(path) || length(path)<=3
                    break;
                elseif isequal([home filesep 'Documents'], path)
                    prefix=['<b>~' filesep 'Documents</b>'];
                    break;
                elseif isequal([home filesep 'Desktop'], path)
                    prefix=['<b>~' filesep 'Desktop</b>'];
                    break;
                elseif isequal(home, path)
                    prefix=['<b>~</b>'];
                    break;
                end
            end
            if cnt-2>0
                out=[prefix filesep '<i><font color=''blue''>' ...
                    String.Pluralize2('folder', ...
                    cnt-2) '</font></i>' filesep out ];
            else
                out=[prefix filesep out];
            end
        end
        
        function [bytes, map]=DirRecursive(folder, spec, parent, map, cur, pu)
            if nargin<6
                pu=[];
                if nargin<5
                    cur=[];
                    if nargin<4
                        %map=java.util.TreeMap;
                        map=[];
                        if nargin<3
                            parent=[];
                        end
                    end
                end
            end
            if isempty(map)
                %map=java.util.TreeMap;
                map=TreeMapOfMany;
            end
            bytes=0;
            fsp=fullfile(folder,spec);
            diskItems=dir(fsp);
            N=length(diskItems);
            if isempty(cur)
                cur=now;
            end
            if isempty(parent)
                scan=true;
            else
                idx=String.LastIndexOf(folder, filesep);
                if idx>0
                    lastPath=String.SubString(folder, idx+1);
                    scan=strcmp(lastPath, parent);
                else
                    scan=false;
                end
            end
            if scan
                if ~isempty(pu)
                    pu.setText2([String.Pluralize2('file', N) ...
                        ' in ' File.Html(folder) ]);
                    if File.DEBUG>0
                        disp([String.Pluralize2('file', N) ...
                        ' in ' File.Html(folder) ])
                    end
                end
                for i=1:N
                    fl=diskItems(i);
                    fname=fl.name;
                    if ~fl.isdir
                        elapsed=cur-fl.datenum;
                        path=fullfile(folder, fname);
                        if File.DEBUG>1 && map.containsKey(elapsed)
                            fprintf(['%s size=%s, elapsed=%s, '...
                                'prior=%s\n'], fname, ...
                                String.encodeInteger(fl.bytes), ...
                                String.encodeRounded(elapsed, 3), ...
                                map.get(elapsed));
                        end
                        map.put(elapsed, path);
                        bytes=bytes+fl.bytes;                        
                    end
                end
                if ~isempty(pu) && pu.cancelled
                    return;
                end
            end
            if ~isempty(parent)
                path=fullfile(folder, parent);
                if exist(path, 'dir')
                    bytes=bytes+...
                        File.DirRecursive(path, spec, parent, map, cur, pu);
                end
            end
            fsp=folder;
            diskItems=dir(fsp);
            N=length(diskItems);
            for i=1:N
                fl=diskItems(i);
                fname=fl.name;
                if fl.isdir
                    if ~strcmp('.', fname) && ~strcmp('..', fname)
                        if isempty(parent) || ~strcmp(fname, parent)
                            path=fullfile(folder, fname);
                            bytes_=File.DirRecursive(path, spec, parent, map, cur, pu);
                            bytes=bytes+bytes_;
                        end
                    end
                    if ~isempty(pu) && pu.cancelled
                        return;
                    end
                end
            end            
        end
        
        function [bytes, items]=DiskUsage(folder, spec, olderThanDays, pu, cur)
            if nargin<5
                cur=now;
                if nargin<4
                    pu=[];
                    if nargin<3
                        olderThanDays=0;
                    end
                end
            end
            bytes=0;
            items=0;
            fsp=fullfile(folder,spec);
            diskItems=dir(fsp);
            N=length(diskItems);
            if ~isempty(pu)
                pu.setText2([String.Pluralize2('file', N) ...
                    ' in ' File.Html(folder) ]);
                if File.DEBUG>0
                    disp([String.Pluralize2('file', N) ...
                        ' in ' File.Html(folder) ])
                end
            end
            for i=1:N
                fl=diskItems(i);
                if ~fl.isdir
                    if olderThanDays==0
                        bytes=bytes+fl.bytes;
                        items=items+1;
                    else
                        if cur-fl.datenum>olderThanDays
                            bytes=bytes+fl.bytes;
                            items=items+1;
                        end
                    end
                end
            end
            if ~isempty(pu) && pu.cancelled
                return;
            end            
            fsp=folder;
            diskItems=dir(fsp);
            N=length(diskItems);
            for i=1:N
                fl=diskItems(i);
                fname=fl.name;
                if fl.isdir
                    if ~strcmp('.', fname) && ~strcmp('..', fname)
                        path=fullfile(folder, fname);
                        [bytes_, items_]=File.DiskUsage(...
                            path, spec, olderThanDays, pu, cur);
                        bytes=bytes+bytes_;
                        items=items+items_;
                    end
                    if ~isempty(pu) && pu.cancelled
                        return;
                    end
                end
            end            
        end
        
        function subFolders=ParseLsStdOut(std)
            subFolders={};
            xxx=regexp(std, '([^ \t\n\r]*)', 'tokens');
            N=length(xxx);
            for i=1:N
                fn=xxx{i}{1};
                if ~String.EndsWith(fn, './')
                    subFolders{end+1}=fn;
                end
            end
        end
        
        function subFolders=GetSubFolders(folder)
            subFolders={};
            if ismac
                [~, std1]=system(['ls -a -d ' folder '/*/']);
                [~, std2]=system(['ls -a -d ' folder '/.*/']);
                subFolders=[ File.ParseLsStdOut(std1) ...
                    File.ParseLsStdOut(std2)];
            end
        end
        
        function Touch(fn)
            java.io.File(fn).setLastModified(java.lang.System.currentTimeMillis);
        end
        
        function [status, stdout]=Spawn(scriptLines, scriptFile, ...
                terminalName, runInBackground, showSpawnWindow)
            if nargin<5
                showSpawnWindow=false;
                if nargin<4
                    runInBackground=true;
                end
            end
            wnd=get(0, 'currentFig');
            [~, scriptF, scriptExt]=fileparts(scriptFile);
            doneFile=fullfile(File.Home, [scriptF scriptExt '.done']);
            if exist(doneFile, 'file')
                delete(doneFile);
            end
            makeDone=['echo > ' String.ToSystem(doneFile)];
            if ismac         
                setTerminalName=['echo -n -e "\033]0;' terminalName '\007"'];
                closeTerminal=['osascript -e ''tell application '...
                    '"Terminal" to close (every window whose name '...
                    'contains "' terminalName '")'' &'];
                if iscell(scriptLines)
                    strs=[setTerminalName, scriptLines, makeDone, ...
                        closeTerminal, 'exit'];
                else
                    strs={setTerminalName, scriptLines, makeDone, ...
                        closeTerminal, 'exit'};
                end
                File.WriteTextFile(scriptFile, strs);
                scriptCmd=String.ToSystem(scriptFile);
                system(['chmod 777 ' scriptCmd]);
                if showSpawnWindow || runInBackground
                    cmd=['open -b com.apple.terminal ' scriptCmd];
                else
                    cmd=scriptCmd;
                end
            else
                if iscell(scriptLines)
                    N=length(scriptLines);
                    strs=cell(1,N+2);
                    for i=1:N
                        strs{i}=[scriptLines{i} ' < nul'];
                    end
                    strs{end-1}=makeDone;
                    strs{end}='exit';
                else
                    strs={[scriptLines ' < nul'], makeDone, 'exit'};
                end
                File.WriteTextFile(scriptFile, strs);
                scriptCmd=String.ToSystem(scriptFile);
                if showSpawnWindow || runInBackground
                    cmd=[scriptCmd ' &'];
                else
                    cmd=scriptCmd;
                end
            end
            if ~runInBackground 
                disp(terminalName);
            end
            [status, stdout]=system(cmd);
            if ~runInBackground
                if showSpawnWindow
                    File.Wait(doneFile, wnd, [], terminalName);
                end
            end            
        end
        
        function Wait(outFile, fig, btn, txt)
            if isempty(fig)
                fig=0;
            end
            setappdata(fig, 'canceling', 0);
            pu=PopUp(txt, 'west',...
                'Stand by, patience ....', false, ...
                @(h,e)cancel);
            set(pu.dlg, 'WindowClosingCallback', @(h,e)closeWnd());
            
            set(btn, 'visible', 'off');
            drawnow;
            done=false;
            was=pause('query');
            pause('on');
            canceled=0;
            closed=false;
            while ~done
                if ~ishandle(fig)
                    return;
                end
                try
                    if exist(outFile, 'file')
                        done=true;
                    elseif getappdata(fig, 'canceling')
                        canceled=canceled+1;
                        break;
                    elseif canceled>3 % hung?
                        break;
                    end
                    if closed
                        pu.dlg.dispose;
                        break;
                    end
                catch ex
                    ex.getReport
                    pu.close;
                    return;
                end
                pause(2);
            end
            pause(was);
            setappdata(fig, 'canceling', 0);
            set(btn, 'visible', 'on');
            pu.close;
            
            function cancel
                if canceled==0
                    pu.setText('Cancelling this process');
                end
                setappdata(fig, 'canceling', 1);
            end
            function closeWnd
                if canceled==0
                    cancel;
                else
                    closed=true;
                end
            end
        end
        

        function DeleteOld(folder, spec, olderThanDays, pu, fileCnt)
            txt1=sprintf('&gt; %s days old', String.encodeInteger(olderThanDays));
            txt2=sprintf('<html>Finding %s files %s in %s</html>', spec, txt1, File.Html(folder));
            if nargin<4
                pu=PopUp(txt2, 'center', 'Deleting ...', true, true);
            end
            if File.DEBUG>0
                disp(txt2)
            end
            [bytes, tm]=File.DirRecursive(folder, spec);
            if nargin<5
                fileCnt=tm.size;
            end
            if File.DEBUG>0
                fprintf('%d files occupying %s in subfolders below %s\n',...
                    fileCnt, String.encodeBytes(bytes), folder);
            end
            File.DeleteByAge(tm, 0, bytes, olderThanDays, pu, fileCnt);
            if nargin<4
                pu.close;
            end
        end
        
        function [fileCnt, bytes]=DeleteByAge(tm, byteLimit, ...
                totalBytes, olderThanDays, pu, totalFiles)
            it=tm.map.descendingKeySet.iterator;
            fileCnt=0;
            txtByteLimit=String.encodeBytes(byteLimit);
            if nargin>4
                if nargin==5
                    totalFiles=tm.size;
                end
                pu.initProgress(totalFiles);
                reportAt=ceil(totalFiles/25);
            end
            while it.hasNext
                days=it.next;
                if olderThanDays>0 && days<=olderThanDays
                    return;
                end
                fileIt=tm.getIterator(days);
                while fileIt.hasNext
                    if byteLimit >0 && totalBytes<byteLimit
                        return;
                    end
                    path=fileIt.next;
                    diskItem=dir(path);
                    bytes=diskItem(1).bytes;
                    totalBytes=totalBytes-bytes;
                    fileCnt=fileCnt+1;
                    if File.DEBUG>0
                        fprintf(['#%s %s days old, %s bytes leaving ' ...
                            '%s/%s\n\t\t--> ..%s\n'], ...
                            String.encodeInteger(fileCnt), ...
                            String.encodeRounded(days),...
                            String.encodeBytes(bytes),...
                            String.encodeBytes(totalBytes),...
                            txtByteLimit,...
                            String.RemoveXml(File.Html(path)));
                    end
                    if nargin>4
                        if mod(fileCnt, reportAt)==0
                            pu.incrementProgress(reportAt);
                        end
                    end
                    delete(path);
                end
            end
        end
        
        function ExportToolBar(tb, gt, labels, data, filterCols, ...
                file, subFolder, property, jtable)
            if nargin<6
                jtable=[];
            end
            ToolBarMethods.addButton(tb, ...
                'save16.gif', 'Export to excel and other file formats', ...
                @(h,e)export(h));
            
            function export(h)
                
                if ~isempty(jtable)
                    filterRows=SortTable.ModelRows(jtable);
                    filterCols2=SortTable.ModelCols(jtable);
                    if ~isempty(filterCols)
                        l=ismember(filterCols2, filterCols);
                        filterCols2=filterCols2(l);
                    end
                    labels2=labels(filterCols2);
                    data2=data(filterRows, filterCols2);
                elseif ~isempty(filterCols)
                    labels2=labels(filterCols);
                    data2=data(:,filterCols);
                else
                    labels2=labels;
                    data2=data;
                end
                File.Export([labels2;data2], gt, file, subFolder, property);
            end
        end
        
        function file=Export(tabLines, gtOrRootFldr, suggestedFile, ...
                subFolder, property, hasImgs, colSplit, rowSplit, cmp)
            if nargin<9
                cmp=[];
                if nargin<8
                    rowSplit=1;
                    if nargin<7
                        colSplit=1;
                        if nargin<6
                            hasImgs=false;
                            if nargin<5
                                property='statsFolder';
                                if nargin<4
                                    subFolder='exports';
                                    if nargin<3
                                        suggestedFile='exported';
                                        if nargin<2
                                            gtOrRootFldr=[];
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if ~isa(gtOrRootFldr, 'GatingTree')
                rootFolder=gtOrRootFldr;
                app=CytoGate.Get;
                props=app;
            else
                gt=gtOrRootFldr;
                rootFolder=gt.app.stripIfUnderUrlFolder(gt.rootFolder, ...
                    CytoGate.Folder);
                props=gt.multiProps;
                app=gt.app;
            end
            if iscell(tabLines)
                tabLines=CellBasics.ToTabLines(tabLines);
            end
            app.currentJavaWindow=cmp;
            CytoGate.setHelp('AG_Export');
            exportType=props.getNumeric(File.PROP_EXPORT, 1);
            types=File.EXPORT_TYPES;
            N=size(types, 1);
            types2={};            
            types2(end+1,:)=types(exportType,:);
            for i=1:N
                if i~=exportType
                    types2(end+1,:)=types(i,:);
                end
            end
            if isempty(rootFolder)
                rootFolder=File.Home;
            end
            File.mkDir(fullfile(rootFolder, subFolder));
            file=File.GetOutputFile(rootFolder, subFolder, suggestedFile,...
                props, property, types2, 'Export to ...');
            if ~isempty(file)
                if String.EndsWithI(file, '.csv')
                    stObj=javaObjectEDT('edu.stanford.facs.swing.Basics');
                    data=StringArray.Cell(stObj.tabToCsv(tabLines));
                    props.set(File.PROP_EXPORT, '2');
                else
                    data=StringArray.Cell(tabLines);
                    props.set(File.PROP_EXPORT, '1');
                end
                if String.EndsWithI(file, '.xls')
                    props.set(File.PROP_EXPORT, '3');
                    if hasImgs
                        embedGateImg=props.is(File.PROP_XLS_IMG, false);
                        [embedGateImg, cancelled]=askYesOrNo(...
                            'Embed images into the worksheet?',...
                            'Confirm...', 'North', embedGateImg, ...
                            [File.PROP_XLS_IMG 'Rid']);
                        if cancelled
                            app.currentJavaWindow=[];
                            return;
                        end                         
                        props.setBoolean(File.PROP_XLS_IMG, embedGateImg);
                        if embedGateImg
                            gateImgFolder=fullfile(rootFolder, CytoGate.Folder);
                        else
                            gateImgFolder=[];
                        end
                    else
                        gateImgFolder=[];
                    end
                    file2=[file '.txt'];
                    File.WriteTextFile(file2, data);
                    if exist(file, 'file')
                        dflt=props.getNumeric(File.PROP_XLS_ADD, 1);
                        choices={'Add a new worksheet to the workbook', ...
                            'Overwrite (erase & recreate) the workbook'};
                        [answ,  cancelled]=mnuMultiDlg(...
                            'This workbook file exists ...', 'Caution...', ...
                            choices, dflt, true, true);
                        if cancelled
                            return;
                        end
                        if answ==2
                            delete(file);
                        end
                        props.set(File.PROP_XLS_ADD, num2str(answ));
                    end
                    makeXlsExternal(file2, file, colSplit, rowSplit, gateImgFolder);
                    %delete(file2);
                    openXlsNow=props.is(File.PROP_XLS_OPEN, false);
                    app.currentJavaWindow=[];
                    if openXlsNow
                        if ismac
                            system(['open ' String.ToSystem(file)]);
                        else
                            system(String.ToSystem(file));
                        end
                    end
                    delete([file '.txt']);
                else
                    File.WriteTextFile(file, data);
                end
                File.OpenFolderWindow(file);
            end
            app.currentJavaWindow=[];     
        end

        function OpenFolderWindow(file)
            if ~isempty(file)
                app=CytoGate.Get;
                if askYesOrNo(Html.WrapHr(['Open folder window now?'...
                        '<br><br>' app.smallStart file app.smallEnd ...
                        ')']), 'Confirm...', 'center', true, 'openFolder')
                    fldr=fileparts(file);
                    fldr=String.ToSystem(fldr);
                    if ispc
                        system(['explorer ' fldr]);
                    else
                        system(['open ' fldr]);
                    end
                end
            end
        end
        
        function [file, folder, fileName]=...
                GetOutputFile(rootFolder, suggestedSubFolder, fileName, ...
                properties, property, types, tip)
            if nargin==0
                rootFolder=File.Home;
                suggestedSubFolder='outputs';
                fileName='output';
                properties=[];
            end
            app=CytoGate.Get;
            parentFolder=app.stripIfUnderUrlFolder(rootFolder, '.autoGate');
            if ~isempty(suggestedSubFolder)
                suggestedSubFolder=fullfile(parentFolder, suggestedSubFolder);
            else
                suggestedSubFolder=parentFolder;
            end
            if ~isempty(properties)
                fullFolder=properties.get(property, suggestedSubFolder);
            else
                fullFolder=suggestedSubFolder;
            end
            fullFile=fullfile(fullFolder, fileName);
            fullFile=File.SwitchExtension(fullFile, types{1,1});
            File.mkDir(suggestedSubFolder);
            [fileName, folder]=uiputfile(types, tip, fullFile);
            file=fullfile(folder,fileName);
            if ~isnumeric(folder) && ~isempty(properties)
                properties.set(property, folder);
            else
                file=[];
                folder=[];
                fileName=[];
            end
        end
        
        function names=CsvNames(t)
            names={};
            try
                names=t.Properties.VariableNames;
                N=length(names);
                vd=t.Properties.VariableDescriptions;
                if length(vd)==N
                    prefix='Original column heading: ''';
                    idx=length(prefix)+1;
                    for i=1:N
                        dsc=vd{i};
                        if String.StartsWith(dsc, prefix)
                            str=dsc(idx:end-1);
                            N2=length(str);
                            for j=1:N2
                                firstChar=double(str(j));
                                if firstChar>=30 && firstChar<=128
                                    break;
                                end
                            end
                            if j>1 && j<=N2
                                str=str(j:end);
                            end
                            names{i}=str;
                        else
                            names{i}=dsc;
                        end
                    end
                end
            catch ex
                ex.getReport
            end
        end
        
        function [inData, columnNames]=ReadCsv(csvFile, mustBeNumbers)
            if nargin<2
                mustBeNumbers=true;
            end
            warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
            t=readtable(csvFile, 'ReadVariableNames', true);
            last1Bad=false;
            try
                inData=table2array(t);
                if verLessThan('matlab', '9.6')
                    if all(isnan(inData(:,end)))
                        C=size(inData,2);                        
                        if isequal(['Var' num2str(C)], ...
                                t.Properties.VariableNames{end}) ...
                                && isequal('', ...
                                t.Properties.VariableDescriptions{end})
                            last1Bad=true;
                            inData(:,end)=[];
                        end
                    end
                end
            catch ex
                if verLessThan('matlab', '9.6')
                    inData=[];
                    columnNames=[];
                    ex.getReport
                    return;
                end
                last1Bad=false;
                t=removevars(t, t.Properties.VariableNames{end});
                inData=table2array(t);
            end
            columnNames=File.CsvNames(t);
            if last1Bad
                columnNames(end)=[];
            end
            if iscell(inData) && mustBeNumbers
                inData=[];
                columnNames=[];
                warning('comma separated file does NOT contain numbers');
            end
        end
        
        function out=Time
            out=char(datetime('now','Format','d-MMM-y HH_mm_ss'));
        end
        
        
        function [folder, file]=PutFile(dfltFolder, dfltFile, props, ...
                property, ttl, longTtl, ext)
            if nargin<7
                ext=[];
                if nargin<6
                    longTtl='Save to which folder & file?';
                    if nargin<5
                        ttl='Save file...';
                        if nargin<4
                            property=[];
                            if nargin<3
                                props=[];
                            end
                        end
                    end
                end
            end
            privFldr=fileparts(tempname);
            [lastFolder, fn, fe]=fileparts(dfltFile);
            if isempty(ext)
                ext=fe;
            end
            File.mkDir(dfltFolder);
            if isempty(lastFolder) || isequal(privFldr, lastFolder)
                fn=File.Time;
                if ~isempty(props) && ~isempty(property)
                    lastFolder=props.get(property, dfltFolder);
                else
                    lastFolder=dfltFolder;
                end
            end
            lIdx=String.LastIndexOf(lastFolder, File.Home);
            if lIdx>2
                lastFolder=lastFolder(lIdx:end);
            end
            if endsWith(lastFolder, '.autoGate')
                lastFolder=fileparts(lastFolder);
            end
            done=false;
            if ismac
                jd=Gui.MsgAtTopScreen(longTtl, 25);
            else
                jd=[];
            end
            while ~done
                done=true;
                [file, folder]=uiputfile(['*' ext], ttl, ...
                    fullfile(lastFolder, [fn ext]));
                if ~isempty(jd)
                    jd.dispose;
                end
                if isempty(folder) || isnumeric(folder)
                    folder=[];
                    file=[];
                    if isequal(dfltFolder, lastFolder)
                        return;
                    end
                    if isequal([dfltFolder filesep], lastFolder)
                        return;
                    end
                    if isequal(dfltFolder, [lastFolder filesep])
                        return;
                    end
                    app=BasicMap.Global;
                    checkDflt=true;
                    if ~isempty(props) && ~isempty(property)
                        propFldr=props.get(property, dfltFolder);
                        if ~isequal(propFldr, dfltFolder)
                            choices{1}=dfltFolder;
                            choices{end+1}=propFldr;
                            answer=Gui.Ask('Check other locations?', choices, ...
                                [property '.other'], 'Confirm...');
                            if isempty(answer)
                                folder=[];
                                file=[];
                                return;
                            end
                            checkDflt=false;
                            dfltFolder=choices{answer};
                        end
                    end
                    if checkDflt
                        [~,yes]=questDlg(...
                            ['<html><center>Look in default folder?<br>'...
                            '(<b>' app.smallStart dfltFolder '</b>)' app.smallEnd...
                            '<hr></center></html>']);
                        if ~yes
                            return;
                        end
                    end
                    [file, folder]=uiputfile(['*' ext], ...
                        'Save to which folder & file?', ...
                        fullfile(dfltFolder, [fn ext]));
                    if isempty(folder)|| isnumeric(folder)
                        folder=[];
                        file=[];
                        return;
                    end
                end
            end
            if ~isempty(props) && ~isempty(property)
                props.set(property, folder);
            end
            file=File.SwitchExtension2(file, ext);
        end
    end
end
