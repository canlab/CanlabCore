classdef ColorsByName <handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
    
    properties(Constant)
        DEFAULT_COLOR=[.9 .9 .9];
    end
    
    properties
        file;
        props;
        synonyms;
        canon;
        map;
        fnTranslate;
    end
    methods
        function this=ColorsByName(file, fnTranslate)
            if nargin<2
                fnTranslate=[];
                if nargin<1
                    propFile='colorsByName.properties';
                    file=fullfile(UmapUtil.LocalSamplesFolder, propFile);
                    if exist(file, 'file')
                    elseif exist('~/.autoGate/colorsByName.properties', 'file')
                        file='~/.autoGate/colorsByName.properties';
                    else
                        url=WebDownload.ResolveUrl(propFile);
                        WebDownload.Get({url}, {file}, false, false, 'south');
                    end
                end
            end
            
            this.file=file;
            if isempty(fnTranslate)
                this.fnTranslate=@translate;
            else
                this.fnTranslate=fnTranslate;
            end
            synonyms=TreeMapOfMany;
            canon=java.util.HashMap;
            map=BasicMap;
            this.synonyms=synonyms;
            this.canon=canon;
            this.map=map;
            colorsByName=ColorsByName.ReadProperties(file);
            if ~isempty(colorsByName)
                this.props=colorsByName;
                it=colorsByName.keySet.iterator;
                while it.hasNext
                    name=it.next;
                    try
                        isCanonicalName=startsWith(name, '*');
                        strColor=char(colorsByName.get(name));
                        if isCanonicalName
                            name=name(2:end);
                        end
                        if ~isempty(strColor)
                            strColor=str2num(strColor);
                            color=strColor/256;
                            strColor=String.Num2Str(strColor, ' ');
                            if length(color)==3
                                if any(color>1|color<0)
                                    warning(['Color(s) > 256 for "' name '": [' strColor ']']);
                                    color(color>1)=1;
                                    color(color<0)=0;
                                end
                                map.set(lower(name), color);
                            else
                                warning(['Expecting 3 color numbers for "' name '": [' strColor ']']);
                                map.set(lower(name), ColorsByName.DEFAULT_COLOR);
                            end
                            strColor=java.lang.String(strColor);
                            synonyms.put(strColor, name);
                            if isCanonicalName || ~canon.containsKey(strColor)
                                canon.put(strColor, name);
                            end
                            if synonyms.size(strColor)>1
                                disp([name ' has ' num2str(synonyms.size(strColor)) ' synonyms']);
                            end
                        else
                            warning(['Null color for "' name '": [' strColor ']']);
                            map.set(lower(name), ColorsByName.DEFAULT_COLOR);
                        end
                    catch 
                        fprintf('Name #%d "%s" is null?\n', i, name);
                    end
                end
            else
                this.props=java.util.Properties;
            end
            function out=translate(in)
                out=in;
            end
        end
        
        function [synonyms, clr, canonicalName, originalName, nameKey]=...
                getSynonyms(this, externalKey)
            originalName=feval(this.fnTranslate, externalKey);
            nameKey=lower(originalName);
            clr=String.Num2Str((floor(256*this.map.get(nameKey))), ' ');
            synonyms=this.synonyms.getCell(clr);
            if nargout>2
                canonicalName=this.canon.get(clr);
            end
        end
        
        function color=get(this, key, idx, nIdxs)
            if isempty(this.props) || this.props.size==0
                color=[];
                return;
            end
            color=this.map.get(key);
            if isempty(color)
                color=this.map.get(lower(key));
                if isempty(color)
                    name=lower(feval(this.fnTranslate, key));
                    color=this.map.get(name);
                    if ~isempty(color)
                        this.map.set(key, color);
                    elseif nargin>2
                        color=Gui.HslColor(idx, nIdxs);
                    end
                end
            end
        end
        
        function yes=isUsing(this, externalKey)
            yes=this.map.containsKey(java.lang.String(externalKey));
        end
        
        function cnt=update(this, externalKey, color)
            cnt=0;
            if isempty(this.props) || this.props.size==0 
                return;
            end
            [names, oldColor, ~, originalName]=...
                getSynonyms(this, externalKey);
            newColor=java.lang.String(String.Num2Str(floor(256*color), ' '));
            if isequal(newColor, oldColor)
                return;
            end
            N=length(names);
            if N>0
                for i=1:N
                    name=names{i};
                    this.map.set(lower(name), color);
                    this.props.put(name, newColor);
                    cnt=cnt+1;
                end
                this.synonyms.renameKey(oldColor, newColor);
                if this.canon.containsKey(oldColor)
                    v=this.canon.remove(oldColor);
                    this.canon.put(newColor, v);
                    this.props.put(['*' v], newColor);
                end
            else
                if this.synonyms.containsKey(newColor)
                    warning(['Color already in use for "' ...
                        originalName '" by "' ...
                        char(this.synonyms.get(newColor)) '"']);
                else
                    cnt=1;
                    this.map.set(lower(originalName), color);
                    this.synonyms.put(newColor, originalName);
                    this.canon.put(newColor, originalName);
                    this.props.put(originalName, newColor);
                end
            end
            this.save;
        end
        
        function save(this)
            File.SaveProperties2(this.file, this.props)
        end
        
        function [yes, synonyms, canonicalName, color, htmlSymb, strColor]=...
                isUsed(this, colorOrKey)
            yes=false;
            synonyms=[];
            canonicalName=[];
            htmlSymb=[];
            strColor=[];
            if isempty(this.props) || isempty(colorOrKey)
                return;
            end
            if isnumeric(colorOrKey) && length(colorOrKey)==3
                color=colorOrKey;
            else
               if ~ischar(colorOrKey)
                   colorOrKey=char(colorOrKey);
               end
               color=this.get(colorOrKey);
               if length(color)~=3
                   return;
               end
            end            
            strColor=String.Num2Str(floor(256*color), ' ');
            yes=this.synonyms.containsKey(strColor);
            if nargout>1
                synonyms=this.synonyms.getCell(strColor);
                if nargout>2
                    canonicalName=this.canon.get(strColor);
                    if nargout>4
                        htmlSymb=Html.Symbol(color, 30, false);
                    end
                end
            end
        end
    end
    
    methods(Static)
        function Update(lblMap, colorFile)
            colorsByName=File.GetProperties(colorFile);
            it=lblMap.keySet.iterator;
            while it.hasNext
                key=char(it.next);
                if ~endsWith(key, '.color')
                    name=String.RemoveTex(lblMap.get(java.lang.String(key)));
                    if ~isempty(name) 
                        clr=lblMap.get([key '.color']);
                        if ~isempty(clr)
                            colorsByName.put(java.lang.String(name), clr);
                        end
                    end
                end
            end
            File.SaveProperties2(colorFile, colorsByName)
        end
        
        function props=ReadProperties(colorFile)
            if ~isempty(colorFile)
                fldr=fileparts(colorFile);
                if isempty(fldr)
                    fl=fullfile(BasicMap.Global.contentFolder, colorFile);
                    if ~exist(fl, 'file')
                        fl=fullfile(fileparts(...
                            BasicMap.Global.contentFolder), 'umap', ...
                            colorFile);
                        if ~exist(fl, 'file')
                            fl=fullfile(BasicMap.Global.contentFolder, ...
                                'umap', colorFile);
                        end
                    end
                    colorFile=fl;
                    
                end
                props=File.ReadProperties(colorFile);
            else
                props=java.util.Properties;
            end
        end
        function cnt=Override(lblMap, colorFile, beQuiet)
            cnt=0;
            if ~isempty(colorFile)
                colorsByName=ColorsByName.ReadProperties(colorFile);
                if ~isempty(colorsByName)
                    c=StringArray.Cell(colorsByName.keySet);
                    N=length(c);
                    for i=1:N
                        key=c{i};
                        try
                            clr=colorsByName.get(key);
                            if startsWith(key, '*') %canonical name
                                colorsByName.put(lower(key(2:end)), clr);
                            else
                                colorsByName.put(lower(key), clr);
                            end
                        catch 
                            if nargin<3 || ~beQuiet
                                fprintf('key #%d "%s" is null?\n', i, key);
                            end
                        end
                    end
                    it=lblMap.keySet.iterator;
                    while it.hasNext
                        key=char(it.next);
                        if ~endsWith(key, '.color')
                            name=lower(String.RemoveTex(lblMap.get(...
                                java.lang.String(key))));
                            if colorsByName.containsKey(name)
                                clr=colorsByName.get(name);
                                lblMap.put([key '.color'], clr);
                                cnt=cnt+1;
                            else
                                if nargin<3 || ~beQuiet
                                    fprintf('No override for %s=%s\n',...
                                        key, name);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
