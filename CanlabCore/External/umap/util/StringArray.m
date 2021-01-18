%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

classdef StringArray < handle
    properties(SetAccess=private)
        strings;
        N;
    end
    
    methods
        function this=StringArray(strings)
            this.strings=strings;
            this.N=length(strings);
        end
        function idx=lastIndexOf(this, search)
            idx=0;
            for i=this.N:-1:1
                if strcmp(this.strings{i}, search)
                    idx=i;
                    break;
                end
            end
        end
        
        function ok=contains(this,search)
            ok=this.indexOf(search)>0;
        end
        
        function idx=indexOf(this, search)
            idx=0;
            for i=1:this.N
                if strcmp(this.strings{i}, search)
                    idx=i;
                    break;
                end
            end
        end
        
        function idx=endsWith(this, with)
            idx=0;
            for i=1:this.N
                if endsWith(this.strings{i}, with)
                    idx=i;
                    break;
                end
            end
        end
        
        function indices=indicesOf(this, search)
            indices=[];
            for i=1:this.N
                if strcmp(this.strings{i}, search)
                    indices(end+1)=i;
                end
            end
        end

        function idx=indexOfI(this, search)
            idx=0;
            for i=1:this.N
                if strcmpi(this.strings{i}, search)
                    idx=i;
                    break;
                end
            end
        end
        
    end
    
    methods(Static)
        function ok=IsOk(x)
            ok=true;
            if iscell(x)
                N=length(x);
                for i=1:N
                    if ~ischar(x{i})
                        ok=false;
                        break;
                    end
                end
            elseif ischar(x)
            else
                ok=false;
                warning('Expecting string or string array');
            end
        end
        
        function strs=Trim(strs, limit)
            if length(strs)>limit
                strs{limit}=[num2str(length(strs)-limit) ' more...'];
                strs(limit+1:end)=[];
            end
        end
        function c=Num2Str(nums, prefix)
            if nargin<2
                prefix='#';
            end
            N_=length(nums);
            c=cell(1, N_);
            for i_=1:N_
                c{i_}=[prefix num2str(nums(i_))];
            end
        end
        
        function out=Intersection(left, right)
            if isempty(left)
                out=right;
            elseif isempty(right)
                out=left;
            else
                out={};
                r=StringArray.Set(right);
                N=length(left);
                for i=1:N
                    try
                        if r.contains(java.lang.String(left{i}))
                            out{end+1}=left{i};
                        end
                    catch ex
                    end
                end
            end
        end
        
        function out=Sort(in, canon)
            if nargin==1
                out=StringArray.Cell(StringArray.Set(in));
            else
                N=length(in);
                this=StringArray(canon);
                map=java.util.TreeMap;
                for i=1:N
                    str=in{i};
                    idx=this.indexOf(str);
                    if idx>0
                        map.put(idx, str);
                    else
                        map.put(N+i, str);
                    end
                end
                out=StringArray.Cell(map.values);
            end
        end
        
        function c=Cell(javaCollection)
            N=javaCollection.size;
            c=cell(1,N);
            it=javaCollection.iterator;
            for i=1:N
                c{i}=it.next;
            end
        end

        function out=JavaString2Chars(in)            
            N=length(in);
            out=cell(1,N);
            for i=1:N
                out{i}=char(in{i});
            end
        end

        function c=JavaArrayToCell(javaArr)
            if iscell(javaArr)
                c=javaArr;
            else
                N=length(javaArr);
                c=cell(1,N);
                for i=1:N
                    c{i}=char(javaArr(i));
                end
            end
        end

        function lst=List(strs)
            lst=java.util.ArrayList;
            N=length(strs);
            for i=1:N
                lst.add(strs{i});
            end
        end
        
        function lst=Vector(strs)
            lst=java.util.Vector;
            N=length(strs);
            for i=1:N
                lst.add(strs{i});
            end
        end
                
        function lst=ListFromJavaArray(strs)
            lst=java.util.ArrayList;
            N=length(strs);
            for i=1:N
                lst.add(strs(i));
            end
        end
        function set=Set(strs, set)
            if nargin<2
                set=java.util.TreeSet;
            end
            N=length(strs);
            for i=1:N
                if ~isempty(strs{i})
                    set.add(java.lang.String(strs{i}));
                end
            end
        end
        
        function out=RemoveDuplicates(strs, set)
            if nargin<2
                set=StringArray.Set(strs);
            else
                set=StringArray.Set(strs, set);
            end
            N_=set.size;
            out=cell(1,N_);
            it=set.iterator;
            for i=1:N_
                out{i}=it.next;
            end
        end
        
        function [out, dupCnt]=ForceUnique(strs, prefix, suffix)
            dupCnt=0;
            map=java.util.HashMap;
            N_=length(strs);
            out=cell(1,N_);
            for i=1:N_
                if map.containsKey(strs{i})
                    cnt=map.get(strs{i})+1;
                else
                    cnt=1;
                end
                map.put(strs{i}, cnt);
                if cnt>1
                    out{i}=[strs{i} prefix num2str(cnt) suffix];
                    dupCnt=dupCnt+1;
                else
                    out{i}=strs{i};
                end
            end
        end
        
        function ok=Contains(searchedStrings, search)
            sa=StringArray(searchedStrings);
            ok=sa.indexOf(search)>0;
        end
        function idx=IndexOf(searchedStrings, search)
            sa=StringArray(searchedStrings);
            idx=sa.indexOf(search);
        end

        function idxs=IndexOfEmpties(searchedStrings)
            N=length(searchedStrings);
            idxs=[];
            for i=1:N
                if isempty(searchedStrings{i})
                    idxs(end+1)=i;
                end
            end
        end

        function idx=IndexOfIgnoreCase(searchedStrings, search)
            sa=StringArray(searchedStrings);
            idx=sa.indexOfI(search);
        end

        function idxs=FirstListIndexes(searched, searches)
            N=length(searches);
            idxs=zeros(1,N);
            for i=1:N
                idxs(i)=searched.indexOf(searches{i});
            end
            idxs=idxs+1;
        end
        
        function [idxs, foundAll]=EndsWith(searchedStrings, searches)
            idxs=[];
            sa=StringArray(searchedStrings);
            N=length(searches);
            for i=1:N
                idx=sa.endsWith(searches{i});
                if idx>0
                    idxs=[idxs idx+1];
                else
                    disp('huh');
                end
            end
            foundAll=length(idxs)>=length(searchedStrings);
        end

        %Not sure if a faster search than this is available at
        %https://www.mathworks.com/matlabcentral/answers/335930-most-efficient-way-to-search-in-text-arrays
        %I think that the the problem with find(strcmp(a,s),1) is that
        %strcmp does not STOP after its first successful strcmp
        %the Containers.map approach WILL ... but has a high startup
        %cost for under 1000 strings
        function idxs=IndexesOf(searchedStrings, searches)
            idxs=[];
            sa=StringArray.List(searchedStrings);
            N=length(searches);
            for i=1:N
                idx=sa.indexOf(searches{i});
                if idx>=0
                    idxs=[idxs idx+1];
                else
                    disp('huh');
                end
            end
        end

        function idxs=IndexesOf2(searchedStrings, searches)
            sa=StringArray.List(searchedStrings);
            N=length(searches);
            idxs=zeros(1,N);
            for i=1:N
                idxs(i)=sa.indexOf(searches{i})+1;
            end
        end

        function idxs=NotFound(searchedStrings, searches)
            idxs=[];
            sa=StringArray.List(searchedStrings);
            N=length(searches);
            for i=1:N
                idx=sa.indexOf(searches{i});
                if idx<0
                    idxs=[idxs i];
                end
            end
        end

        function [ok, sameButReordered]=AreSameOrEmpty(d1, d2)
            ok=isequal(d1, d2);
            sameButReordered=false;
            if ~ok
                if length(d1)==length(d2)
                    nf=StringArray.NotFound(d2, d1);
                    if isempty(nf)
                        sameButReordered=true;
                    elseif length(nf)==1
                        % is FMO?
                        if isempty(d1{nf}) || isempty(d2{nf})
                            d1_=d1;
                            d2_=d2;
                            d1_{nf}='dummy';
                            d2_{nf}='dummy';
                            ok=isequal(d1_,d2_);
                        end
                    end
                end
            end
        end
        
        function [allFound, idxs]=Find(d1, d2, allowOneEmptyLikeFMO)
            if nargin<3
                allowOneEmptyLikeFMO=false;
            end
            allFound=true;
            N=length(d1);
            idxs=zeros(1,N);
            
            sa2=StringArray(d2);
            for i=1:N
                ii=sa2.indexOf(d1{i});
                if ii<1
                    if allowOneEmptyLikeFMO
                        allowOneEmptyLikeFMO=false;
                        [emptyCount, ii]=StringArray.GetEmpties(d2);
                        if emptyCount ~= 1
                            [emptyCount, i2]=StringArray.GetEmpties(d1);
                            if emptyCount ~= 1 
                                allFound=false;
                                ii=0;
                            else
                                ii=0;
                                if i2==i
                                    sa1=StringArray(d1);
                                    N2=length(d2);
                                    for j=1:N2
                                        jj=sa1.indexOf(d2{j});
                                        if jj<1
                                            ii=j;
                                            break;
                                        end
                                    end
                                end
                                if ii<1
                                    allFound=false;
                                end
                            end
                        end
                    else
                        allFound=false;
                    end
                end
                if ~isempty(ii)
                    idxs(i)=ii(1);
                end
            end
        end

        
        function str=toString(stringArray, delimiter, useQuotes)
            if nargin<3
                useQuotes=false;
                if nargin<2
                    delimiter=', ';
                end
            end
            N=length(stringArray);
            str='';
            for i=1:N
                v=stringArray{i};
                if ~ischar(v)
                    v=String.toString(v);
                end
                if useQuotes
                    v=['"' v '"'];
                end
                if (i>1)
                    str=[str delimiter v];
                else
                    str=v;
                end
            end
        end
        
        function ok=isEqual(stringArray1, stringArray2)
            N1=length(stringArray1);
            N2=length(stringArray2);
            ok=false;
            if N1==N2
                ok=true;
                for i=1:N1
                    if ~strcmp(stringArray1{i}, stringArray2{i})
                        ok=false;
                        break;
                    end
                end
            end
        end
    
        
        function a=Join(a, a2)
            N=length(a2);
            for i=1:N
                a{end+1}=a2{i};
            end
        end
        
        function stringArray=Subset(strings, idxs)
            stringArray=StringArray( strings(idxs) );
        end
        
        function [searchResultIdxs, unSearchedStrings, unFoundStrings]=...
                SearchSubset(searchedStrings, subsettingIndices, searchStrings)
            N=length(searchStrings);
            unFoundStrings={};
            searchResultIdxs=1:N;
            subset=StringArray.Subset(searchedStrings, subsettingIndices);
            unSearchedStrings=subset.strings;
            for i=N:-1:1                
                idx=subset.indexOf(searchStrings{i});
                searchResultIdxs(i)=idx;
                if idx>0
                    %remove found string
                    unSearchedStrings(idx)=[];
                else 
                    unFoundStrings{end+1}=searchStrings{i};
                end
            end
        end
        
        function [searchResultIdxs, unsearchedStrings, unFoundStrings]=...
                Search2Subsets(searchedStrings1, searchedStrings2,...
                subsettingIndices, searchStrings)
            N=length(searchStrings);
            
            unFoundStrings={};
            searchResultIdxs=1:N;
            subset1=StringArray.Subset(searchedStrings1, subsettingIndices);
            subset2=StringArray.Subset(searchedStrings2, subsettingIndices);
            temp=subset1.strings;
            for i=N:-1:1       
                searchString=searchStrings{i};
                idx=subset1.indexOf(searchString);
                if idx==0
                    idx=subset2.indexOf(searchString);
                    
                end
                searchResultIdxs(i)=idx;
                if idx>0
                    %remove found string
                    temp{idx}='';
                else 
                    unFoundStrings{end+1}=searchString;
                end
            end
            unsearchedStrings={};
            for i=1:length(temp)
                if ~isempty(temp{i})
                    unsearchedStrings{end+1}=temp{i};
                end
            end
        end
        
        
         function[files]=getField(handles, fileField)
             files={};
             if isfield(handles, fileField)
                 temp=getfield(handles,fileField);
                 if ~strcmp('cell', class(temp))
                     files={temp};
                 else
                     files=temp;
                 end
             end
         end
         
         function array=toArray(something)
             if ~strcmp('cell', class(something))
                 array={something};
             else
                 array=something;
             end
         end
         
         function s=encodeCell(c, quotes)
             if nargin==1
                 quotes=false;
             end
             N=length(c);
             s='';
             for i=1:N
                 if i>1
                     s=[s ', '];
                 end
                 if quotes
                     s=[s '''' c{i} ''''];
                 else
                     s=[s c{i}];
                 end
             end
             if quotes
                 s=['{' s '}'];
             end
         end
         
         function c=decodeCell(s)
             if String.StartsWith(s, '{''')
                 c=eval(s);
             else
                 c=regexp(s, '[^,\s]*', 'match');
             end
         end
         
         function idxs=StartsWith(sa, str)
             idxs=[];
             N=length(sa);
             for i=1:N
                 if startsWith(str, sa{i})
                     idxs(end+1)=i;
                 end
             end
         end
         
         function ok=EqualsIndexed(these, those, idxs)
             try
                 N_=length(idxs);
                 N2=length(these);
                 N3=length(those);
                 ok=true;
                 for i=1:N_
                     idx=idxs(i);
                     if idx>N2 || idx>N3
                         continue;
                     end
                     if ~strcmp(these{idx}, those{idx})
                         ok=false;
                         break;
                     end
                 end
             catch exception
                 warning(exception.message);
                 ok=false;
             end
         end
         
         function ok=Equals(these, those)             
             N1=length(these);
             N2=length(those);
             if N1 ~= N2
                 ok=false;
             else
                 ok=StringArray.EqualsIndexed(these, those, [1:N1]);
             end
         end
         
         function found=HasAny(superSet, subSet, superSetValues)
             found=false;
             for i=1:length(subSet)
                 idx=StringArray.IndexOf(superSet, subSet{i});
                 if idx>0
                     if nargin<2 || ~isempty(superSetValues{idx})
                         found=true;
                         break;
                     end
                 end
             end
         end
         
         function found=HasAll(superSet, subSet, superSetValues)
             found=true;
             for i=1:length(subSet)
                 idx=StringArray.IndexOf(superSet, subSet{i});
                 if idx>0
                     if nargin>1 && isempty(superSetValues{idx})
                         found=false;
                         break;
                     end
                 else
                     found=false;
                     break;
                 end
             end
         end
     
         function [arrayOfArrayIdx, arrayIdx]=FindAofA(arrayOfArrays, searchValue)
             arrayOfArrayIdx=0;
             arrayIdx=0;
             N=length(arrayOfArrays);
             for i=1:N
                 j=StringArray.IndexOf(arrayOfArrays{i}, searchValue);
                 if j>0
                     arrayIdx=j;
                     arrayOfArrayIdx=i;
                     return;
                 end
             end
         end
         
         function ok=ContainsInAofA(arrayOfArrays, arrayIdx)
             ok=false;
             N=length(arrayOfArrays);
             if N>=arrayIdx
                 j=StringArray.IndexOf(arrayOfArrays{arrayIdx}, searchValue);
                 if j>0
                   ok=true;                   
                 end
             end
         end
         
         function [lbl, sz]=GetLongest(strs)
             N=length(strs);
             sz=0;
             lbl=[];
             for i=1:N
                 ln=length(strs{i});
                 if ln>sz
                     sz=ln;
                     lbl=strs{i};
                 end
             end
         end
         
         function out=ToLower(strs)
             N=length(strs);
             out=cell(1,N);
             for i=1:N
                 out{i}=lower(strs{i});
             end
         end
         
         function [cnt, idxs]=GetEmpties(strs)
             cnt=0;
             idxs=[];
             N=length(strs);
             for i=1:N
                 if isempty(strs{i})
                     cnt=cnt+1;
                     idxs(end+1)=i;
                 end
             end
         end
         
         function ok=IsValid(in)
             if iscell(in)
                 ok=true;
                 N=length(in);
                 for i=1:N
                     if ~ischar(in{i})
                         ok=false;
                         break;
                     end
                 end
             else
                 ok=false;
             end
             
         end
    end
end
