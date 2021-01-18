%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause

classdef TreeMapOfMany < handle
    properties
        map=[];
        useSet=false;
    end
    methods
        function this=TreeMapOfMany(map)
            if nargin<1
                this.map=java.util.TreeMap;
            else
                this.map=map;
            end
        end
        
        function clear(this)
            this.map.clear;
        end
        
        function putAll(this,k,values)
            N=length(values);
            for i=1:N
                this.put(k, values{i});
            end
        end
        
        function set(this, k, v)
            this.put(k, v);
        end
        
        function k=keys(this)
            k=CellBasics.Java(this.map.keySet);
        end
        function put(this, k, v)
            if ischar(k) && length(k)==1
                k=java.lang.String(k);
            end
            if this.map.containsKey(k)
                l=this.map.get(k);
            else
                if this.useSet
                    l=java.util.HashSet;
                else
                    l=java.util.ArrayList;
                end
                this.map.put(k,l);
            end
            if ischar(v) && length(v)==1
                v=java.lang.String(v);
            end
            l.add(v);
        end
        
        function it=getIterator(this,k)
            if ischar(k) && length(k)==1
                k=java.lang.String(k);
            end
            if this.map.containsKey(k)
                l=this.map.get(k);
                it=l.iterator;
            else
                l=java.util.ArrayList;
                it=l.iterator;
            end
        end
        
        function n=keySize(this)
            n=this.map.size;
        end
        
        function n=size(this, k)
            if nargin<2
                n=0;
                it=this.map.keySet.iterator;
                while (it.hasNext)
                    k=it.next;
                    if ischar(k) && length(k)==1
                        k=java.lang.String(k);
                    end
                    l=this.map.get(k);
                    n=n+l.size;
                end
            else
                n=this.valueCount(k);
            end
        end
        
        function n=valueCount(this, k)
            l=this.map.get(java.lang.String(k));
            if ~isempty(l)
                n=l.size;
            else
                n=0;
            end
        end
        
        function values=getUniqueValues(this)
            values=StringArray.Cell(this.getUniqueValueSet);
        end

        function count=getUniqueValueCount(this)
            count=this.getUniqueValueSet.size;
        end

        function c=getKeyCountStrings(this, hasWord, itemWord)
            if nargin<3
                itemWord='item';
                if nargin<2
                    hasWord='has';
                end
            end
            c={};
            it=this.map.keySet.iterator;
            while it.hasNext
                k=char(it.next);
                c{end+1}=[k ' ' hasWord ' ' ...
                    String.Pluralize2(itemWord, this.valueCount(k))];
            end
        end
        
        function hs=getUniqueValueSet(this, key)
            if nargin>1
                hs=java.util.HashSet;
                it=this.getIterator(key);
                while it.hasNext
                    hs.add(it.next);
                end
                return;
            end
            hs=java.util.HashSet;
            it=this.map.descendingKeySet.iterator;
            n=this.size;
            i=1;
            while it.hasNext
                k=it.next;
                if ischar(k) && length(k)==1
                    k=java.lang.String(k);
                end
                l=this.map.get(k);
                it2=l.iterator;
                while it2.hasNext
                    hs.add(it2.next);
                end
            end
        end

        function all=getAll(this, ascending)
            if nargin<2
                ascending=true;
            end
            if ascending
                it=this.map.keySet.iterator;
            else
                it=this.map.descendingKeySet.iterator;
            end
            n=this.size;
            all=cell(n,2);
            i=1;
            while it.hasNext
                k=it.next;
                if ischar(k) && length(k)==1
                    k=java.lang.String(k);
                end
                l=this.map.get(k);
                it2=l.iterator;
                while it2.hasNext
                    all{i,1}=k;
                    all{i,2}=it2.next;
                    i=i+1;
                end
            end
        end
        
        function [key, value]=getNth(this, nTh, ascending)
            key=[];
            value=[];
            if nargin<3
                ascending=true;
            end
            if ascending
                it=this.map.keySet.iterator;
            else
                it=this.map.descendingKeySet.iterator;
            end
            i=1;
            while it.hasNext
                k=it.next;
                if ischar(k) && length(k)==1
                    k=java.lang.String(k);
                end
                l=this.map.get(k);
                it2=l.iterator;
                while it2.hasNext
                    v=it2.next;
                    if i==nTh
                        key=k;
                        value=v;
                        return;
                    end
                    i=i+1;
                end
            end
        end
        
        function ok=containsKey(this, key)
            ok=this.map.containsKey(key);
        end
        
        function c=getCell(this, key)
            c=this.map.get(java.lang.String(key));
            if isempty(c)
                c={};
            else
                c=CellBasics.Java(c);
            end
        end
        
        function [obj,l]=get(this, key)
            l=this.map.get(java.lang.String(key));
            if isempty(l)
                obj=[];
                if nargout>1
                    if this.useSet
                        l=CytoGate.Get.emptySet;
                    else
                        l=CytoGate.Get.emptyList;
                    end
                end
            else
                obj=l.toString;
            end
        end
        
        function ok=hasString(this, key, str)
            [~,hs]=this.get(key);
            if hs.size>0
                ok=hs.contains(java.lang.String(str));
            else
                ok=false;
            end
        end
        function ok=hasAnyString(this, key, strs)
            [~,hs]=this.get(key);
            if hs.size>0
                N=length(strs);
                for i=1:N
                    if hs.contains(java.lang.String(strs{i}))
                        ok=true;
                        return;
                    end
                end
            end
            ok=false;
        end

        function ok=hasAllStrings(this, key, strs)
            [~, hs]=this.get(key);
            if hs.size>0
                N=length(strs);
                for i=1:N
                    if ~hs.contains(java.lang.String(strs{i}))
                        ok=false;
                        return;
                    end
                end
                ok=true;
            else
                ok=false;
            end
        end

        function all=values(this, ascending)
            if nargin<2
                ascending=true;
            end
            if ascending
                it=this.map.keySet.iterator;
            else
                it=this.map.descendingKeySet.iterator;
            end
            n=this.size;
            all=cell(1, n);
            i=1;
            while it.hasNext
                k=it.next;
                if ischar(k) && length(k)==1
                    k=java.lang.String(k);
                end
                l=this.map.get(k);
                it2=l.iterator;
                while it2.hasNext
                    all{i}=it2.next;
                    i=i+1;
                end
            end
        end

        function all=valueList(this, ascending)
            if nargin<2
                ascending=true;
            end
            if ascending
                it=this.map.keySet.iterator;
            else
                it=this.map.descendingKeySet.iterator;
            end
            all=java.util.ArrayList;
            while it.hasNext
                k=it.next;
                if ischar(k) && length(k)==1
                    k=java.lang.String(k);
                end
                all.addAll(this.map.get(k));
            end
        end

        
        function out=remove(this, k, v)
            out=[];
            if ischar(k) && length(k)==1
                k=java.lang.String(k);
            end
            if this.map.containsKey(k)
                l=this.map.get(k); 
                if ischar(v) && length(v)==1
                    v=java.lang.String(v);
                end
                out=l.remove(v);
            end
        end
        
        function renameKey(this, k, newKey)
            if ischar(k) && length(k)==1
                k=java.lang.String(k);
            end
            l=this.map.remove(k);
            if ~isempty(l)
                k=newKey;
                if ischar(k) && length(k)==1
                    k=java.lang.String(k);
                end
                this.map.put(k,l);
            end
        end
        
    end
end