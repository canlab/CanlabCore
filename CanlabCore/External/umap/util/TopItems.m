%
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

classdef TopItems <handle
    properties(SetAccess=private)
        tm;
        top=10;
        highestWins=true;
    end
    methods
        function this=TopItems(topNumber, highestWins)
            if nargin>1
                this.highestWins=highestWins;
            end
            this.tm=java.util.TreeMap;
            if nargin>0
                this.top=topNumber;
            end
        end
        
        function N=size(this) 
            N=this.tm.size;
        end
        
        function add(this, item, idx)
            this.tm.put(item, idx);
            if this.tm.size>this.top
                this.tm.remove(this.worst);
            end
        end
        
        function [item, idx]=worst(this)
            if this.highestWins
                item=this.tm.firstKey;
            else
                item=this.tm.lastKey;
            end
            if nargout>1
                idx=this.tm.get(item);
            end
        end
        
        function [item, idx]=best(this)
            if this.highestWins
                item=this.tm.lastKey;
            else
                item=this.tm.firstKey;
            end
            if nargout>1
                idx=this.tm.get(item);
            end
        end 
        
        function [idxs, items]=all(this)
            N=this.tm.size;
            idxs=zeros(1, N);
            items=cell(1,N);
            i=1;
            if this.highestWins
                it=this.tm.descendingKeySet.iterator;
            else
                it=this.tm.keySet.iterator;
            end
            while it.hasNext
                items{i}=it.next;
                idxs(i)=this.tm.get(items{i});
                i=i+1;
            end
        end
    end
end