%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef MultiProps < handle
    
    properties(SetAccess=private)
        instances;
        N;
    end
    
    methods
        function this=MultiProps(varargin)
            this.N=length(varargin);
            this.instances=varargin;
        end
        
        function isTrue=is(this, name, defaultIsTrue)
            isTrue=defaultIsTrue;
            for i=1:this.N
                isTrue=this.instances{i}.is(name, isTrue);
            end
        end

        function num=getNumeric(this, name, defaultNumber)
            num=defaultNumber;
            for i=1:this.N
                num=this.instances{i}.getNumeric(name, num);
            end
        end
        
        function num=getNumbers(this, name)
            num=[];
            for i=1:this.N
                temp=this.instances{i}.get(name, []);
                if ~isempty(temp)
                    num=temp;
                end
            end
            if isempty(num)
                num=[];
            else
                num=str2num(num);
            end
        end
        
        function ok=has(this, name)
            for i=1:this.N
                ok=this.instances{i}.has(name);
                if ok
                    return;
                end
            end
        end
    end
    
    methods(Static)
        function n=Name(sampleNode, gateNode, name)
            p1='';
            p2='';
            if ~isempty(sampleNode)
                p1=[sampleNode '.' ];
            end
            if ~isempty(gateNode)
                p2=[gateNode '.'];
            end
            n=[p1 p2 name]; 
        end
    end
    
    methods
        function ok=hasNode(this, sampleNode, gateNode, name)
            if this.N>1
                name=MultiProps.Name(sampleNode, gateNode, name);
                ok=this.instances{2}.has(name);
            else
                ok=false;
            end
        end
        
        function num=getNodeNumeric(this, sampleNode,gateNode, name, defaultNumber)
            num=defaultNumber;
            num=this.instances{1}.getNumeric(name, num);
            if this.N>1
                name=MultiProps.Name(sampleNode, gateNode, name);
                num=this.instances{2}.getNumeric(name, num);
            end
        end

        function setNodeNumeric(this, sampleNode, gateNode, name, value)
            num=num2str(value);
            this.instances{1}.set(name, num);
            if this.N>1
                name=MultiProps.Name(sampleNode, gateNode, name);
                this.instances{2}.set(name, num);
            end
        end
        
        function value=get(this, name, defaultValue)
            value=defaultValue;
            for i=1:this.N
                value=this.instances{i}.get(name, value);
            end
        end
        
        function value=remove(this, name)
            for i=1:this.N
                value=this.instances{i}.remove(name);
            end
        end
        

        function setBoolean(this, name, isTrue, defaultToo)
            if isTrue
                isTrue='true';
            else
                isTrue='false';
            end
            if nargin<4 || defaultToo
                for i=1:this.N
                    this.instances{i}.set(name, isTrue);
                end
            else
                this.instances{end}.set(name, isTrue);
            end
        end
        
        function values=getAll(this, name)
            for i=1:this.N
                values{i}=this.instances{i}.get(name);
            end
        end
        
        function setAll(this, name, values)
            for i=1:this.N
                if isempty(values{i})
                    this.instances{i}.remove(name);
                else
                    this.instances{i}.set(name, values{i});
                end
            end
        end
        
        function setNumeric(this, name, num, defaultToo)
            num=num2str(num);
            if nargin<4 || defaultToo
                for i=1:this.N
                    this.instances{i}.set(name, num);
                end
            else
                this.instances{end}.set(name, num);
            end
        end
        
        function set(this, name, value, defaultToo)
            if nargin<4 || defaultToo
                for i=1:this.N
                    this.instances{i}.set(name, value);
                end
            else
                this.instances{end}.set(name, value);
            end
        end
    end
end