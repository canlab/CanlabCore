%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

% this wrapper around MatLab internal map is meant for 
% polymorphic function signature compatability with other
% classdef and JAVA class in CytoGenie AutoGate
classdef Map < handle

    properties(Constant)
        PROP_ID='nextId';    
    end
    
    properties
        propertyFile;
    end
    
    properties(SetAccess=private)
        map;
    end
    
    methods(Static)
        function SetOnOffBoolean(props, prop, value)
            if strcmp('on', value)
                props.set(prop, 'true');
            else
                props.set(prop, 'false');
            end
        end

        function SetBoolean(props, prop, value)
            if value
                props.set(prop, 'true');
            else
                props.set(prop, 'false');
            end
        end

        
        function this=Global(closeIfTrueOrMap)
            persistent singleton;
            if nargin>0 && islogical(closeIfTrueOrMap) && closeIfTrueOrMap
                clear singleton;
                singleton=[];
                disp('Resetting global Map');
                this=[];
            else
                if nargin==1
                    try
                        priorMap=closeIfTrueOrMap;
                        %test method compatibility
                        priorMap.size
                        prop='IsBasicMap';
                        if ~priorMap.has(prop)
                            priorMap.set(prop, 'false')
                            priorMap.get(prop, 'true')
                            priorMap.remove(prop);
                        else
                            priorMap.get(prop, 'true')
                        end
                        singleton=priorMap;
                        this=singleton;
                        return;
                    catch ex
                        ex.getReport
                    end
                end
                if isempty(singleton) 
                    singleton=Map;
                    singleton.highDef=Gui.hasPcHighDefinitionAProblem(2000, 2500, false);
                    Map.SetHighDef(singleton, singleton.highDef);
                end     
                this=singleton;
            end
        end
        
        function path=Path
            path=Map.Global.contentFolder;
        end
        

        
        function nums=GetNumbers(props, name)
            nums=props.get(name, []);
            if ~isempty(nums)
                nums=str2num(nums);
            else
                nums=[];
            end
        end
        
        function SetStruct(this, id, structure)
            names=fieldnames(structure)';
            N=length(names);
            for i=1:N
                name=names{i};
                value=getfield(structure, name);
                this.set([id '.' name], value);
            end
        end
    end
    
    methods
        function this=Map(keysOrFileName, values)
            if nargin<2
                values={};
                if nargin<1
                    keysOrFileName=[];
                end
            end
            if ~isempty(values)
                this.map=containers.Map(keysOrFileName, values);
            elseif ~isempty(keysOrFileName) && ischar(keysOrFileName)
                this.propertyFile=keysOrFileName;
                this.load;
            else
                this.map=containers.Map;
            end
        end
        
        function cnt=size(this)
            cnt=this.map.length;
        end

        function reset(this)
            this.clear;
        end
        
        function clear(this)
            remove(this.map, this.map.keys);
        end

        function priorValue=remove(this, name)
            priorValue=this.get(name);
            remove(this.map, name);
        end
        
        function value=getIdentified(this, id, name)
            value=this.get([id '.' name]);
        end
        
        function set(this, name, value)
            this.map(name)=value;
        end
        
        function prior=put(this, name, value)
            prior=this.get(name);
            this.map(name)=value;
        end
        
        %for value char compatibility with outside map classes
        function setBoolean(this, name, isTrue)
            if isTrue
                isTrue='true';
            else
                isTrue='false';
            end
            this.map(name)=isTrue;
        end
        
        % for value char compatibility with outside map classes
        function setNumeric(this, name, num)
            this.map(name)=num2str(num);
        end

        function isTrue=is(this, name, defaultIsTrue)
            if this.map.isKey(name)
                value=this.map(name);
                if ischar(value)
                    value=strcmpi('yes', value) || strcmpi('true', value);
                end 
                isTrue=value;
            else
                if nargin<3
                    isTrue=false;
                else
                    isTrue=defaultIsTrue;
                end
            end
        end
        
        function value=get(this, name, defaultValue)
            if this.map.isKey(name)
                value=this.map(name);
            else
                if nargin<3
                    value=[];
                else
                    value=defaultValue;
                end
            end
        end
        
        
        function value=getNumeric(this, property, defaultValue)
            value=this.get(property);
            if isempty(value) 
                if nargin>2
                    value=defaultValue;
                else
                    value=0;
                end
            else
                if isa(value,'char')
                    value=str2double(value);
                    if isnan(value)
                        if nargin>2
                            value=defaultValue;
                        else
                            value=0;
                        end
                    end
                elseif ~isnumeric(value) || any(isnan(value)) 
                    if nargin>2
                        value=defaultValue;
                    else
                        value=0;
                    end
                end
            end
        end
        
        function ok=has(this, name)
            ok=this.map.isKey(name);
        end
        
        function ok=containsKey(this, name)
            ok=this.map.isKey(name);
        end
        
        function save(this, propertyFile)
            if nargin>1
                this.propertyFile=propertyFile;
            end
            if ~isempty(this.propertyFile)
                basicMap=this;
                save(this.propertyFile, 'basicMap');
            end
        end
         
        function load(this, propertyFile)
            if nargin>1
                this.propertyFile=propertyFile;
            end
            if ~isempty(this.propertyFile) && exist(this.propertyFile, 'file')
                load(this.propertyFile, 'basicMap');
                this.map=basicMap.map;
            end
        end
         
        
        
        function k=keys(this)
            k=this.map.keys;
        end
        
        function this=addAll(this, that, keys)
            if nargin<3
                it=that.map.keySet.iterator;
                while(it.hasNext)
                    name=it.next;
                    this.set(name, that.props{that.map.get(name)});
                end
            else
                N=length(keys);
                for i=1:N
                    if that.has(keys{i})
                        this.set(keys{i}, that.get(keys{i}));
                    end
                end
            end
        end
        
        function idx=indexOf(this, name, value)
            idx=-1;
            N=this.countIfMultipleElse0(name);
            if N==0 && this.has(name)
                v=this.get(name);
                if strcmp(v, value)
                    idx=1;
                end
            else
                for i=1:N
                    v=this.get([name '.' num2str(i)]);
                    if strcmp(v, value)
                        idx=i;
                        break;
                    end
                end
            end            
        end
        
        function wasMissing=addIfMissing(this, name, value)
            idx=this.indexOf(name,value);
            if idx<=0
                wasMissing=true;
                this.add(name,value);
            else
                wasMissing=false;
            end
        end
        
        function wasMissing=putIfMising(this, name, value)
            if ~this.has(name)
                wasMissing=true;
                this.put(name,value);
            else
                wasMissing=false;
            end
        end
        
        function debug(this)
            k=this.map.keys;
            N=length(k);
            for i=1:N
                name=k{i};
                v=this.get(name);
                if ischar(v)
                    fprintf('%s=%s\n', name, v);
                elseif isnumeric(v)
                    fprintf('%s=%d\n', name, v);
                else
                    ss=class(v);
                    fprintf('%s=[%s object]\n', name,  ss);
                end
            end
        end

        function add(this, name, value)
            name2=this.increment(name);
            this.set(name2, value);
        end
        
        function setAll(this, name, values)
            this.remove(name);
            N=length(values);
            for i=1:N
                this.add(name,values{i});
            end
        end
        
        function c=getAll(this, name)
            N=this.countIfMultipleElse0(name);
            if N==0 && this.has(name)
                c{1}=this.get(name);
            else
                c=cell(1,N);
                for i=1:N
                    c{i}=this.get([name '.' num2str(i)]);
                end
            end            
        end
        
        function name2=increment(this, name)
            N=this.countIfMultipleElse0(name)+1;
            if N==1 && this.has(name)
                N=N+1;
                value=this.get(name);
                this.remove(name);
                this.set([name '.1'], value);
            end
            name2=[name '.' num2str(N)];
            this.set([name '.count'], num2str(N));
        end

        function N=countIfMultipleElse0(this, name)
            N=0;
            if ischar(name)
                num=this.get([name '.count']);
                if ~isempty(num)
                    N=str2num(num);
                end
            end
        end
        
        function id=newId(this)
            id=num2str(this.getNumeric(Map.PROP_ID, 0)+1);
            this.set(Map.PROP_ID, id);
            this.set([id '.createdWhen'], char(datetime));
            this.save;
        end
    end
end