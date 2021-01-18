classdef JavaProperties < handle
    properties(SetAccess=private)
        p;
        fileName;
    end
    
    methods
        function this=JavaProperties(fileName)
            if nargin<1
                this.p=java.util.Properties;
            else
                this.p=File.ReadProperties(fileName, true);
                this.fileName=fileName;
            end
        end
        
        function load(this, fileName)
            this.p=File.ReadProperties(fileName);
        end
        
        function save(this, fileName)
            if nargin<2
                fileName=this.fileName;
            end
            File.SaveProperties2(fileName, this.p);
        end
        
        function cnt=size(this)
            cnt=this.p.size;
        end

        function reset(this)
            this.p.clear;
        end
        
        function clear(this)
            this.p.clear;
        end

        function priorValue=remove(this, name)
            priorValue=this.get(name);
            this.p.remove(java.lang.String(name));
        end
        
        function value=getIdentified(this, id, name)
            value=this.get([id '.' name]);
        end
        
        function priorValue=set(this, name, value)
            priorValue=this.p.setProperty(java.lang.String(name), java.lang.String(value));
        end
        
        function prior=put(this, name, value)
            prior=this.set(name, value);
        end
        
        %for value char compatibility with outside map classes
        function prior=setBoolean(this, name, isTrue)
            if isTrue
                isTrue='true';
            else
                isTrue='false';
            end
            prior=this.set(name, isTrue);
        end
        
        % for value char compatibility with outside map classes
        function prior=setNumeric(this, name, num)
            prior=this.set(name, num2str(num));
        end

        function isTrue=is(this, name, defaultIsTrue)
            if this.containsKey(name)
                value=char(this.get(name));
                isTrue=strcmpi('yes', value) || strcmpi('true', value);
            else
                if nargin<3
                    isTrue=false;
                else
                    isTrue=defaultIsTrue;
                end
            end
        end
        
        function value=get(this, name, defaultValue)
            if this.containsKey(name)
                value=char(this.p.getProperty(java.lang.String(name)));        
            else
                if nargin<3
                    value=[];
                else
                    value=defaultValue;
                end
            end
        end
        
        
        function value=getNumeric(this, name, defaultValue)
            value=this.get(name);
            if isempty(value) 
                if nargin>2
                    value=defaultValue;
                else
                    value=0;
                end
            else
                value=str2double(value);
                if isnan(value)
                    if nargin>2
                        value=defaultValue;
                    else
                        value=0;
                    end
                end
            end
        end
        
        function ok=has(this, name)
            ok=this.containsKey(name);
        end
        
        function ok=containsKey(this, name)
            ok=this.p.containsKey(java.lang.String(name));
        end
        
        
        function k=keys(this)
            k=CellBasics.Java(  this.p.keySet );
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
        
        function wasMissing=putIfMissing(this, name, value)
            if ~this.has(name)
                wasMissing=true;
                this.put(name,value);
            else
                wasMissing=false;
            end
        end
        
        function debug(this)
            k=this.keys;
            N=length(k);
            for i=1:N
                name=k{i};
                v=this.get(name);
                fprintf('%s=%s\n', name, v);
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