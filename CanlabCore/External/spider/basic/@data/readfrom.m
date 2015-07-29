function dret = readfrom(d,format,filename)

switch(lower(format))   
case {'arff'}
    dret=readarff(filename);

case {'libsvm'}
    X=[];
    Y=[];
    
    fid=fopen(filename,'rt');
    index=1;
    
    while 1
        tline = fgetl(fid);
        
        if ~ischar(tline), break, end
        % read label
        [token,rem] = strtok(tline);
        Y(index,1)=str2num(token);
        [token,rem] = strtok(rem);
        
        while(~isempty(rem))
            a=sscanf(token,'%d:%e');
            X(index,a(1))=a(2);
            [token,rem] = strtok(rem);            
        end 
        index=index+1;
        if (mod(index,50)==0)
            fprintf('Line: %d\n',index)
        end
    end
    
    fclose(fid);
    dret=data(X,Y);
otherwise 
    warning(['Unknown format ',format])
end


function d=readarff(fname)
% Matlab Arff reader for the Spider package
% Author: gb   gb@tuebingen.mpg.de
% Date: 20.2.2006


f=fopen(fname);
if f>0
    dataset_name='';

    X=[];

    dataset_structured=0;
    attributs={};

    %     try
    state='scanningforrelation';

    while 1
        tline = fgetl(f);


        if ~ischar(tline), break, end
        % parseline
        switch state
            case 'scanningforrelation'
                if( ~isempty(strfind(lower(tline),'@relation')))
                    [tok,tline]=strtok(tline);
                    dataset_name=tline;
                    state='relation';
                    fprintf('Reading in %s\n', dataset_name)
                    %                 disp( tline)
                end
            case 'relation'
                if( ~isempty(strfind(lower(tline),'@attribute')))
                    [tok,tline]=strtok(tline);
                    atr=[];
                    % ============================================
                    % get attribute name
                    % ============================================


                    tline=strtrim(tline);
                    if(tline(1)=='''')
                        i=find(tline(2:end)=='''');
                        attrib_name=tline(2:i);
                        tline=tline(i+2:end);
                    else
                        [tok,tline]=strtok(tline);
                        attrib_name=tok;
                    end

                    atr.name=attrib_name;
                    tline=strtrim(tline);
                    % ============================================
                    % get attribute type(s)
                    % ============================================

                    %check if numeric
                    if( ~isempty(strfind(lower(tline), 'real')) |...
                            ~isempty(strfind(lower(tline), 'integer')))
                        atr.type='numeric';
                    elseif( ~isempty(strfind(lower(tline), 'string')) )
                        atr.type='string';
                        dataset_structured=1;
                    elseif( tline(1)=='{')
                        atr.type='nominal';
                        [tok,tline]=consume(tline,'{');

                        atr.list={};
                        while (~isempty(tline))
                            [tok,tline2]=consume(tline,',');
                            if(isempty(tline2))
                                [tok,tline]=consume(tline,'}');
                                tline=strtrim(tline);
                            else
                                tline=strtrim(tline2);
                            end

                            tok=tok(1:end-1);

                            atr.list{length(atr.list)+1}=strtrim(tok);

                        end
                    end

                    attributs{length(attributs)+1}=atr;

                elseif( ~isempty(strfind(lower(tline),'@data')))
                    state='readingdata';
                end

            case 'readingdata'
                % we parse tline dependening on attributs
                tline=strtrim(tline);
                if(tline(1)~='%')


                    if dataset_structured==1
                        x={};
                    else
                        x=[];
                    end
                    % is this a sparse file?
                    if(tline(1)=='{')
                        % parse sparse
                        error('Sorry, sparse arff files are not supported!')
                    else

                        for k=1:length(attributs)
                            switch attributs{k}.type
                                case 'numeric'


                                    [tok,tline]=readnext(tline);
                                    if(strfind(tok,'?'))
                                        x=addnan(x,k);
                                    else
                                        x=add(x,k,sscanf(tok,'%f'));
                                    end

                                case 'string'
                                    [tok,tline]=readnext(tline);
                                    if(tok(end)==',')
                                        tok=tok(1:end-1);
                                    end


                                case 'nominal'
                                    [tok,tline]=readnext(tline);

                                    if(tok(end)==',')
                                        tok=tok(1:end-1);
                                    end

                                    if(strfind(tok,'?'))
                                        x=addnan(x,k);
                                    else

                                        pos=get_in_list(attributs{k}.list, strtrim(tok));

                                        x=add(x,k,pos);
                                    end
                            end
                        end

                        X=[X;x];

                    end
                end

        end


    end
    %
    %     catch
    %         lasterr
    %     end


    fclose(f);
else
    warning(sprintf('Cant find file: %s\n',fname));
end

d=data(X);


function [tok,tline]=readnext(tline)
[tok,tline2]=consume(tline,',');
if( isempty(tok))
    tok=tline;
    tline =[];
else
    tline=tline2;
end



function [tok,tline]=consume(tline,token)
i=strfind(tline,token);
if isempty(i)
    tline=[];
    tok=[];
else
    tok=tline(1:i);
    tline=tline(i+1:end);
end


function  pos=get_in_list(alist, tok)
pos=[];
for k=1:length(alist)
    if(strcmp(alist{k},tok)==1)
        pos=k; break;
    end
end


function x=add(x,k,pos)
if(iscell(x))
    x{k}=pos;
else
    x(k)=pos;
end

function x=addnan(x,k)
if(iscell(x))
    x{k}=NaN;
else
    x(k)=NaN;
end


