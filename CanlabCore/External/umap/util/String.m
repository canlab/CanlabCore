%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%


classdef String
    properties(SetAccess=private)
        js;
    end
    
    properties(Constant)
        SHELL='[\\\|&\(\)< >'':\`\*;"]';
    end
    methods
        
        function string=String(str)
            string.js=java.lang.String(str);
        end
        
        function ok=endsWith(this, ending)
            ok=this.js.endsWith(ending);
        end
        
        function ok=startsWith(this, ending)
            ok=this.js.startsWith(ending);
        end
        
        function ok=contains(this, str)
            ok=this.js.indexOf(str)>=0;
        end
        
        function i=indexOf(this, str)
            i=this.js.indexOf(str)+1;
        end
        
        function i=lastIndexOf(this, str)
            i=this.js.lastIndexOf(str)+1;
        end
        
        function str=capitalize(this)
            str=char(this.js);
            if length(str)>0
                str(1,1)=upper(str(1,1));
            end
            
        end
        function str=uncapitalize(this)
            str=char(this.js);
            if length(str)>0
                str(1,1)=lower(str(1,1));
            end
            
        end
        
        function out=subString(this, idx)
            out=char(this.js.substring(idx-1));
        end
        
        function out=subString2(this, idx1, idx2)
            out=char(this.js.substring(idx1-1, idx2-1));
        end
        
    end 
    
    
    methods(Static)
        function str=AddSuffix(str, suffix)
            N=length(str);
            if N>13
                if strcmpi('</html>', str(end-6:end))
                    str=[str(1:end-7) suffix '</html>'];
                    return;
                end
            end
            str=[str suffix];
        end
        function out=SubField(str, strus, idxs, fldOffset, ...
                fieldName, newFieldValues)
            N=length(strus);
            if N==0
                out=str;
            else
                out='';
                start=1;
                
                for i=1:N
                    stop=idxs(i)+fldOffset;
                    fld=getfield(strus(i), fieldName);
                    out=[out str(start:stop) newFieldValues{i} ];
                    start=stop+length(fld)+1;
                end
                if start<length(str)+1
                    out=[out str(start:end)];
                end
            end
        end
        
        function str=EncodeStrs(strs)
            N=length(strs);
            str=strs{1};
            for i=2:N
                str=[str ';' strs{i}];
            end
        end
        function strs=DecodeStrs(str)
            strs=strsplit(str, ';');
        end
        function ttl=EscapeHtmlTex(in)
            ttl=String.ToLaTex(char(...
                edu.stanford.facs.swing.Basics.RemoveXml(in)));
        end

        function in=RemoveTex(in)
            if any( strfind(in, '^{'))
                in=strrep(in, '^{', '');
                in=strrep(in, '}','');
            end
            in=strrep(strrep(in, '_', ''), '^', '');
        end

        function num=Rank(num)
            if num>3 && num<14
                prefix='th';
            else
                sNum=num2str(num);
                if sNum(end)=='1'
                    prefix='st';
                elseif sNum(end)=='2'
                    prefix='nd';
                elseif sNum(end)=='3'
                    prefix='rd';
                else
                    prefix='th';
                end
            end
            num=[num2str(num) prefix];
        end
        
        function str=Num(num, sigDigits, decDigits, delimiter)
            if nargin<3
                decDigits=0;
                if nargin<2
                    sigDigits=0;
                end
            end
            N=length(num);
            if N>1
                if nargin<4
                    delimiter=', ';
                end
            end
            str='';
            for i=1:N
                if isequal('k', decDigits)
                    nStr=String.encodeK(num(i));
                elseif decDigits>0
                    nStr=String.encodeNumber(num(i), decDigits);
                else
                    nStr=num2str(num(i));
                end
                str=[str String.Pad(nStr, sigDigits)];
                if i<N
                    str=[str delimiter];
                end
            end
        end

        function str=Num2Str(num, delimiterOrPadding)
            if nargin<2
                delimiterOrPadding='_';
            end
            isDelimiter=ischar(delimiterOrPadding);
            N=length(num);
            str='';
            for i=1:N
                str=[str num2str(num(i)) ];
                if isDelimiter
                    if i<N
                        str=[str delimiterOrPadding];
                    end
                end
            end
            if ~isDelimiter % is padding
                str=String.Pad(str, delimiterOrPadding);
            end
        end
        
        function str=Pad(str, size, ch)
            if nargin<3
                ch=' ';
            end
            pads=size-length(str);
            for i=1:pads
                str=[ch str];
            end
        end

        function str=PadHtml(str, size)            
            pads=size-length(str);
            for i=1:pads
                str=['&nbsp;' str];
            end
            str=['<html>' str '</html>'];
        end

        function str=PadEnding(str, size, ch)
            if nargin<3
                ch=' ';
            end
            for i=1:size
                str=[str ch];
            end
        end

        function txt=TimeEstimate(secs)
            if secs<3600
                mins=secs/60;
                mins2=floor(mins);
                if mins2<1 && secs<40
                    txt=String.Pluralize2('sec', secs);
                else
                    secs2=secs-(mins2*60);
                    if secs2>50 
                        strMins=String.Pluralize2('min', mins2+1);
                        txt=['about ' strMins];
                    elseif secs2<10
                        strMins=String.Pluralize2('min', mins2);
                        txt=['about ' strMins];
                    elseif secs2>30
                        strMins=String.Pluralize2('min', mins2+1);
                        txt=['less than ' strMins];
                    else
                        txt=String.MinutesSeconds(secs);
                    end
                end
            else
                hours=secs/3600;
                hours2=floor(hours);
                strHours=String.Pluralize2('hour', hours2);
                secs=secs-(hours2*3600);
                if secs>5
                    if secs<60
                        txt=[strHours ', ' num2str(secs) ' secs'];
                    else
                        strMinutes=String.encodeRounded(secs/60,1);
                        txt=[strHours ', ' strMinutes ' mins'];
                    end
                else
                    txt=strHours;
                end
            end
        end


        function txt=TimeReport(secs, lastSecs, isEstimate, ...
                reportingInterval)
            if nargin<4
                reportingInterval=6;
                if nargin<3
                    isEstimate=false;
                end
            end
            dif=abs(secs-lastSecs);
            %fprintf('%s between %s %s %d', num2str(dif), num2str(lastSecs), num2str(secs), reportingInterval);
            if lastSecs==0 || dif>reportingInterval
                if isEstimate
                    txt=String.TimeEstimate(secs);
                else
                    txt=String.HoursMinutesSeconds(secs);
                end
                %fprintf(' reporting\n');
            else
                txt=[];
                %fprintf(' NO REPORT\n');
            end
        end
        
        function txt=HoursMinutesSeconds(secs)
            if secs<3600
                txt=String.MinutesSeconds(secs);
            else
                hours=secs/3600;
                hours2=floor(hours);
                strHours=String.Pluralize2('hour', hours2);
                secs=secs-(hours2*3600);
                if secs>5
                    if secs<60
                        txt=[strHours ', ' num2str(secs) ' secs'];
                    else
                        strMinutes=String.encodeRounded(secs/60,1);
                        txt=[strHours ', ' strMinutes ' mins'];
                    end
                else
                    txt=strHours;
                end
            end
        end
        function secs=ToSecs(strMinSecs)
            s=string(strMinSecs);
            try
                secs=str2double(s.extractBefore(' min'))*60;
                secs=secs+str2double(s.extractBetween(',', ' sec'));
            catch 
                try
                    secs=str2double(s.extractBefore(' sec'));
                catch
                  secs=0;
                end
            end
        end
        function txtSecs=MinutesSeconds(secs)
            if secs<60
                if secs>9
                    txtSecs=[String.encodeRounded(secs,1) ' secs'];
                else
                    txtSecs=[String.encodeRounded(secs,2) ' secs'];
                end
            else
                mins=floor(secs/60);
                secs=mod(ceil(secs), 60);
                if secs<1
                    txtSecs=[String.Pluralize2('min', mins, 'mins')];
                else
                    mins=String.Pluralize2('min', mins, 'mins');
                    secs=String.Pluralize2('sec', secs);
                    txtSecs=[mins ', ' num2str(secs)];
                end
            end
        end
        
        function ok=IsEmpty(str)
            ok=true;
            if ~isempty(str)
                if ~isempty(strtrim(str))
                    ok=false;
                end
            end
        end
        function c=Low2High(low, high, step, decimalPlaces, suffix)
            c={};
            n=low:step:high;
            N=length(n);
            for i=1:N
                s=String.encodeRounded(n(i), decimalPlaces, true);
                c{end+1}=[s suffix];
            end
        end
        function letter=toLetter(number)
            alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
            num=mod(number,26);
            if (num==0)
                num=26;
            end
            letter=alphabet(num);            
            if number>26
                letter=strcat(letter, num2str(ceil(number/26)));
            end
        end
        
        function out=GetSuffix(string, lastDelimiter)
            li=String.LastIndexOf(string, lastDelimiter);
            if li>0
                out=String.SubString(string, li+1);
            else
                out='';
            end
        end
        
        function out=PruneSuffix(string, suffix)
            if isequal(string(end-length(suffix)+1:end), suffix)
                out=string(1:end-length(suffix));
            else
                out=string;
            end
        end
        
        function out=PruneSuffixSlow(string, suffix)
            str=String(string);
            if str.endsWith(suffix)
                N2=length(suffix);
                N1=length(string);
                out=str.subString2(1, N1-N2+1);
            else
                out=string;
            end
        end
        
        function out=PrunePrefix(string, prefix)
            str=String(string);
            if str.startsWith(prefix)
                N2=length(prefix);
                N1=length(string);
                out=str.subString2(N2+1, N1+1);
            else
                out=string;
            end
        end
        
        function str=Pluralize2(singularItem, N, pluralItem)            
            if nargin>2
                str=String.Pluralize(singularItem, N, pluralItem);
            else
                str=String.Pluralize(singularItem, N, '');
            end
            str=sprintf('%s %s', String.encodeInteger(N), str);
        end
        
        function str=Pluralize(singularItem, N, pluralItem)
            if N>1 || N==0
                if nargin>2 && ~isempty(pluralItem)
                    str=pluralItem;
                else
                    str=[singularItem 's'];
                end
            else
                str=singularItem;
            end
        end

        function str=PluralCntItem_(singularItem, N, pluralItem)            
            str=sprintf('%d %s', N, String.PluralItem_(singularItem, N, pluralItem));
        end
        
        function str=PluralItem_(singularItem, N, pluralItem)
            if N>1 || N==0
                str=pluralItem;
            else
                str=singularItem;
            end
        end

        function str=PluralCntItem(singularItem, N)  
            str=sprintf('%d %s', N, String.PluralItem(singularItem, N));
        end
        
        function str=PluralItem(singularItem, N)
            if N>1 || N==0
                str=[singularItem 's'];
            else
                str=singularItem;
            end
        end

        function str=NoBulletIfZero(singularForm, N, pluralForm)
            if nargin<3
                pluralForm='';
            end
            if N~=0
                str=String.Pluralize2(singularForm, N, pluralForm);
                str=['<li>' str]; 
            else
                str='';
            end
        end
        
        function str=Capitalize(str)
            str=String(str);
            str=str.capitalize();
        end
        
        function ok=EndsWith(str, ending)
            ok=java.lang.String(str).endsWith(ending);
        end
        
        function ok=StartsWith(str, beginning)
            ok=java.lang.String(str).startsWith(beginning);
        end
        
        
        function ok=EndsWithI(str, ending)
            ok=java.lang.String(str).toLowerCase.endsWith(lower(ending));
        end
        
        function ok=StartsWithI(str, beginning)
            ok=java.lang.String(str).toLowerCase.startsWith(lower(beginning));
        end
        
        function i=IndexOf(str, test)
            i=java.lang.String(str).indexOf(test)+1;
        end
        
        function i=LastIndexOf(str, test)
            i=java.lang.String(str).lastIndexOf(test)+1;
        end
        
        function out=SubString(str, idx)
            out=String(str).subString(idx);
        end
        
        function ok=Contains(str, search)
            ok=~isempty(strfind(str, search));
        end
        
        function ok=ContainsI(str, search)
            ok=~isempty(strfind(lower(str), lower(search)));
        end
        
        function out=SubString2(str, idx1, idx2)
            out=String(str).subString2(idx1, idx2);
        end
        
        function [ prefix ]=Prefix( string, endOfPrefixIfFound )
            %UNTITLED Summary of this function goes here
            %   Detailed explanation goes here
            
            a=strfind(string, endOfPrefixIfFound);
            if isempty(a)
                prefix='';
            else
                lastIdx=a(1,size(a,2));
                prefix=string(1, [1: lastIdx-1]);
            end
        end

        function str=Uncapitalize(str)
            if ~isempty(str)
                str(1,1)=lower(str(1,1));
            end
        end

        function out= ToHtml(in)
            out=strrep(in, '&', '&amp;');
            out=strrep(out, '<', '&lt;');            
            out=strrep(out, '>', '&gt;');            
        end

        function out=ToHtmlSupFromTex(out)
            if contains(out, '^{')
                app=BasicMap.Global;
                out=strrep(out, '^{', app.supStart);
                out=strrep(out, '}', app.supEnd);
            end
            
        end     
        
        function output=ToFile(input, escapeHtmlTex)
            if nargin>1 && escapeHtmlTex
                input=String.EscapeHtmlTex(input);
            end
            if ~isempty(regexp(input,  String.SHELL))
                output=regexprep(input,  String.SHELL, '_');
            else
                output=input;
            end
            if ~isempty(regexp(output, ','))
                output=regexprep(output, ',', '');
            end
            if ~isempty(regexp(output, '/'))
                output=regexprep(output, '/', '_');
            end
        end

        function output=ToLaTex(input)
            output=strrep(strrep(strrep(strrep(...
                strrep(input, '\', '\\'), '_', '\_'), '^', '\^'), ...
                '{', '\{'), '}', '\}');
        end
        
        function output=ToSystem(input)
            if ~isempty(regexp(input, String.SHELL))
               if ispc
                   output=['"' input '"'];
               else
                   output=regexprep(input, String.SHELL, '\\$0');
               end
            else
                output=input;
            end
        end

        function out=RemoveXml(in)
            if isempty(in)
                out=[];
            else
                if length(in)>6 && strcmp( in(1:6), '<html>')
                    out=regexprep(in, '<[^>]*', '');
                    out=regexprep(out, '>*', '');
                    out=strrep(out, '&lt;', '');
                    out=strrep(out, '&gt;', '');
                else
                    out=in;
                end
            end
        end
        
        function str=encodeNum(num)
            str=String.encodeRounded(num, 3, true);
        end
        
        function str=encodeBank(num, prfx, sfx)
            str=String.encodeRounded(num, 2, true);
            if length(str)>3
                if isequal('0.', str(1:2)) 
                    str=str(2:end);
                elseif isequal('< 0.', str(1:4))
                    str=['<' str(4:end)];
                end
            end
            if nargin>1
                str=[prfx str sfx];
            end
        end
        
        function [str]=encodeNumbers(num)
            str = arrayfun(@(x) String.encodeNumber(x) , num, 'UniformOutput', false) ;
        end
        
        function [str]=encodeNumber(num, decimalPlaces)
            if nargin<2
                decimalPlaces=2;
            end
            div=10^decimalPlaces;
            num=floor(num*div)/div;
            isNeg=0>num;
            num=abs(num);
            str=['%18.' num2str(decimalPlaces) 'f'];
            str=strtrim(sprintf(str, num));
            k=find(str == '.', 1);
            isInt=false;
            if(isempty(k))
                isInt=true;
                tail='.';
                if isinf(decimalPlaces)
                    decimalPlaces=0;
                end
                for i=1:decimalPlaces
                    tail=strcat(tail, '0');
                end
                str=[str,tail];
            end
            
            
            FIN = min(length(str),find(str == '.')-1);
            for i = FIN-2:-3:2
                str(i+1:end+1) = str(i:end);
                str(i) = ',';
            end
            if isInt
                str= String.Prefix(str, tail);
            end
            if isNeg
                str=['-' str];
            end
        end
        
        function num=encodeInteger(num)
            num=String.encodeNumber(floor(num));
            num=String.PruneSuffix(num,'.00');
        end

        function s=encodeCntAndPercent(numerator, denominator)
            s={ [ String.encodeK(numerator) '/' ...
                String.encodeK(denominator)], ...
                [String.encodePercent(numerator, denominator)]};
        end
        
        function out=encodePadHtml(in, sig, dec, suffix)
            if nargin<4
                suffix='';
            end
            lenSuffix=length(suffix);
            w=sig+dec+1;
            if isnan(in)
                out='<html> N/A';
                w=w-3;
                for i=1:w
                    out=[out '&nbsp;'];
                end
                out=[out '</html>'];
            else
                s=String.encodeRounded(in, dec);
                if s(1) == '<'
                    w=w-(length(s)-1);
                    out='<html>&lt;';
                    for i=1:w
                        out=[out '&nbsp;'];
                    end
                    out=[out s(3:end) suffix '</html>'];
                elseif dec==0
                    out=String.PadHtml([s suffix], w-1+lenSuffix);
                else
                    out=String.PadHtml([s suffix], w+lenSuffix);
                end
            end
        end
        
        function out=encodePercent(numerator, denominator, decimals)
            if nargin<3
                decimals=2;
                if nargin<2
                    denominator=1;
                end
            end
            in=numerator/denominator*100;
            if in==0
                out='0%';
            elseif in<1
                if in<10^-decimals && in <10^-(decimals+1)
                    out=num2str(10^-(decimals+1));
                    out(1)='<';
                    out=[out '%'];
                else
                    out=[String.encodeRounded(in,decimals,true) '%'];
                end
            elseif in<10
                out=[String.encodeRounded(in,1,true) '%'];
            else
                if in<100 && in>98.49
                    if in>99.9
                        out='>99.9%';
                    else
                        out=[String.encodeRounded(in, 1) '%'];
                    end
                elseif isnan(in)
                    out='NaN';
                else
                    out=[String.encodeRounded(in, 0) '%'];
                end
            end
        end

        function num=encodeCount(num, pref)
            if abs(num)<1000
                num=String.encodeInteger(num);
            else
                if nargin==1
                    pref=BasicMap.Global.getNumeric('cntPref', 2);
                end
                if pref==0
                    num=String.encodeInteger(num);
                elseif pref==1
                    num=[String.encodeRounded(num/1000, 2, true) 'k'];
                elseif pref==3
                    num=String.encodeK(num);
                else
                    num=[String.encodeRounded(num/1000, 1, true) 'k'];
                end
            end
        end
        
        function num=encodeK(num, min, decimalPlaces)
            if nargin<3
                decimalPlaces=0;
            end
            if isempty(num) || isnan(num)
                num='n/a';
                return;
            end
            negative=num<0;
            if negative
                num=abs(num);
            end
            if nargin<2 || isempty(min)
                min=2^10;
            end
            if num<min
                num=String.encodeInteger(num);
            else
                if decimalPlaces==0
                    num=[String.encodeInteger(round(num/(2^10))) 'k'];
                else
                    num=[String.encodeRounded(num/(2^10), decimalPlaces) 'k'];
                end
            end
            if negative
                num=['-' num];
            end
        end
        
        function num=encodeMb(num, min, decimalPlaces)
            if nargin<3
                decimalPlaces=0;
            end
            if isempty(num) || isnan(num)
                num='n/a';
                return;
            end
            negative=num<0;
            if negative
                num=abs(num);
            end
            if nargin<2 || isempty(min)
                min=2^20;
            end
            if num<min
                num=String.encodeK(num, [], decimalPlaces);
            else
                if decimalPlaces==0
                    num=[String.encodeInteger(round(num/(2^20))) 'MB'];
                else
                    num=[String.encodeRounded(num/(2^20), decimalPlaces) 'MB'];
                end
            end
            if negative
                num=['-' num];
            end
        end
        
        function num=encodeBytes(num, decimalPlaces)
            if nargin<2
                decimalPlaces=1;
            end
            if isempty(num) || isnan(num)
                num='n/a';
                return;
            end
            negative=num<0;
            if negative
                num=abs(num);
            end
            if num<2^10
                num=[String.encodeInteger(num) ' bytes'];
            elseif num<2^20
                num=[String.encodeRounded(num/(2^10), decimalPlaces, true) ' KB'];
            else
                num=[String.encodeRounded(num/(2^20), decimalPlaces, true) ' MB'];
            end
            if negative
                num=['-' num];
            end
        end
        
        function num=encodeEMD(num)
            if ischar(num)
                num=str2double(num);
            end
            num=String.encodeRounded(num, 3);
        end
        
        function num=encodeVectorRounded(nums, decimalPlaces, delimiter)
            if nargin<3
                delimiter=' ';
            end
            N=length(nums);
            num=[];
            for i=1:N
                if i>1
                    num=[num delimiter String.encodeRounded(nums(i), decimalPlaces, true)];
                else
                    num=String.encodeRounded(nums(i), decimalPlaces, true);
                end
            end
            
        end
        function num=encodeRounded(num, decimalPlaces, zeroIfInteger, ...
                num2StrThreshold, lessThanSign)
            if num==0
                num='0';
                return
            end
            if nargin==1
                decimalPlaces=1;
            end
            if nargin>3 && ~isempty(num2StrThreshold) && ...
                    decimalPlaces>=num2StrThreshold
                num=num2str(num);
                return;
            end
            f=10^decimalPlaces;
            num=round(num*f);
            num=floor(num);
            num=num/f;
            if nargin>2 && zeroIfInteger
                if (num>=1 || num<=-1 )&& round(num)==num
                    num=String.encodeInteger(num);
                    return;
                end
            end
            if num==0 
                if nargin<5 || lessThanSign
                    num=['< ' String.encodeNumber(1/f,decimalPlaces)];
                else
                    num=String.encodeNumber(1/f,decimalPlaces);
                end
            else
                num=String.encodeNumber(num, decimalPlaces);
            end
        end
        function name = getNowName( prefix )
            %UNTITLED4 Summary of this function goes here
            %   Detailed explanation goes here
            t=today;
            n=now;
            name=[prefix '_' num2str( year(t), '%4d')  ];
            name = [name num2str( month(t), '%2d')];
            
            name = [name num2str( day(t), '%2d') '_' ];
            name = [name num2str( hour(n), '%2d') ];
            name = [name num2str( minute(n),'%2d') ];
            seconds=second(n);
            seconds=round(seconds);
            name = [name num2str( seconds,'%2d')];
        end
        

        function idx=Longest(strs)
            max=0;
            idx=0;
            N=length(strs);
            if N==1
                idx=1;
            else
                for i=1:N
                    if length(strs{i}) > max
                        max=length(strs{i});
                        idx=i;
                    end
                end
            end
        end
        
        function num=GetSuffixNum(str, lastDelim)
            str=String(str);
            i=str.lastIndexOf(lastDelim);
            if i>0
                sub=str.subString(i+1);
                num=str2double(sub);
                if isnan(num)
                    num=[];
                end
            else
                num=[];
            end
        end
        
        function xmlNode=GetXmlNode(str, nodeName)
            xmlNode='';
            if ischar(str)
                b=String.IndexOf(str, ['<' nodeName '>' ]);
                e=String.IndexOf(str, ['</' nodeName '>' ]);
                if b>-1 && e >b
                    xmlNode=String.SubString2(str, b+length(nodeName)+2, e);
                end
            end
        end

        function op=PlusOrMinus(percentage)
            if percentage<.09
                op='--';
            elseif percentage>.75
                op='++';
            elseif percentage<.25
                op='-';
            elseif percentage>.45
                op='+';
            else
                op='';
            end
        end
        
        function out=JavaArray(in)
            if ischar(in)
                out=in;
            else
                out='';
                if in.length>0
                    out=in(1);
                    for i=2:in.length
                        out=[out ', ' in(i)];
                    end
                end
            end
        end
        function c = linewrap(s, maxchars)
            %LINEWRAP Separate a single string into multiple strings
            %   C = LINEWRAP(S, MAXCHARS) separates a single string into multiple
            %   strings by separating the input string, S, on word breaks.  S must be a
            %   single-row char array. MAXCHARS is a nonnegative integer scalar
            %   specifying the maximum length of the broken string.  C is a cell array
            %   of strings.
            %
            %   C = LINEWRAP(S) is the same as C = LINEWRAP(S, 80).
            %
            %   Note: Words longer than MAXCHARS are not broken into separate lines.
            %   This means that C may contain strings longer than MAXCHARS.
            %
            %   This implementation was inspired a blog posting about a Java line
            %   wrapping function:
            %   http://joust.kano.net/weblog/archives/000060.html
            %   In particular, the regular expression used here is the one mentioned in
            %   Jeremy Stein's comment.
            %
            %   Example
            %       s = 'Image courtesy of Joe and Frank Hardy, MIT, 1993.'
            %       c = linewrap(s, 40)
            %
            %   See also TEXTWRAP.

            % Steven L. Eddins
            % $Revision: 1.19 $  $Date: 2020/11/28 23:16:31 $

            error(nargchk(1, 2, nargin));

            bad_s = ~ischar(s) || (ndims(s) > 2) || (size(s, 1) ~= 1);
            if bad_s
               error('S must be a single-row char array.');
            end

            if nargin < 2
               % Default value for second input argument.
               maxchars = 80;
            end

            % Trim leading and trailing whitespace.
            s = strtrim(s);

            % Form the desired regular expression from maxchars.
            exp = sprintf('(\\S\\S{%d,}|.{1,%d})(?:\\s+|$)', maxchars, maxchars);

            % Interpretation of regular expression (for maxchars = 80):
            % '(\\S\\S{80,}|.{1,80})(?:\\s+|$)'
            %
            % Match either a non-whitespace character followed by 80 or more
            % non-whitespace characters, OR any sequence of between 1 and 80
            % characters; all followed by either one or more whitespace characters OR
            % end-of-line.

            tokens = regexp(s, exp, 'tokens').';

            % Each element if the cell array tokens is single-element cell array 
            % containing a string.  Convert this to a cell array of strings.
            get_contents = @(f) f{1};
            c = cellfun(get_contents, tokens, 'UniformOutput', false);

            % Remove trailing whitespace characters from strings in c.  This can happen
            % if multiple whitespace characters separated the last word on a line from
            % the first word on the following line.
            c = deblank(c);
        end
         
        function [s,reduce]=ReduceSyllable(s, reduce)
            n4=length(s)-2;
            if n4<reduce
                s=[s(1) '*'];
                reduce=reduce-n4;
            else
                s=[s(1:floor(reduce)) '*'];
                reduce=0;
            end
        end
        
        function word=Toggle(b, yes, no)
            if nargin==1
                if b
                    word='Hide';
                else
                    word='Show';
                end
            else
                if b
                    word=yes;
                else
                    word=no;
                end
                
            end
        end
        
        function uri=ToFileUrl(s)
            if ispc
                uri=char(java.io.File(s).toURI);
            else
                uri=String.ToSystem(char(java.io.File(s).toURI));
            end
        end
        
        function str=Replace1(str, prev, next)
            idx=String.IndexOf(str, prev);
            if idx>0
                str=[str(1:idx-1) next str(idx+length(prev):end)];
            end
        end
        
        function where=DifWhere(s1, s2, start)
            if nargin<3
                start=1;
            end
            n1=length(s1);
            n2=length(s2);
            if n1>n2
                N=n2;
            else
                N=n1;
            end
            for i=start:N
                if s1(i)~=s2(i)
                    where=i;
                    return;
                end
            end
            where=[];
        end
        
        function str=toString(anyValue, nameIfStruct, delimiterIfStruct)
            if isempty(anyValue)
                str='';
            elseif ischar(anyValue)
                str=anyValue;
            elseif isnumeric(anyValue)
                str=num2str(anyValue);
            elseif islogical(anyValue)
                str=num2str(anyValue);
            elseif iscell(anyValue)
                if ischar(anyValue{1})
                    str=StringArray.toString(anyValue);
                elseif isnumeric(anyValue{1})
                    str=num2str(cell2mat(anyValue));
                end
            elseif isstruct(anyValue)
                c=struct2cell(anyValue);
                if nargin<3
                    delimiterIfStruct=',';
                    if nargin<2
                        nameIfStruct=false;
                    end
                end
                N=length(c);
                str='';
                if nameIfStruct
                    n=fieldnames(anyValue);
                    for i=1:N
                        str=[str n{i} delimiterIfStruct ...
                            String.toString(c{i}, nameIfStruct) ...
                            delimiterIfStruct ];
                    end
                    
                else
                    for i=1:N
                        str=[str String.toString(c{i}, nameIfStruct) ...
                            delimiterIfStruct ];
                    end
                end
            else
                error('Unconvertable value');
            end
        end
    end
    
end
