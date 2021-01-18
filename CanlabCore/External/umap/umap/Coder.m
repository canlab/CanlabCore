%
%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
classdef Coder < handle
    methods(Static)
        function out=NumMatLab(name, matrix, lims, eol, isInt)
            if nargin<5
                isInt=false;
                if nargin<4
                    eol=[];
                    if nargin<3
                        lims=[nan nan];
                    end
                end
            end
            cnt=1;
            [R, C]=size(matrix);
            if length(lims)>=1 && ~isnan(lims(1))
                if lims(1)<R
                    R=lims(1);
                end
            end
            if length(lims)>=2 && ~isnan(lims(2))
                if lims(2)<C
                    C=lims(2);
                end
            end
            if isempty(eol)
                eol=R*C+1;
            end
            if C==1
                flipped='''';
                C=R;
                R=1;
                matrix=matrix';
            else
                flipped='';
            end
            all=R*C;
            sb=java.lang.StringBuilder(all);
            if all>1
                sb.append([name '=[']);
            else
                sb.append([name '=']);
            end
            for r=1:R
                for c=1:C
                    if isInt
                        num=floor(matrix(r,c));
                    else
                        num=matrix(r,c);
                    end
                    s=num2str(num);
                    if isempty(find(s=='.',1))
                        sb.append(int64(num));
                    else
                        sb.append(num);
                    end
                    if mod(cnt, eol)==0
                        sb.append(sprintf("\n   "));
                    else
                        sb.append(' ');
                    end     
                    cnt=cnt+1;
                end
                if R>1
                    sb.append('; ');
                end
            end
            if all>1
                sb.append([']' flipped ';']);
            else
                sb.append([flipped ';']);
            end
            out=char(sb.toString);
        end
        function out=NumJava(name, matrix, lims, eol, isInt)
            if nargin<5
                isInt=false;
                if nargin<4
                    eol=[];
                    if nargin<3
                        lims=[nan nan];
                    end
                end
            end
            cnt=1;
            [R, C]=size(matrix);
            if length(lims)>=1 && ~isnan(lims(1))
                if lims(1)<R
                    R=lims(1);
                end
            end
            if length(lims)>=2 && ~isnan(lims(2))
                if lims(2)<C
                    C=lims(2);
                end
            end
            if isempty(eol)
                eol=R*C+1;
            end
            if C==1
                C=R;
                R=1;
                matrix=matrix';
            end
            all=R*C;
            sb=java.lang.StringBuilder(all);
            if isInt
                sb.append('int ');
            else
                sb.append('double ');
            end
            
            if all>1
                if R>1
                    sb.append('[][] ');
                else
                    sb.append('[] ');
                end
                sb.append([name '={' char(sprintf("\n\t"))]);
            else
                sb.append([name '=']);
            end
            justDidEol=false;
            for r=1:R
                if r>1 && ~justDidEol
                    sb.append(', ');
                end
                if R>1
                    sb.append('{');
                end
                for c=1:C
                    if isInt
                        num=floor(matrix(r,c));
                    else
                        num=matrix(r,c);
                    end
                    s=num2str(num);
                    if isempty(find(s=='.',1))
                        sb.append(int64(num));
                    else
                        sb.append(num);
                    end
                    if mod(cnt, eol)==0
                        justDidEol=true;
                        if c==C && R>1
                            sb.append('}');
                        end
                        sb.append(',');
                        sb.append(sprintf("\n\t"));
                    else
                        justDidEol=false;
                        if c<C
                            sb.append(', ');
                        end
                    end     
                    cnt=cnt+1;
                end
                if ~justDidEol && R>1
                    sb.append('}');
                end
            end
            if all>1
                sb.append('};');
            else
                sb.append(';');
            end
            out=char(sb.toString);
        end
    end
end
