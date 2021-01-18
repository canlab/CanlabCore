%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%

classdef Html
    properties(Constant)
        SORT_REGEX='(?<=<Sort name=")([^"]+)" type="([^"]+)" value="([^"]+)">';
        SORT_REGEX_FMT='(?<=<Sort name=")%s" type="([^"]+)" value="([^"]+)">';
        YELLOW= '#FFFF00';
        RED=    '#FF0000'
        MAROON= '#8000000';
        BLACK=  '#000000';
        GRAY=   '#808080';
        OLIVE=  '#808000';
        GREEN=  '#00FF00';
        AVACODO='#008000';
        CYAN=   '#00FFFF';
        TEAL=   '#008080';
        BLUE=   '#0000FF';
        NAVY=   '#000080';
        MAGENTA='#FF00FF';
        PURPLE= '#800080';
        ORANGE= '#FFBF00';
        GOLD=   '#FF8000';
        COLORS={Html.NAVY, Html.BLACK, Html.BLUE, Html.ORANGE, Html.RED, ...
            Html.MAGENTA, Html.OLIVE, Html.GRAY, Html.MAROON, ...
            Html.YELLOW, Html.CYAN, Html.GREEN, Html.PURPLE, Html.AVACODO,...
            Html.TEAL, Html.GOLD};
        NCOLORS=length(Html.COLORS);
    end
    
    methods(Static)
        function [values, types]=DecodeSortValues(str, name)
            fmt=sprintf(Html.SORT_REGEX_FMT, name);
            if isa(str, 'javax.swing.JPanel')
                try
                    str=char(str.getComponent(0).getText);
                catch
                    str='';
                end
            end
            tokens=regexp(str, fmt, 'tokens');
            N=length(tokens);
            types=cell(1,N);
            values=cell(1,N);
            for i=1:N
                values{i}=tokens{i}{2};
                types{i}=tokens{i}{1};
            end
        end
        
        function [keys, types, values]=DecodeSort(str)
            if isa(str, 'javax.swing.JPanel')
                try
                    str=char(str.getComponent(0).getText);
                catch
                    str='';
                end
            end
            tokens=regexp( str, Html.SORT_REGEX, 'tokens');
            N=length(tokens);
            keys=cell(1,N);
            types=cell(1,N);
            values=cell(1,N);
            for i=1:N
                keys{i}=tokens{i}{1};
                types{i}=tokens{i}{2};
                values{i}=tokens{i}{3};
            end
        end
        
        function str=EncodeSort(key, value, type)
            if nargin<3
            	if isnumeric(value)
                    type='N';
                else
                    type='C';
                end
            end
            if isnumeric(value)
                value=num2str(value);
            end
            str=['<Sort name="' key '" type="' type '" value="'...
                value '">'];
        end
        
        function html=RotatedStyle
            html=[...
                '<STYLE>'...
                'th.rotate {'...
                '  /* Something you can count on */'...
                ' height: 140px;'...
                ' white-space: nowrap;'...
                '}'...
                'th.rotate > div {'...
                '  transform: '...
                '    /* Magic Numbers */'...
                '    translate(0px, 51px)'...
                '    rotate(315deg);'...
                '  width: 30px;'...
                '}'...
                'th.rotate > div > span {'...
                '  border-bottom: 1px solid #ccc;'...
                '  padding: 5px 10px;'...
                '}</STYLE>'];
            
            
        end
        
        function html=Rotate(txts)
            html='';
            if iscell(txts)
            N=length(txts);
            for i=1:N
                html=[html '<th class="rotate"><div><span>'...
                    txts{i} '</span></div></th>'];
            end
            else
                html=['<th class="rotate"><div><span>'...
                    txts '</span></div></th>'];
            end
        end
        
        function html=TestRotate
            rows={{'Killarney', 'Black', '9'}, ...
                {'Fergus', 'Gold', '5'},...
                {'Pepper', 'Gold', '13'}};
            cols={'Doggie name', 'Color', '# of years'};
            html=['<html>' Html.RotatedStyle '<table>'];
            html=[html '<thead>' Html.Rotate(cols) '</thead>'];
            html=[html Html.TableBody(rows, false) '</table></html>'];
            
        end
        
        function html=TableBody(rows, headerDone)
            if nargin<2
                headerDone=true;
            end
            if headerDone
                html='<table>';
            else
                html='';
            end
            nRows=length(rows);
            for i=1:nRows
                row=rows{i};
                nCols=length(row);
                html=[html '<tr>'];
                for j=1:nCols
                    html=[html '<td>' row{j} '</td>'];
                end
                html=[html '</tr>'];
            end
            if headerDone
                html=[html '</table>'];
            end
        end
        function RotateTable(cols, rows)
            style=Html.RotatedStyle;
            hdr=Html.Rotations(cols);
        end
        
        function html=Vertical(txt, pfx, sfx, limit)
            if nargin<4
                limit=9;
                if nargin<3
                    pfx='';
                    sfx='';
                end
            end
            txt=char(edu.stanford.facs.swing.Basics.RemoveXml(txt));
            N=length(txt);
            if N>limit
                txt=[txt(1:limit) '*'];
                N=limit+1;
            end
            html='<table cellspacing="0" cellpadding="0" border="0"><tr><td align="center">';
            for i=1:N
                html=[html pfx txt(i) sfx '<br>'];
            end
            html=[html '</td></tr></table>'];
        end
        
        function html=Verticals(txts, pfx, sfx)
            if nargin<2
                pfx='';
                sfx='';
            end
            N=length(txts);
            html='<table valign="bottom" cellpadding="0" cellspacing="0"><tr>';
            for i=1:N
                html=[html '<td>' Html.Vertical(txts{i}, pfx, sfx) '</td>'];
            end
            html=[html '</tr></table>'];
        end

        function str=H2(str, hrToo)
            app=BasicMap.Global;
            if nargin<2 || hrToo
                hr='<hr>';
            else
                hr='';
            end
            str=[app.h2Start str app.h2End hr];
        end
        
        function str=Small2(str)
            str=Html.Small(String.ToHtml(str));
        end
        
        function str=Small(str)
            app=BasicMap.Global;
            str=[app.smallStart str app.smallEnd];
        end
        
        function str=Wrap(str)
            str=['<html>' str '</html>'];
        end
        
        function str=WrapC(str)
            str=['<html><center>' str '</center></html>'];
        end
        
        function str=WrapSmall(str)
            app=BasicMap.Global;
            str=['<html>' app.smallStart str app.smallEnd '</html>'];
        end
        
        function str=WrapSmallTags(str)
            app=BasicMap.Global;
            str=[app.smallStart str app.smallEnd ];
        end
        
        function str=WrapHr(str)
            str=['<html><table border="0"><tr><td></td><td align="center">' ...
                str '</td><td></td></tr></table><hr></html>'];
        end

        function str=Hr(str)
            str=['<html><center>' str '</center><hr></html>'];
        end

        function str=WrapClr(words, clr, otherFontAttributes)
            if nargin<3
                str=['<html><font ' Html.Color(clr) '>' ...
                    words '</font></html>'];
            else
                str=['<html><font ' Html.Color(clr) otherFontAttributes ...
                    ' >' words '</font></html>'];
            end
        end

        function str=WrapClr2(words, clr, otherFontAttributes)
            if nargin<3
                str=['<font ' Html.Color(clr) '>'  words '</font>'];
            else
                str=['<font ' Html.Color(clr) otherFontAttributes ...
                    ' >' words '</font>'];
            end
        end

        function str=WrapSm(str, app)
            if nargin<2
                app=BasicMap.Global;
            end
            str=['<html>' app.smallStart str app.smallEnd '</html>'];
        end

        function str=WrapBoldSmall(str, app)
            if nargin<2
                app=BasicMap.Global;
            end
            str=['<b>' app.smallStart str app.smallEnd '</b>'];
        end

        function colors=ColorIds(ids)
            u=unique(ids);
            N=length(ids);
            colors=cell(1, N);
            for i=1:N
                idx=find(u==ids(i));
                colors{i}=Html.Color(idx);
            end
        end
        
        function s=Color(ml)
            if length(ml)==1
                idx=mod(ml, Html.NCOLORS)+1;
                s=Html.COLORS{idx};
                return;
            end
            ml=floor(ml*255);
            try
                s=sprintf('color="#%s%s%s"', hex(ml(1)), ...
                    hex(ml(2)), hex(ml(3)));
            catch ex
                s='color="black"';
            end
            function c=hex(c)
                c=dec2hex(c);
                if length(c)==1
                    c=['0' c];
                end
            end
        end

        function Browse(html)
            web( File.SaveTempHtml(html), '-browser');
        end
        
        function html=TitleMatrix(m, decimals, prefix, num2strThreshold)
            [R,C]=size(m);
            html=prefix;
            for r=1:R
                for c=1:C
                    html=[html '&#009;' ...
                        String.encodeRounded(m(r,c), decimals, true, num2strThreshold)];
                end
                html=[html '&#013;&#010;'];
            end
        end
        
        function img=ImgFile(f)
            img=['<img src=''file:/' ...
                        edu.stanford.facs.swing.Basics.EncodeFileUrl(f)...
                        '''>'];
        end
        
        function html=MatrixColored(rowHdrs, colHdrs, data, colors, ...
                emphasizeColumn, encode)
            if nargin<6
                encode=@encode_;
                if nargin<5
                    emphasizeColumn=-1;
                    if nargin<4
                        colors={};
                    end
                end
            end            
            html='<table cellpadding="5"><tr><td></td>';
            nCols=length(colHdrs);
            nRows=length(rowHdrs);
            for col=1:nCols
                html=[html '<td bgcolor="#AABDFF" align="center">' ...
                    encodeHead(colHdrs{col}) '</td>'];
            end
            html=[html '</tr>'];
            for row=1:nRows
                html=[html '<tr><td bgcolor="#FFFFDD">' ...
                    encodeHead(rowHdrs{row}) '</td>'];
                if isempty(colors)
                    clr='black';
                else
                    clr=colors{row};
                end
                for col=1:nCols
                    if col==emphasizeColumn
                        pre='<i><b>';
                        suf='</b></i>';
                    else
                        pre='';
                        suf='';
                    end
                    str=feval(encode, row, col, data(row, col));
                    html=[html '<td align="right"><font color="'...
                        clr '">' pre String.Pad(str, 9) ...
                        suf '</font></td>'];
                end
                html=[html '</tr>'];
            end
            html=[html '</table>'];
            function num=encode_(row, col, num)
                num=String.encodeRounded(num, 2, true);
            end
            function value=encodeHead(value)
                if isnumeric(value)
                    if size(value,1)>1
                        value=num2str(value');
                    else
                        value=num2str(value);
                    end
                end
            end
                    
        end
        

        function html=Matrix(rowHdrs, colHdrs, data, encode)
            if nargin<4
                encode=@encode_;
            end
            [R, C]=size(data);
            if isempty(rowHdrs) 
                rowHdrs=StringArray.Num2Str(1:R);
            end
            if isempty(colHdrs)
                colHdrs=StringArray.Num2Str(1:C);
            end
            html='<table><tr><td></td>';
            nCols=length(colHdrs);
            nRows=length(rowHdrs);
            for col=1:nCols
                html=[html '<td bgcolor="#AABDFF" align="center">' ...
                    String.ToHtml(colHdrs{col}) '</td>'];
            end
            html=[html '</tr>'];
            for row=1:nRows
                html=[html '<tr><td bgcolor="#FFFFDD">' rowHdrs{row} '</td>'];
                for col=1:nCols
                    str=feval(encode, row, data(row, col));
                    html=[html '<td align="right">'...
                        String.Pad(str, 9) '</td>'];
                end
                html=[html '</tr>'];
            end
            html=[html '</table>'];
            function num=encode_(col, num)
                num=String.encodeRounded(num,2,true);
            end

        end
        
        function img=Samusik(scale, forBrowser)
            file=Gui.SamusikIconFile;
            folder=Gui.SamusikIconFolder;
            img=Html.Img(file, folder, scale, forBrowser);
        end

        function img=Organizer(scale, forBrowser, filePath)
            if nargin<2
                forBrowser=false;
                if nargin<1
                    scale=.12;
                end
            end
            %file='plottype-dendrogram.png';
            %folder=fullfile(matlabroot,'/toolbox/matlab/icons/');
            file='tree2.png';
            folder=BasicMap.Global.contentFolder;
            if nargin<3 || ~filePath
                img=Html.Img(file, folder, scale, forBrowser);
            else
                img=fullfile(folder,file);
                img=edu.stanford.facs.swing.Basics.GetResizedImg(...
                    java.io.File(img), scale, ...
                    java.io.File(BasicMap.Global.appFolder));
            end
        end

        function img=Microscope(scale, forBrowser)
            file='microScope.png';
            folder=BasicMap.Global.contentFolder;
            img=Html.Img(file, folder, scale, forBrowser);
        end

        function img=Img(file, folder, scale, forBrowser)
            if nargin<4
                forBrowser=false;
                if nargin<3
                    if ispc && BasicMap.Global.highDef
                        scale=.2;
                    else
                        scale=.11;
                    end
                    if nargin<2 
                        folder=BasicMap.Global.contentFolder;
                    end
                end
            end
            if isempty(folder)
                folder=BasicMap.Global.contentFolder;
            end
            if scale==1
                f=fullfile(folder,file);
                if forBrowser
                    img=Html.ImgForBrowser(f,'');
                else
                    img=['<img src=''file:/' ...
                        edu.stanford.facs.swing.Basics.EncodeFileUrl(f)...
                        '''>'];
                end
            else
                img=edu.stanford.facs.swing.Html.ImgSized3(...
                    file, folder, scale, forBrowser);
            end
        end
        
        function img=ImgXy(file, folder, scale, forBrowser)
            if nargin<4
                forBrowser=false;
                if nargin<3
                    scale=1;
                    if nargin<2 
                        folder=BasicMap.Global.contentFolder;
                    end
                end
            end
            if isempty(folder)
                folder=BasicMap.Global.contentFolder;
            end
            if scale==1
                f=fullfile(folder,file);
                if forBrowser
                    img=Html.ImgForBrowser(f,'');
                else
                    img=['<img src=''file:/' ...
                        edu.stanford.facs.swing.Basics.EncodeFileUrl(f)...
                        '''>'];
                end
            else
                img=edu.stanford.facs.swing.Html.ImgSizedXy3(...
                    file, folder, scale, forBrowser);
            end
        end
        function s=HexColor(ml)
            s=Gui.HtmlHexColor(ml);
        end
        function img=ImgForBrowser(file, folder)
            if nargin<2
                folder=BasicMap.Global.contentFolder;
            end
            f=fullfile(folder, file);
            q='''';
            if String.Contains(f, '''')
                q='"';
            end
            f=regexprep(f, '#', '%23');
            img=['<img src=' q 'file:' f q '>'];
        end
        
        function html=To2Lists(strs1, strs2, type, hdr1, hdr2, ...
                highlightDifferences, limit)
            if nargin==7
                strs1=StringArray.Trim(strs1, limit);
                strs2=StringArray.Trim(strs2, limit);
            end
            if nargin>5 && highlightDifferences
                list1=Html.ToListDiff(strs1, strs2, type, false);
            else
                list1=Html.ToList(strs1, type, false);
            end
            if nargin>5 && highlightDifferences
                list2=Html.ToListDiff(strs2, strs1, type, false);
            else
                list2=Html.ToList(strs2, type, false);
            end
            html=['<tr><td valign="top">' list1 '</td><td valign="top">'...
                list2 '</td></tr>'];
            if nargin>3
                html=['<tr><td align="center" bgcolor="white"><b><u>'...
                    hdr1 '</u></b></td><td align="center" '...
                    'bgcolor="white"><b><u>' hdr2 '</u></b></td></tr>' html];
            end
            html=['<table bgcolor="white" >' html '</table>'];
        end
        
        function str=ToListDiff(strs, strs2, type, noListIf1, ...
                noListStart, noListEnd)
            if nargin<4
                if nargin<3
                    type='ul';
                    noListIf1=true;
                else
                    noListIf1=false;
                end
            end
            N=length(strs);
            if N==1 && noListIf1
                if nargin<6
                    noListStart='<b>';
                    noListEnd='</b>';
                end
                str=[noListStart String.ToHtml(strs{1}) noListEnd];
            else
                str=['<' type '>'];
                for i=1:N
                    if StringArray.Contains(strs2, strs{i})
                        str=[str '<li>' String.ToHtml(strs{i})];
                    else
                        str=[str '<li><font color="red"><b>' ...
                            String.ToHtml(strs{i}) '</b></font>'];
                    end
                end
                str=[str '</' type '>' ];
            end
        end
        
        function str=ToLimitedList(strs, limit, listType)
            if nargin<3
                listType='ol';
            end
            str=Html.ToList(strs, listType, true, '<b>', '</b>', limit);
        end
        
        function str=ToList(strs, listType, noListIf1, noListStart,...
                noListEnd, limit)
            if nargin<3
                if nargin<2
                    listType='ul';
                    noListIf1=true;
                else
                    noListIf1=false;
                end
            end
            N=length(strs);
            if N==1 && noListIf1
                if nargin<5 || isempty(noListStart)
                    noListStart='<b>';
                    noListEnd='</b>';
                end
                str=[noListStart String.ToHtml(strs{1}) noListEnd];
            else
                str=['<' listType '>'];
                if nargin>=6 && N>limit
                    for i=1:limit
                        str=[str '<li>' String.ToHtml(strs{i})];
                    end
                    str=[str '<li>' num2str(N-limit) ' more...'];
                else
                    for i=1:N
                        str=[str '<li>' String.ToHtml(strs{i})];
                    end
                end
                str=[str '</' listType '>' ];
                
            end
        end
        
        function str=ToListWithHtml(strs, type, noListIf1, noListStart, noListEnd)
            if nargin<3
                if nargin<2
                    type='ul';
                    noListIf1=true;
                else
                    noListIf1=false;
                end
            end
            N=length(strs);
            if N==1 && noListIf1
                if nargin<5
                    noListStart='<b>';
                    noListEnd='</b>';
                end
                str=[noListStart strs{1} noListEnd];
            else
                str=['<' type '>'];
                for i=1:N
                    str=[str '<li>' strs{i}];
                end
                str=[str '</' type '>' ];
            end
        end
        
        function list=WrapList(list, prefix, suffix)
            if nargin<2
                prefix='<small>';
                suffix='</small>';
            end
            N=length(list);
            for i=1:N
                list{i}=[prefix list{i} suffix];
            end
        end
        
        function str=ToList2(strs, type, noListIf1, noListStart, noListEnd)
            if nargin<3
                if nargin<2
                    type='ul';
                    noListIf1=true;
                else
                    noListIf1=false;
                end
            end
            N=length(strs);
            if N==1 && noListIf1
                if nargin<5
                    noListStart='<b>';
                    noListEnd='</b>';
                end
                str=[noListStart strs{1} noListEnd];
            else
                str=['<' type '>'];
                for i=1:N
                    str=[str '<li>' strs{i}];
                end
                str=[str '</' type '>' ];
            end
        end

        
        function uri=ToFileUrl(s)
            if ispc
                uri=char(java.io.File(s).toURI);
            else
                file=java.io.File(s);
                if ~file.isAbsolute
                    fullpath=fullfile(pwd, s);
                    file=java.io.File(fullpath);
                end
                uri=String.ToSystem(char(file.toURI));
            end
        end
        
        function BrowseFile(fileName, convert)
            if nargin>1 && convert && ismac
                str=File.ReadTextFile(fileName);
                str=strrep(str, 'file:/%2F', 'file:/');
                str=strrep(str, '%2F', filesep);
                fileName=File.SaveTempHtml(str);
            end
            web(Html.ToFileUrl(fileName), '-browser');
        end
        
        function ok=BrowseString(lines)
            try
                file=[tempname '.html'];
                fid=fopen(file, 'wt');
                fprintf(fid, '%s\n', char(lines));
                fclose(fid);
                Html.BrowseFile(file);
                ok=true;
            catch ex
                disp(ex.message);
                ok=false;
            end
        end
        
        function out=remove(in)
            in=strtrim(in);
            if length(in)>13 && ...
                strcmpi(in(1:6), '<html>') && strcmpi(in(end-6:end), '</html>')
                    out=in(7:end-7);
            else
                out=in;
            end
        end
       
        function str=Symbol(clr, sz, wrapHtml)
            if nargin<3
                wrapHtml=true;
            end
            str=['<font  ' Gui.HtmlHexColor(clr)...
                '>&bull;</font>'];
            sz=((sz-7)/4)+2;
            if sz>10
                sz=10;
            end
            sz=String.encodeInteger(ceil(sz) );
            if wrapHtml
                str=['<html><font size="' sz '">' str '</font></html>'];
            else
                str=['<font size="' sz '">' str '</font>'];
            end
        end

        function img=TempImg(fig, forBrowser, scale)
            if nargin<3
                scale=.8;
                if nargin<2
                    forBrowser=true;
                end
            end
            file=[tempname '.png'];
            try
                saveas(fig, file);
            catch 
            end
            [fldr, fn, ext]=fileparts(file);
            img=Html.ImgXy([fn ext], fldr, scale, forBrowser);
        end
        
        function out=H2Small2(h2, small)
            out=Html.WrapC([Html.H2(h2) Html.Small2(small)]);
        end
        
    end
    
end
