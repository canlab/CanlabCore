function []=make_html_contents
disp('note that we assume to be in a directoy called spider');
'going..'

 ss=spider_path;

f=fopen('core_algorithms.txt','rt');
fseek(f,0,'eof');
l=ftell(f);
fseek(f,0,'bof');
core_algos=char(fread(f,l)');
fclose(f);

 
 fin = fopen([ss 'Contents.m'],'r')
 F = fread(fin);
 s = char(F');
 fclose(fin);
 
% fid = fopen(['c:\spider/doc/web/objects.html'],'w')
 fid = fopen(['objects.html'],'w')
 fprintf(fid,'<HTML>\n');
 fprintf(fid,'<HEAD>\n');
 fprintf(fid,'<HEAD>\n');
 a='<BODY BGCOLOR="white">\n'; fprintf(fid,a);
a='<BODY TEXT="blue">\n'; fprintf(fid,a);
a='<BODY LINK="red">\n'; fprintf(fid,a);
a='<BODY VLINK="red">\n'; fprintf(fid,a);
a='<BODY ALINK="blue">\n'; fprintf(fid,a);
a='<BLOCKQUOTE> \n'; fprintf(fid,a);
a=' \n'; fprintf(fid,a);
a='<title> SPIDER </title>\n'; fprintf(fid,a);
a=' \n'; fprintf(fid,a);
a='<font color="Black">\n'; fprintf(fid,a);
a='<font size=+5 face="Impact"> The Spider Objects\n'; fprintf(fid,a);
a='</font></font> \n'; fprintf(fid,a);
a='<font size=-2 face="Verdana">\n'; fprintf(fid,a);
a='</br>\n'; fprintf(fid,a);
a=' \n'; fprintf(fid,a);
a='<div align=left>\n'; fprintf(fid,a);
a='</div>\n'; fprintf(fid,a);
a='<p>\n'; fprintf(fid,a);
a='</font>\n'; fprintf(fid,a);
a=' \n'; fprintf(fid,a);
a='<font color="blue" size="+0">\n'; fprintf(fid,a);
a='  \n'; fprintf(fid,a);
a='<pre>\n'; fprintf(fid,a);

f=[find(s=='%')]; s(f)=' '; % purge comments
s=s(find(s~=10)); 
f=[find(double(s)==13)]; % find line feeds
for i=2:length(f)-1
   process(s(f(i):f(i+1)-1),fid,core_algos);
end    
a='</pre>\n'; fprintf(fid,a);
fprintf(fid,'<font color="Black"> Objects marked with <font color="Red">-</font> are available in the extra package. </font>')
fclose(fid);
        

function a=process(s,fid,core_algos)
%------------------------
a=[s(1:length(s)) ]; doit=1;


%a; keyboard

%%------ search for title ----------
if findstr('objects',a) & not(s(4)==' ')
    a=['<TABLE BORDER=0 CELLSPACING=0 WIDTH="550" align=center> <TR> <TD WIDTH="100" BGCOLOR="#000022"><I><B> </B></I></TD> <TD WIDTH="490" BGCOLOR="#000022"><B> <FONT COLOR="white" size="+1">  ' char(s) '  </FONT></B> </TD></TR> </TABLE>'];
    fprintf(fid,'\n');
    fprintf(fid,'%s',a);
    doit=0;
end
%%----------------------------------
    
if doit  %% -----add font colors
    f=min(find(a=='-'));
    if not(isempty(f))
        a=s(2:f-1); a2=s(f+1:length(s));
        a=strtrim(a);
        a2=strtrim(a2);
        f1=min(find(not(a==' ')));
        f2=max(find(not(a==' ')));
        %f1=min(find((a>='a' & a<='z')));
        %f2=max(find((a>='a' & a<='z')));
        a_bef=a(1:f1-1); a_aft=a(f2+1:length(a));
        a=a(f1:f2);
        obj_a=a;
        
%         make_html_help(a);
        a=[a_bef '<a href="help_' a '.html">' a ' </a>' a_aft ];
        s=['<TABLE BORDER=0 CELLSPACING=0 WIDTH="550" align=center colspan=2> <TR> <td style="width: 30%;"> <font color="red" size="+0"><b>' a '</font></td> <td>' ispartofcore(obj_a,core_algos) '</td> <td style="width: 65%;"> <font color="green" face="Verdana" size="0" >' a2 char(13) '</font> </td> </tr></TABLE>'];
        fprintf(fid,'%s',s);
    end 
%    fprintf(fid,'<br>\n');
end    

a=1; return;



function s=ispartofcore(a,core_algos)

if(isempty(strfind(core_algos,a)))
    s='-';
else
    s='+';
end