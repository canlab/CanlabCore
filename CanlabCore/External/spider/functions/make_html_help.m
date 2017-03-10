function []=make_html_help(fin)

while 1
if(fin(1)==10 | fin(1)==32) fin=fin(2:end); else break; end;
end
    
fin=strtrim(fin)

 ss=spider_path;
 s=help(fin);
 
% fid = fopen(['c:\windows/desktop/personal/web/spider/help_' fin '.html'],'w')
 fid = fopen(['help_' fin '.html'],'w')
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

f=[1 find(double(s)==10)]; % find line feeds

for i=1:length(f)-1
  if (f(i+1)-f(i))>5
      break;
  end   
end 
take=s(f(i):f(i+1));
if sum(take=='=')>30
  i=i+1; take=s(f(i):f(i+1));
end

if findstr(fin,take) 
    take=[char(fin + ('A'-'a')) ' object']; 
    take(take=='?')='_';
end;
ss=take;
a=['<TABLE BORDER=0 CELLSPACING=0 WIDTH="550" align=center> <TR> <TD WIDTH="100" BGCOLOR="red"><I><B> </B></I></TD> <TD WIDTH="490" BGCOLOR="red"><B> <FONT COLOR="white" size="+1">  ' ss '  </FONT></B> </TD></TR> </TABLE><br><br><pre>'];
fprintf(fid,'%s',a);

jump=i+1;
for i=jump:length(f)-1
  take=s(f(i):f(i+1)-1);
  if sum(take=='=')<30  %% ignore lines of '=' signs
    process(take,fid);
  end
end    
a='</pre></table>\n'; fprintf(fid,a);
fclose(fid);
        

function a=process(s,fid)
%------------------------
a=[s(1:length(s)) ]; doit=1;


%%------ search for reference,author,etc. ----------
if length(s)>12
  if ~isempty(findstr('Reference',s))
    a=char(10);
    fprintf(fid,'<br>%s',a);  %% make extra space
  end
  if findstr('Author',s)
end
 if ~isempty(findstr('Reference',s)) | ~isempty(findstr('Author',s)) | ~isempty(findstr('Link',s))
   ss=s; f=min(find(s==':'));
   if f+5<length(s)
    if  ~isempty(findstr('Link',s)) %% try to index actual website
       f=min(find(s==':')); 
       ss=[s(1:f) ' <a href= ' s(f+1:end) '  target="_blank"> ' s(f+1:end) ' </a>'];    
    end
    a=['</pre></pre><TABLE BORDER=0 CELLSPACING=0 WIDTH="550" align=center> <TR> <TD WIDTH="100" BGCOLOR="black"><I><B> </B></I></TD> <TD WIDTH="490" BGCOLOR="black"><B> <FONT COLOR="white" size="-1">  ' ss '  </FONT></B> </TD></TR>'];
    fprintf(fid,'%s\n',a);
     if  ~isempty(findstr('Link',s)) %% try to index actual website
           fprintf(fid,'</table>\n');  
     end
   end
   doit=0;
 end
end
%%----------------------------------


if doit  %% -----add font colors
    a(double(a)==10)=' ';aa=a; a=[a char(10)];
    ff=findstr('--',a);        
    f1=min(find(not(aa==' ')));

    if 0 %not(isempty(ff)) | f1>20
        f=min(find(a=='-')); f2=f+1;
        if isempty(ff) f=min(find(not(a==' '))); f2=f; end;
        a=s(2:f-1); a2=s(f2:length(s));
        if not(isempty(ff)) a=[a '-']; end;
        a2(a2=='-')=' ';
        s=['<font color="red" size="+0"><b>' a '</font></b><font color="green" face="Verdana" size="0" >' a2 '</font><br><br>'];
        fprintf(fid,'%s',s);
    else    
        %if findstr(a,'p1') 
         %a    
         %   keyboard; 
         %end;
        
        f1=min(find(not(aa==' ')));
        f2=max(find(not(aa==' ')));
        col='"red" ';
        if f1<20 col='"black"'; end; 
        if not(isempty(find(aa==':'))) | not(isempty(findstr(aa,'Methods'))) | not(isempty(findstr(aa,'Model'))) | not(isempty(findstr(aa,'Hyperparameters')))
            f=min(find(not(aa==' ')));
            if aa(f)>='A' & aa(f)<='Z'   col='"blue"'; end; 
        end; 
        if findstr(aa,'--') col='"red"'; end;
        if strcmp(col,'"blue"') a=['<b>' a '</b>'];end;
        s=['<font color=' col '  size="+0"  >' a '</font>'];
        fprintf(fid,'%s',s);
            
    end 
end    


a=1; return;


