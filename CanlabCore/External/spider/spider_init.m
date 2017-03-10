

disp(' ');
disp('SPIDER : a machine learning toolbox for Matlab(R).');
disp(' ');
disp('This program is free software; you can redistribute it and/or');
disp('modify it under the terms of the GNU General Public License');
disp('as published by the Free Software Foundation; either version 2');
disp('of the License, or (at your option) any later version.');
disp(' ');
disp('This program is distributed in the hope that it will be useful,');
disp('but WITHOUT ANY WARRANTY; without even the implied warranty of');
disp('MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the');
disp('GNU General Public License for more details:');
disp('http://www.gnu.org/copyleft/gpl.html');
disp(' ');
disp(' ');
  
if isunix  
  
a=[ '~/matlab'];  
if exist(a)==0  
  disp('[matlab directory does not exist, creating it]');  
  a=unix(['mkdir ~/matlab']);  
end  
a=unix([' echo "path(' char(39)  pwd  char(39) ',path); spider_path;" >> ~/matlab/startup.m']);  
if a==0 disp('Installed!'); end;  
  
else   
  
str=[matlabroot '/toolbox/local/startup.m'];  
  
fid=fopen(str,'a');  
str=['path(' char(39)  pwd  char(39) ',path); spider_path;'];  
fprintf(fid,'%s\n',str);  
fclose(fid);  
  
disp('Installed!');  
  
end   
path(pwd,path); spider_path;   
