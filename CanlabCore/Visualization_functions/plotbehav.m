function [blockbehavdata slope]=plotbehav(tmp,numcorr,replace)

if nargin<3
    replace=0;
end

if nargin<2;
    numcorr=1;
end

figure('color','w');

for b=1:numcorr;
subplot(numcorr,1,b);
numsub=length(tmp);
behavdata=zeros(numsub,size(tmp{1},1));
for n=1:numsub;
    behav=tmp{n};
    behavdata(n,1:length(behav))=behav(:,b);
end

blockbehavdata=behavdata;

%mean replace if necessary
if replace==1;
    % find subjects with zeros & replace
    for b=1:size(behavdata,2);
          wheregood=find(behavdata(:,b)~=0);
          wherebad=find(behavdata(:,b)==0);
          mbv=squeeze(mean(behavdata(wheregood,b)));
          if ~isempty(wherebad)
           disp(['replacing...']);
          blockbehavdata(wherebad,b)=squeeze(dyadic(mbv,size(wherebad,1)));
          end
    end
end

%blockbehavdata=ztransf(blockbehavdata);


steplot(blockbehavdata,'p');
mbb=squeeze(mean(blockbehavdata));
p=polyfit(1:length(mbb),mbb,1);
p1=polyval(p,1:length(mbb));
hold on;
plot(p1,'linewidth',3,'color','k');
XLim([0 length(mbb)+1])
end

%corblockbehavdata=blockbehavdata-(squeeze(dyadic(p1,16)));
%figure;steplot(corblockbehavdata,'p');


slope=p1;       %information to be partialled out of correlations: the slope
%slope=squeeze(mean(blockbehavdata));       % or the global mean?

fid=fopen('blockcorrect.txt','w');
for s=1:numsub;
    fprintf(fid, '%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n',squeeze(blockbehavdata(s,:)));
end
fclose (fid);
