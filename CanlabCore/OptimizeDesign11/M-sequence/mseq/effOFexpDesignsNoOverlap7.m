% effOFexpDesignsNoOverlap7
% ER-fMRI data analysis
% script:
% throughout we assume TR=2.
clf; colormap(gray);

loadDataYes=1; % =1 -> loads existing data file

%cd /home/giedrius/giedrius/projects/matlab/erfmri
pluginCorrYes=0;
overlapYes=0;
nReps=10000; %
nVals=7;
n2=24; % assumed # TRs to cover HDR fn
pwRange=6:10; %10;

daleYes=1; % 1-usual Dale efficiency. 0-Fisher
% event matrices are defined based on vectors with the average
% energy removed
plotMtrx=0;

cutEndPad=1; %cutting reduces efficiency!

if ~loadDataYes,
rndEffAll=[];
mEffAll=[];

tic
for pwr=pwRange,
% creating event vectors
ms=m2bin([mseq(2,pwr,0,1)])'; % fixing to the length of 2^n
n=length(ms)  % scan duration in TRs
if pluginCorrYes
  b=[0.406;0.8825]; % parameters from SPG fMRI data
  fittedACorr=autocorrFnct(b,1:n+(n2-1)*rem(cutEndPad+1,2));
  Cninv=inv(toeplitz(fittedACorr));
else
  Cninv=eye(n+(n2-1)*rem(cutEndPad+1,2));
end;

sh =(n2:2*n2)+1;
sh1=(n2:2*n2)+2+n2;
mEff=zeros(length(sh),length(sh1));

i=0;
for shift =sh,
i=i+1;
j=0;
for shift1=sh1,
j=j+1;
        %shift=floor(n2+(n-2*n2)*rand(1)); 
	%shift1=shift;
        %while abs(shift-shift1)<n2,
	%  shift1=floor(n2+(n-2*n2)*rand(1));
	%end;
  
        %shift=i+n2; 
	%shift1=i*2+n2;
	
	foo=[ms(shift+1:end) ms(1:shift)];%shifted version
	foo1=[ms(shift1+1:end) ms(1:shift1)];%shifted version
	foo=ms*4+foo*2+foo1;

	mEvent=zeros(nVals,n);
	for k=1:nVals,
   	  ind=find(foo==k);
   	  mEvent(k,ind)=1;
	end;

	% defining event matrix as convolution matrix
	%eventMatrix=makeEventMtrx(mEvent-(ones(length(mEvent),1)*...
	%sum(mEvent'))'/length(mEvent),n2); 
	
	eventMatrix=makeEventMtrx(mEvent,n2); 
	
        if cutEndPad, eventMatrix=eventMatrix(1:n,:); end;
	eventMatrix=eventMatrix-ones(size(eventMatrix,1),1)* ...
	    sum(eventMatrix)/size(eventMatrix,1);
	if plotMtrx, imagesc(mEvent); drawnow; pause(.5); end;
	
   if daleYes,
	  designEff=1/trace(inv(eventMatrix'*Cninv*eventMatrix));
	else
	  designEff=trace(eventMatrix'*Cninv*eventMatrix);
	end;
	
	mEff(i,j)=designEff;
end;
end;

%mProbab=sum(mEvent'); mProbab=mean(mProbab)/n;
imagesc(mEff); colorbar;
set(gca,'FontSize',18);
xlabel('shift1');
ylabel('shift');
drawnow;

%pause
mEffAll=[mEffAll, max(max(mEff))];

% for simply randomized designs:
nVals=size(mEvent,1);
%n=n+1
if pluginCorrYes
  b=[0.406;0.8825]; % parameters from SPG fMRI data
  fittedACorr=autocorrFnct(b,1:n+(n2-1)*rem(cutEndPad+1,2));
  Cninv=inv(toeplitz(fittedACorr));
else
  Cninv=eye(n+(n2-1)*rem(cutEndPad+1,2));
end;

if pwr<8, correction=0.02; 
elseif pwr<9, correction=0.005;
elseif pwr==10, correction=0.001;
end;
correction=0;

pEv=nVals/(nVals+1)+correction; %

rndEff=zeros(1,nReps);
for k=1:nReps;
   rndEvent=balancedRnd(n,nVals,pEv,overlapYes);%sum(rndEvent')

   % defining event matrix as convolution matrix
   %eventMatrix=makeEventMtrx(rndEvent-(ones(length(rndEvent),1)*...
   %sum(rndEvent'))'/length(rndEvent),n2); 
	
   eventMatrix=makeEventMtrx(rndEvent,n2); 
   
   if cutEndPad, eventMatrix=eventMatrix(1:n,:); end;
   eventMatrix=eventMatrix-ones(size(eventMatrix,1),1)* ...
       sum(eventMatrix)/size(eventMatrix,1);
   
   if daleYes,
	  designEff=1/trace(inv(eventMatrix'*Cninv*eventMatrix));
   else
	  designEff=trace(eventMatrix'*Cninv*eventMatrix);
   end;
   rndEff(k)=designEff;
end;
   

rndEffAll=[rndEffAll;rndEff];
end;
toc;

else 
  eval(['load dataEffNOV',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh']);
end %fi loadDataYes


nn=2.^[pwRange]-1;

%[b,a]=hist(rndEff,20); 
%bar(a,b/nReps); hold on;
%plot(max(mEff),0,'k*'); 

figure(1); clf;
loglog(nn,mEffAll,'r*');hold on;

% theoretical max -- only for one event type:
%theoMax=(nn+1)/n2/4;
%loglog(nn,theoMax,'g+');hold on;

p99=prctile(rndEffAll',99.9);
p00=prctile(rndEffAll',0.1);
med=median(rndEffAll');
loglog(nn,med,'k.-','LineWidth',2);
loglog([nn;nn],[p00;p99],'k')
loglog([nn*.93;nn*1.1],[p00;p00],'k')
loglog([nn*.93;nn*1.1],[p99;p99],'k')

axis([50 10^3*1.2 10^-2 10^1]);
%set(gca,'YTick',[0.003 .01 .1 .3 1 3 10 30]);
set(gca,'XTick',[2^6 2^7 2^8 2^9 2^10]);
set(gca,'LineWidth',2,'FontSize',12);

xlabel('Sequence length','FontSize',12);
ylabel('Efficiency','FontSize',12);
title(['No Overlap: n_e=',num2str(nVals),' p=',num2str(pEv),' n_h=',num2str(n2)]);

set(gca,'Position',[.2,.15,.7,.7]);
set(gcf,'PaperPosition',[1,1,4,3]);
eval(['print -dpsc2 effNOV',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh']);
disp(['print -dpsc2 effNOV',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh']);

figure(2); clf;
semilogx(nn,mEffAll./max(rndEffAll'),'-ok','LineWidth',2); hold on;
semilogx(nn,mEffAll./med,':+k','LineWidth',2)
set(gca,'Position',[.2,.15,.7,.7]);
xlabel('Sequence length','FontSize',12);
ylabel('m-seq/random efficiency ratio','FontSize',12);
title(['No Overlap: n_e=',num2str(nVals),' p=',num2str(pEv),' n_h=',num2str(n2)]);
set(gca,'XTick',[2^6 2^7 2^8 2^9 2^10]);
set(gca,'LineWidth',2,'FontSize',12);
axis([50 10^3*1.2 .9 1.5]);
set(gca,'Position',[.2,.15,.7,.7]);
set(gcf,'PaperPosition',[1,1,4,3]);
eval(['print -dpsc2 effNOV',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nhRATIO']);
disp(['print -dpsc2 effNOV',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nhRATIO']);

eval(['save dataEffNOV',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh rndEffAll p99 p00 med mEffAll']);
disp(['save dataEffNOV',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh rndEffAll p99 p00 med mEffAll']);
