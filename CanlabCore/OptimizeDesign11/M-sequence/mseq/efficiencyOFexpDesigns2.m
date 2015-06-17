% efficiencyOFexpDesigns2 : OVERLAPPING EVENTS
% ER-fMRI data analysis
% script:
% throughout we assume TR=2.


loadDataYes=1; % =1 -> loads existing data file

pluginCorrYes=0;
daleYes=1; % 1-usual Dale efficiency. 0-Fisher

% event matrices are defined based on vectors with the average
% energy removed

cutEndPad=0; %cutting reduces efficiency!
nReps=10000; %n-2;
n2=24; % assumed # TRs to cover HDR fn
nVals=1;
pwRange=6:10; %7

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

mEff=zeros(1,n-1);

shift=2^(pwr-1);
%shift1=2^(pwr-2)*3;
%XXXXXXXXXXXXXXXXXXXXX
for cycle=1:(2^pwRange(1))-2,
  mEvent=ms; %row 1
  mEvent=[mEvent; [ms(shift+1:end) ms(1:shift)]];%row 3
  mEvent=[mEvent; [ms(cycle+1:end) ms(1:cycle)]];%row 2

  % defining event matrix as convolution matrix
  %eventMatrix=makeEventMtrx(mEvent-(ones(length(mEvent),1)*...
  %sum(mEvent'))'/length(mEvent),n2); 
	
  eventMatrix=makeEventMtrx(mEvent,n2); 
	
  if cutEndPad, eventMatrix=eventMatrix(1:n,:); end;
  eventMatrix=eventMatrix-ones(size(eventMatrix,1),1)* ...
    sum(eventMatrix)/size(eventMatrix,1);

  if daleYes,
    designEff=1/trace(inv(eventMatrix'*Cninv*eventMatrix));
  else
    designEff=trace(eventMatrix'*Cninv*eventMatrix);
  end;
	
  mEff(cycle)=designEff;
end;

%mProbab=sum(mEvent'); mProbab=mean(mProbab)/n;
plot(mEff,'k','LineWidth',2); 
set(gca,'FontSize',12);
xlabel('Cyclical shift of event vector #3','FontSize',16)
ylabel('Efficiency','FontSize',16)
%axis([0 n 0 1])
set(gca,'Position',[.2,.15,.7,.7]);
set(gcf,'PaperPosition',[1,1,5,4]);
%print -dpsc2 cyclingOverlapping3ev
drawnow;
pause

foo=max(mEff);
mEffAll=[mEffAll, foo(1)];

% for simply randomized designs:
nVals=size(mEvent,1);
nOnes=0.5;
overlapYes=1;

rndEff=[];
for k=1:nReps;
   %nOnes=k/(n-1);
      
	rndEvent=balancedRnd(n,nVals,nOnes,overlapYes);
	%if k==1, sum(rndEvent'), end;
        eventMatrix=makeEventMtrx(rndEvent-(ones(length(rndEvent),1)*sum(rndEvent'))'/length(rndEvent),n2);
   if cutEndPad, eventMatrix=eventMatrix(1:n,:); end;
        if daleYes,
	  designEff=1/trace(inv(eventMatrix'*Cninv*eventMatrix));
	else
	  designEff=trace(eventMatrix'*Cninv*eventMatrix);
	end;
	
	rndEff=[rndEff, designEff];
end;
rndEffAll=[rndEffAll;rndEff];
end;
toc;

else 
  eval(['load dataEff',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh']);
end %fi loadDataYes

nn=2.^[pwRange]-1;


figure(1); clf;
%[b,a]=hist(rndEff,20); 
%bar(a,b/nReps); hold on;
%plot(max(mEff),0,'k*'); 

loglog(nn,mEffAll,'r*');hold on;

% theoretical max:
%theoMax=(nn+1)/n2/4;
%loglog(nn,theoMax,'g+');hold on;
p99=max(rndEffAll'); %prctile(rndEffAll',99.9);
p00=min(rndEffAll'); %prctile(rndEffAll',0.1);
med=median(rndEffAll');
loglog(nn,med,'k.-','LineWidth',2);
loglog([nn;nn],[p00;p99],'k')
loglog([nn*.93;nn*1.1],[p00;p00],'k')
loglog([nn*.93;nn*1.1],[p99;p99],'k')

axis([50 10^3*1.2 10^-.5 10^1*3]);
%set(gca,'YTick',[0.003 .01 .1 .3 1 3 10 30]);
set(gca,'XTick',[2^6 2^7 2^8 2^9 2^10]);
set(gca,'LineWidth',2,'FontSize',12);

xlabel('Sequence length','FontSize',12);
ylabel('Efficiency','FontSize',12);
title(['With Overlap: n_e=',num2str(nVals),' p=',num2str(nOnes),' n_h=',num2str(n2)]);

set(gca,'Position',[.2,.15,.7,.7]);
set(gcf,'PaperPosition',[1,1,4,3]);
eval(['print -dpsc2 eff',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh_over']);
disp(['print -dpsc2 eff',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh_over']);



figure(2); clf;
semilogx(nn,mEffAll./max(rndEffAll'),'-ok','LineWidth',2); hold on;
semilogx(nn,mEffAll./med,':+k','LineWidth',2)
set(gca,'Position',[.2,.15,.7,.7]);
xlabel('Sequence length','FontSize',12);
ylabel('m-seq/random efficiency ratio','FontSize',12);
title(['With overlap: n_e=',num2str(nVals),' p=',num2str(nOnes),' n_h=',num2str(n2)]);
set(gca,'Position',[.2,.15,.7,.7]);
set(gcf,'PaperPosition',[1,1,4,3]);
set(gca,'XTick',[2^6 2^7 2^8 2^9 2^10]);
set(gca,'LineWidth',2,'FontSize',12);
axis([50 10^3*1.2 .9 2.5]);
eval(['print -dpsc2 eff',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nhRATIO']);
disp(['print -dpsc2 eff',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nhRATIO']);

eval(['save dataEff',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh rndEffAll p99 p00 med mEffAll']);
disp(['save dataEff',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh rndEffAll p99 p00 med mEffAll']);
