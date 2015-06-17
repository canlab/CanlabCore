% effOFexpDesignsNoOverlap3
% ER-fMRI data analysis
% script:
% throughout we assume TR=2.

loadDataYes=1; % =1 -> loads existing data file

%cd /home/giedrius/giedrius/projects/matlab/erfmri
pluginCorrYes=0;
cutEndPad=1; %cutting reduces efficiency!
overlapYes=0;
nReps=10000; %
shiftYes=0;

daleYes=1; % 1-usual Dale efficiency. 0-Fisher
% event matrices are defined based on vectors with the average
% energy removed

n2=12; % assumed # TRs to cover HDR fn
%pwRange=4:6;
pwRange=3:4;
base=5;

nVals=base-1;
if base==3, fixEv=1; elseif base==5, fixEv=2; end;

if ~loadDataYes,
rndEffAll=[];
mEffAll=[];

tic
for pwr=pwRange,
% creating event vectors
ms=mseq(base,pwr,0,1)'+fixEv; % fixing to the length of 2^n
n=length(ms)  % scan duration in TRs
if pluginCorrYes
  b=[0.406;0.8825]; % parameters from SPG fMRI data
  if ~cutEndPad, nl=n+n2-1; else nl=n; end;
  fittedACorr=autocorrFnct(b,1:nl);
  Cninv=inv(toeplitz(fittedACorr));
else
 if ~cutEndPad, Cninv=eye(n+n2-1);
 else Cninv=eye(n);
 end;
end;

mEff=zeros(1,n-1);

if shiftYes,
shift1=2^(pwr-2)+11;
for shift=1:(2^pwr(1))/2,
	%shift=2^(pwr-1)+5;
	foo=[ms(shift+1:end) ms(1:shift)];%shifted version
	%foo1=[ms(shift1+1:end) ms(1:shift1)];%shifted version
	foo=ms+foo*2; %+foo1*4;

	mEvent=zeros(nVals,n);
	for k=1:nVals,
   	  ind=find(foo==k);
   	  mEvent(k,ind)=1;
	end;

	eventMatrix=makeEventMtrx(mEvent-(ones(length(mEvent),1)*sum(mEvent'))'/length(mEvent),n2); % defining event matrix as convolution matrix
        if cutEndPad, eventMatrix=eventMatrix(1:n,:); end;
	if daleYes,
	  designEff=1/trace(inv(eventMatrix'*Cninv*eventMatrix));
	else
	  designEff=trace(eventMatrix'*Cninv*eventMatrix);
	end;
	
	mEff=[mEff, designEff];
end;


%mProbab=sum(mEvent'); mProbab=mean(mProbab)/n;
plot(mEff,'k','LineWidth',2); 
set(gca,'FontSize',18);
xlabel('Cyclical shift of event vector #3','FontSize',24)
ylabel('Efficiency','FontSize',24)
drawnow;
%pause

else % if no shifts are needed
  foo=ms;
  mEvent=zeros(nVals,n);
  for k=1:nVals,
    ind=find(foo==k);
    mEvent(k,ind)=1;
  end;
  eventMatrix=makeEventMtrx(mEvent-(ones(length(mEvent),1)*sum(mEvent'))'/length(mEvent),n2); % defining event matrix as convolution matrix
  if cutEndPad, 
    eventMatrix=eventMatrix(1:n,:); 
  end;
  if daleYes,
    designEff=1/trace(inv(eventMatrix'*Cninv*eventMatrix));
  else
    designEff=trace(eventMatrix'*Cninv*eventMatrix);
  end;
  mEff=designEff;
    
end;
mEffAll=[mEffAll, max(mEff)];


% for simply randomized designs:
nVals=size(mEvent,1);
%n=n+1

if pwr<8, correction=0.02; 
elseif pwr<9, correction=0.005;
elseif pwr==10, correction=0.001;
end;
pEv=nVals/(nVals+1)+0; %correction; %

rndEff=zeros(1,nReps);
for k=1:nReps;
   rndEvent=balancedRnd(n,nVals,pEv,overlapYes);%sum(rndEvent')
   eventMatrix=makeEventMtrx(rndEvent-(ones(length(rndEvent),1)*sum(rndEvent'))'/length(rndEvent),n2);
   if cutEndPad, eventMatrix=eventMatrix(1:n,:); end;
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
  n=base.^pwRange(end)-1;
  eval(['load dataEffNOV',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh',num2str(n),'n']);
end %fi loadDataYes

nn=base.^[pwRange]-1;

%[b,a]=hist(rndEff,20); 
%bar(a,b/nReps); hold on;
%plot(max(mEff),0,'k*'); 

figure(1); clf;
loglog(nn,mEffAll,'r*');hold on;

% theoretical max:
%theoMax=(nn+1)/n2/4;
%loglog(nn,theoMax,'g+');hold on;

p99=prctile(rndEffAll',99.9);
p00=prctile(rndEffAll',0.1);
med=median(rndEffAll');
loglog(nn,med,'k.-','LineWidth',2);
loglog([nn;nn],[p00;p99],'k')
loglog([nn*.93;nn*1.1],[p00;p00],'k')
loglog([nn*.93;nn*1.1],[p99;p99],'k')

axis([50 10^3*1.2 10^-1 10^1]);
set(gca,'XTick',nn);
set(gca,'LineWidth',2,'FontSize',12);

xlabel('Sequence length','FontSize',12);
ylabel('Efficiency','FontSize',12);
title(['No Overlap: n_e=',num2str(nVals),' p=',num2str(pEv),' n_h=',num2str(n2)]);
hold off;

set(gca,'Position',[.2,.15,.7,.7]);
set(gcf,'PaperPosition',[1,1,4,3]);
eval(['print -dpsc2 effNOV',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh',num2str(n),'n']);
disp(['print -dpsc2 effNOV',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh',num2str(n),'n']);


figure(2); clf;
semilogx(nn,mEffAll./max(rndEffAll'),'-ok','LineWidth',2); hold on;
semilogx(nn,mEffAll./med,':+k','LineWidth',2)
set(gca,'Position',[.2,.15,.7,.7]);
xlabel('Sequence length','FontSize',12);
ylabel('m-seq/random efficiency ratio','FontSize',12);
title(['No Overlap: n_e=',num2str(nVals),' p=',num2str(pEv),' n_h=',num2str(n2)]);
set(gca,'XTick',nn);
set(gca,'LineWidth',2,'FontSize',12);
axis([50 10^3*1.2 .9 2]);
eval(['print -dpsc2 effNOV',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh',num2str(n),'nRATIO']);
disp(['print -dpsc2 effNOV',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh',num2str(n),'nRATIO']);


eval(['save dataEffNOV',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh',num2str(n),'n rndEffAll p99 p00 med mEffAll']);
disp(['save dataEffNOV',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh',num2str(n),'n rndEffAll p99 p00 med mEffAll']);
