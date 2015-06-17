% efficiencyOFexpDesignsCorrMtrx
% gtb: event matrices are defined based on vectors with the average
% energy removed

loadDataYes=1; % =1 -> loads existing data file

cutEndPad=1; %cutting reduces efficiency!
collectEvMtrx=1;
% for simply randomized designs:
pwr=9; %9
ms=m2bin([mseq(2,pwr,0,1)])'; %
n=length(ms);  % scan duration in TRs

plotMeanYes=1;

n2=12;
nVals=3;%1,2,3,4 or 8 
nOnes=0.5;
overlapYes=1;
nReps=1000; %n-2;

b=[0.406;0.8825]; % parameters from SPG fMRI data
rndEvent=balancedRnd(n,nVals,nOnes,overlapYes);
eventMatrix=makeEventMtrx(rndEvent,n2);
fittedACorr=autocorrFnct(b,1:n+(n2-1)*rem(cutEndPad+1,2));
fittedACorr=fittedACorr/sqrt(sum(fittedACorr.^2));
Cninv=inv(toeplitz(fittedACorr));

clear ER;
ERrndEff=zeros(nReps,1);
ERrndEffC=zeros(nReps,1);
if collectEvMtrx, ER.eventMatrix=[]; end;

if ~loadDataYes
tic
for k=1:nReps;

  %nOnes=k/(n-1);
   rndEvent=balancedRnd(n,nVals,nOnes,overlapYes);
   eventMatrix=makeEventMtrx(rndEvent,n2);
   if cutEndPad, eventMatrix=eventMatrix(1:n,:); end;
   eventMatrix=eventMatrix-ones(size(eventMatrix,1),1)* ...
     sum(eventMatrix)/size(eventMatrix,1);
   
   if collectEvMtrx,
     ER(k).eventMatrix=eventMatrix;
   end;  
   designEff=1/trace(inv(eventMatrix'*eventMatrix));
   ERrndEff(k)=designEff;
   designEff=1/trace(inv(eventMatrix'*Cninv*eventMatrix));
   ERrndEffC(k)=designEff;
end;
rndTime=toc


maxREffCind=find(max([ERrndEffC])==[ERrndEffC]);
minREffCind=find(min([ERrndEffC])==[ERrndEffC]);
maxREffind=find(max([ERrndEff])==[ERrndEff]);
minREffind=find(min([ERrndEff])==[ERrndEff]);

if collectEvMtrx,
  figure(2); clf;
  subplot(2,2,1);
  imagesc(ER(maxREffCind).eventMatrix)
  title('max Eff with Corr')
  subplot(2,2,2);
  imagesc(ER(maxREffind).eventMatrix)
  title('max Eff')
  subplot(2,2,3);
  imagesc(ER(minREffCind).eventMatrix)
  title('min Eff with Corr')
  subplot(2,2,4);
  imagesc(ER(minREffind).eventMatrix)
  title('min Eff')
end;


mEffAll=[];
mEffCAll=[];

shift1=floor(1*n/nVals);
shift2=floor(2*n/nVals);
shift3=floor(3*n/nVals);
shift4=floor(4*n/nVals);
shift5=floor(5*n/nVals);
shift6=floor(6*n/nVals);
shift7=floor(7*n/nVals);

for cycle=1:floor(n/nVals)
 mEvent=ms; 
 switch nVals,
  case 1
    mEvent=ms; 
  case 2
    mEvent=ms; 
    mEvent=[mEvent; [ms(shift1+cycle+1:end) ms(1:shift1+cycle)]];
  case 3
    mEvent=ms; 
    mEvent=[mEvent; [ms(shift1+cycle+1:end) ms(1:shift1+cycle)]];
    mEvent=[mEvent; [ms(shift2+1:end) ms(1:shift2)]];
  case 4   
    mEvent=ms; 
    mEvent=[mEvent; [ms(shift1+cycle+1:end) ms(1:shift1+cycle)]];
    mEvent=[mEvent; [ms(shift2+cycle+1:end) ms(1:shift2+cycle)]];
    mEvent=[mEvent; [ms(shift3+cycle+1:end) ms(1:shift3+cycle)]];
  case 7   
    mEvent=ms; 
    mEvent=[mEvent; [ms(shift1+cycle+1:end) ms(1:shift1+cycle)]];
    mEvent=[mEvent; [ms(shift2+cycle+1:end) ms(1:shift2+cycle)]];
    mEvent=[mEvent; [ms(shift3+cycle+1:end) ms(1:shift3+cycle)]];
    mEvent=[mEvent; [ms(shift4+cycle+1:end) ms(1:shift4+cycle)]];
    mEvent=[mEvent; [ms(shift5+cycle+1:end) ms(1:shift5+cycle)]];
    mEvent=[mEvent; [ms(shift6+cycle+1:end) ms(1:shift6+cycle)]];
  case 8   
    mEvent=ms; 
    mEvent=[mEvent; [ms(shift1+cycle+1:end) ms(1:shift1+cycle)]];
    mEvent=[mEvent; [ms(shift2+cycle+1:end) ms(1:shift2+cycle)]];
    mEvent=[mEvent; [ms(shift3+cycle+1:end) ms(1:shift3+cycle)]];
    mEvent=[mEvent; [ms(shift4+cycle+1:end) ms(1:shift4+cycle)]];
    mEvent=[mEvent; [ms(shift5+cycle+1:end) ms(1:shift5+cycle)]];
    mEvent=[mEvent; [ms(shift6+cycle+1:end) ms(1:shift6+cycle)]];
    mEvent=[mEvent; [ms(shift7+cycle+1:end) ms(1:shift7+cycle)]];
  otherwise
   error(['wrong # of events - in m-seq']);
 end;
 
  eventMatrix=makeEventMtrx(mEvent,n2);
  if cutEndPad, eventMatrix=eventMatrix(1:n,:); end;
  eventMatrix=eventMatrix-ones(size(eventMatrix,1),1)* ...
     sum(eventMatrix)/size(eventMatrix,1);
  
  mEff=1/trace(inv(eventMatrix'*eventMatrix));
  mEffC=1/trace(inv(eventMatrix'*Cninv*eventMatrix));

 mEffAll=[mEffAll mEff];
 mEffCAll=[mEffCAll mEffC];
end;
else 
  eval(['load dataEffVSeffC',num2str(nVals),'ev',num2str(nReps),'reps',num2str(n2),'nh',num2str(n),'n']);
end;


ind=find(max(mEffAll)==mEffAll); ind=ind(1);
indC=find(max(mEffCAll)==mEffCAll);indC=indC(1);

figure(1); clf
rndMean=mean([ERrndEff]);
rndMeanC=mean([ERrndEffC]);
maxRndEff=rndMean; %max([ERrndEff]);
maxRndEffC=rndMeanC; %max([ERrndEffC]);

plot([ERrndEff]/rndMean,[ERrndEffC]/rndMeanC,'o','MarkerSize',1);hold on;
plot(mEffAll(ind)/rndMean,mEffCAll(ind)/rndMeanC,'r*');
plot(mEffAll(indC)/rndMean,mEffCAll(indC)/rndMeanC,'r+');

if plotMeanYes,
  rndStd=std([ERrndEff]/rndMean);
  rndStdC=std([ERrndEffC]/rndMeanC);

  %plot(rndMean,rndMeanC,'ob','LineWidth',2);
  plot([1,1+rndStd],[1,1],'w','LineWidth',2); hold on;
  plot([1,1],[1,1+rndStdC],'w','LineWidth',2);
  
end;

title(['n=',num2str(n),' n_e=',num2str(nVals),' n_h=',num2str(n2)],'FontSize',12)
xlabel('Efficiency','FontSize',12)
ylabel('Efficiency with noise','FontSize',12)
minLim=min([ERrndEff]/rndMean);
maxLim=max([[ERrndEff]/rndMean; mEffAll(ind)/rndMean]);
minLimC=min([ERrndEffC]/rndMeanC);
maxLimC=max([[ERrndEffC]/rndMeanC; mEffCAll(indC)/rndMeanC]);
axis([0.9*minLim 1.1*maxLim 0.9*minLimC 1.1*maxLimC]);
axis equal
set(gca,'Position',[.2,.2,.6,.6]);
set(gcf,'PaperPosition',[1,1,3,3]);

eval(['print -dpsc2 effVSeffC',num2str(nVals),'ev',num2str(n2),'nh',num2str(n),'n']);
disp(['print -dpsc2 effVSeffC',num2str(nVals),'ev',num2str(n2),'nh',num2str(n),'n']);

disp(['save dataEffVSeffC',num2str(nVals),'ev',num2str(nReps), ...
      'reps',num2str(n2),'nh',num2str(n),'n ERrndEff ERrndEffC mEffAll mEffCAll']);
eval(['save dataEffVSeffC',num2str(nVals),'ev',num2str(nReps), ...
      'reps',num2str(n2),'nh',num2str(n),'n ERrndEff ERrndEffC mEffAll mEffCAll']);

