function ModelSummary = summary(Designs,Models)

%function ModelSummary = summary(Designs,Models,conditions,freqConditions,ISI,TR,cbalColinPowerWeights,lowerLimit,maxOrder)

% ModelSummary = summary(Designs,Models)

% model summary



% get default values from gui or script

defaults_optimize



H = findobj('Tag', 'E1');                                                                                                                                 

conditions = str2num(get(H, 'String'));     

H = findobj('Tag', 'E2');                                                                                                                                 

freqConditions = str2num(get(H, 'String')); 

H = findobj('Tag', 'E4');                                                                                                                                 

ISI = str2num(get(H, 'String'));  

H = findobj('Tag', 'E5');                                                                                                                                 

TR = str2num(get(H, 'String')); 

H = findobj('Tag', 'E6');                                                                                                                                 

cbalColinPowerWeights = str2num(get(H, 'String'));

H = findobj('Tag', 'E9');                                                                                                                                 

lowerLimit = str2num(get(H, 'String'));  

H = findobj('Tag', 'E10');                                                                                                                                 

maxOrder = str2num(get(H, 'String')); 





if not(size(Designs,2) == size(Models,3))

   error('number of designs and models does not match.')

else

   nmod = size(Designs,2);

end



ModelSummary = zeros(nmod,6);



for nm = 1:nmod

   model = Models(:,:,nm);

   stimList = Designs(:,nm);

   

   disp(['Computing model ' num2str(nm)])

   

   ModelSummary(nm,1) = getCounterBal(stimList, maxOrder,conditions,freqConditions);

   

   colin = abs(corrcoef(model));

   colin(colin == 1) = 0;		% zero diagonals

   colin = max(max(colin));	% RETURNS MAXIMUM CORRELATION COEF BTWN 2 PREDS AS FITNESS 

   ModelSummary(nm,2) = colin;

      

	ModelSummary(nm,3) = getPower(model, TR, ISI, lowerLimit,0);   



   %ModelSummary(nm,4) = cbalColinPowerWeights(1:3) * [ModelSummary(nm,1) ModelSummary(nm,2) ModelSummary(nm,3)]';

   

   for i = 1:niterations

      % Test idealized GA - no response variability

   	[t(i),r(i)] = designsim(model,ISI,noise_var,c);

      if mod(i,50) == 0, disp(['done ' num2str(i)]),end;

   end

   ModelSummary(nm,5) = mean(t);

   ModelSummary(nm,6) = mean(r);

   clear t,r;

end



% Determine overall fitness

% column vector of z scores for cbal, colin, power

clear z;

for i = 1:3

z(:,i) = (ModelSummary(:,i) - mean(ModelSummary(:,i))) / std(ModelSummary(:,i));

end

z = z'



												% weighted

   ModelSummary(:,4) = (cbalColinPowerWeights(1:3) * z)'          





disp('cbal colin power overall avg_t avg_r')









