function reg = modifiedconv(tr,condf,varargin)
% :Usage:
% ::
%
%     model = modifiedconv(tr,condf,heighteq [all opt],delayeq,ttopeakeq,uonseteq)
%
% :Inputs:
%
%   **tr:**
%        repetition time (sampling rate) of scanning, in seconds
%
%   **condf:**
%        condition function
%         an indicator vector of zeros and ones, where ones indicate event
%         onsets
%
% USES nonlinear saturation in height only
% with a guess as to what the decrease in saturation is as a function of
% the time since previous stimulation (exponential model, alpha version)
% 
% :Examples:
% ::
%
%    condf = [1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]';
%    X = modifiedconv(2,condf);
%    % X is convolved predictor
%    plot(X)
%    X2 = conv(condf,spm_hrf(2)./max(spm_hrf(2)));
%    hold on;
%    plot(X2,'r');
%    legend({'Modified' 'Linear'})
%
% :Please see:
% Wager, T. D., Hernandez, L., Vasquez, A., Nichols, T., and Noll, D.
% C. (in press). Accounting for nonlinear BOLD effects in fMRI: Parameter 
% estimates and model for accurate prediction in variable-duration blocked 
% and rapid event-related studies.  Neuroimage.
%
% ..
%    06/20/01 Tor Wager
% ..

%heighteq = []; delayeq = []; peakeq = []; uonseteq = [];

% ---------------------------------------------------------------------
% * defaults
% ---------------------------------------------------------------------
height = 1; delay = 0; peak = 6; uonset = 16;
dispers = 1; udisp = 1; rtou = 6; klength = 32;
rtou = Inf; % get rid of undershoot

p = [peak uonset dispers udisp rtou delay klength];

hrf = spm_hrf(tr,p)./ max(spm_hrf(tr,p));

%heighteq = inline('1.7141.*(exp(-2.1038.*x)) + 0.4932.*(exp(-0.0770.*x))');
%delayeq = inline('-13.4097.*(exp(-1.0746.*x)) + 4.8733.*(exp(-0.1979.*x))');
%peakeq = inline('37.5445.*(exp(-2.6760.*x)) + -3.2046.*(exp(-0.2120.*x)) + 5.6344');%

% new idea: model nonlinear saturation response as a function of 
% a fixed saturation response (that we know is e.g., 1-(.6658 x original) with a
% stimulus 1 s before) x an exponential discount factor for how LONG AGO the
% previous stimulus occurred
% SO THAT: Rnow = Ro - sum(s(t-to)gamm(t-to)
% gam(1) should be 1-.6658, and area under gam should be no more than 1,
% because you can never lose more than 100% of the signal
% BUT we may end up having to take interactions betwn stim into account,
% which is harder: e.g., a stim occurs .5 s before and another occurs .3
% sec before.  The .3 one is going to produce a smaller response and thus
% less saturation, and you can never saturate more than 1 (100%) even w/
% many stimuli occurring in a short time frame (w/i 1 s)
% SO this model, now, can predict negative response values with more than
% 100% sat.  

gam = inline('(1-.6658) * (1./exp(-a))*exp(-a*t)','t','a');
a = .7;  % let's assume this is the exp for now. gives 70% sat w/hist 1:5, 74% w/hist 1:10 


myzeros = zeros(length(condf),1);
mylen = length(myzeros);
numels = 12 ./ tr;  % number of elements to count

reg = myzeros;

whstim = find(condf);   % which elements contain indicators of stimulation

if any(condf(whstim)>1), 
    disp('warning: modifiedconv not valid for onset mag > 1, which you appear to have entered.'), 
end
    
for i = 1:length(whstim)
    
    j = whstim(i);  % index of which element
    
    % for each element, get predicted height

			trialdelta = myzeros;
			trialdelta(j) = 1;
			trialp = p;
			
			% figure out how many of same type came before
            % 1 is "first stim in sequence"
            
            %myc = condf(j-min(29,j)+1:j)
            %myw = timeweights(end-length(myc)+1:end)
            
            myc = condf(j-min(numels,j)+1:j);   % recent events, including current stim
            
            wh = find(myc(1:end-1));            % which elements
            
            times = length(myc) - wh;  % times at which these occurred in elements
            times = times * tr;        % convert to seconds
            
            s = myc(wh);               % stimulus intensity at each time, probably 1 for typical indicator vector
            
            sat = sum(s .* gam(times,a));  % total saturation
            
            height = max(0,1 - sat);           % magnitude of this event
            
            trialdelta = myzeros;
            trialdelta(j) = condf(j);
            
			mytrialpred = conv(height.*hrf,trialdelta);
			mytrialpred = mytrialpred(1:mylen);

			reg = reg + mytrialpred;
            
		
end


return
