function model = modifiedconv_wager2005(tr, indicator_vec, varargin)
% Modified convolution adjusting height, time-to-peak, and dispersion as a function of stumulus history
%
% model = modifiedconv_wager2005(tr,indicator_vec,heighteq [all opt],delayeq,ttopeakeq,uonseteq)
%
% 06/20/01 Tor Wager
%
% tr = repetition time (sampling rate) of scanning, in seconds
% indicator_vec = condition function
%         an indicator vector of zeros and ones, where ones indicate event
%         onsets
%
% variable input functions:
% field followed by equation for scaling
% fields are:
% 'height', 'delay', 'ttopeak', 'uonset'
% height, onset delay, time to peak, undershoot onset delay
% equations are inline functions, as defined by default in the body of this
% function.
%
% NOTE:
% the default equations are only valid for 1 s interstimulus intervals at this point
% AND
% this function gives an approximation based on the number of events
% occurring w/i the last 12 seconds, so will need work to improve accuracy for
% building real jittered event-related designs.
%
% example:
% indicator_vec = [1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]';
% X = modifiedconv(2,indicator_vec);
% X is convolved predictor
%plot(X)
%X2 = conv(indicator_vec,spm_hrf(2)./max(spm_hrf(2)));
%hold on; plot(X2,'r'); legend({'Modified' 'Linear'})
%
% Please see:
% Wager, T. D., Hernandez, L., Vasquez, A., Nichols, T., and Noll, D.
% C. (2005). Accounting for nonlinear BOLD effects in fMRI: Parameter
% estimates and model for accurate prediction in variable-duration blocked
% and rapid event-related studies.  Neuroimage.

heighteq = []; delayeq = []; peakeq = []; uonseteq = [];

% ---------------------------------------------------------------------
% * defaults
% ---------------------------------------------------------------------
height = 1; delay = 0; peak = 6; uonset = 16;
disp = 1; udisp = 1; rtou = 6; klength = 32;

p = [peak uonset disp udisp rtou delay klength];

hrf = spm_hrf(tr,p)./ max(spm_hrf(tr,p));

% These equations use parameters estimated from Wager et al. 2005
heighteq = inline('1.7141.*(exp(-2.1038.*x)) + 0.4932.*(exp(-0.0770.*x))');
delayeq = inline('-13.4097.*(exp(-1.0746.*x)) + 4.8733.*(exp(-0.1979.*x))');
peakeq = inline('37.5445.*(exp(-2.6760.*x)) + -3.2046.*(exp(-0.2120.*x)) + 5.6344');

% ---------------------------------------------------------------------
% * set up arguments
% ---------------------------------------------------------------------

for ind = 1:2:length(varargin)
    if strcmp(varargin{ind},'height')
        heighteq = varargin{ind+1};
        %disp('modifiedconv.m  - height scaling equation is:')
        heighteq
    elseif strcmp(varargin{ind},'delay')
        delayeq = varargin{ind+1};
        %disp('modifiedconv.m  - onset delay scaling equation is:')
        delayeq
    elseif strcmp(varargin{ind},'ttopeak')
        peakeq = varargin{ind+1};
        %disp('modifiedconv.m  - time to peak scaling equation is:')
        peakeq
    elseif strcmp(varargin{ind},'uonset')
        uonseteq = varargin{ind+1};
        %disp('modifiedconv.m  - undershoot onset delay scaling equation is:')
        uonseteq
    end
end

myzeros = zeros(length(indicator_vec),1);
mylen = length(myzeros);
numels = 12 ./ tr;  % number of elements to count

%figure;hold on

for i = 1:max(indicator_vec)
    % For each integer in indicator_vec
    
    reg = myzeros;
    
    for j = 1:length(indicator_vec)
        % for each time point
        
        if indicator_vec(j) == i
            
            trialdelta = myzeros;
            trialdelta(j) = 1;
            trialp = p;
            
            % figure out how many of same type came before
            % 1 is "first stim in sequence"
            
            % need to have a weighted function that falls off with time
            
            %myc = indicator_vec(j-min(29,j)+1:j)
            %myw = timeweights(end-length(myc)+1:end)
            
            myc = indicator_vec(j-min(numels,j)+1:j);
            stimpos = sum(myc == i);
            tmp = find(myc);
            dens = stimpos ./ length(myc(tmp(1):end));  % empty to filled within last stimpos positions
            % this is a scaling factor for stimulus position based on history density
            % if prior history is not full of stimulation,
            % the 'position' estimate is shifted toward 1, but cannot go below 1
            
            %max(1,stimpos .* dens)
            
            if ~isempty(heighteq)
                height = heighteq(max(1,stimpos .* dens));
            end
            
            if ~isempty(delayeq)
                trialp(6) = delayeq(max(1,stimpos .* dens));
            end
            
            if ~isempty(peakeq)
                trialp(1) = peakeq(max(1,stimpos .* dens));
            end
            
            if ~isempty(uonseteq)
                trialp(2) = uonseteq(stimpos);
            end
            
            trialhrf = spm_hrf(tr,trialp);
            trialhrf = height .* (trialhrf ./ max(trialhrf));
            % fprintf('%3.2f ',height)
            
            mytrialpred = conv(trialhrf,trialdelta);
            mytrialpred = mytrialpred(1:mylen);
            
            reg = reg + mytrialpred;
            
            %subplot 121, plot(trialhrf,'Color',rand(1,3)),set(gca,'YLim',[0 1]),title(['cond ' num2str(i) ' evt ' num2str(j)]),pause(.5),
            %subplot 122, plot(reg),set(gca,'YLim',[0 2])
            
        end  % time point
    end % integer
    
    % normalize so that height of unit response (stimpos == 1) is 1
    if ~isempty(heighteq)
        height = heighteq(1);
        model(:,i) = reg ./ height;
    else
        model(:,i) = reg;
    end
    
end


end
