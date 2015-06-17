function [evtonsets,X,bestp,meff,origvec] = ga2(gensize,numgen,conditions,mspec,varargin)
%[evtonsets,X,bestp,meff,origvec] = ga2(gensize,numgen,conditions,mspec,[nooverlap])
%
% everything in conditions should be in number of frames (TRs), not seconds
% model_setup_ui.m automatically creates the correct format
%
% optional last argument is a 1 or 0.  1 Constrains the events so that
% events of different trials do not overlap.
%
% recommended: test conditions setup before running with
% [X,pv] = construct_model(mspec,conditions,[]);
% also [X,pv] = paramvec2onsets(mspec,conditions,[]);
%
% OR
%
% [X,pv] = construct_model(mspec,conditions,[],[1 14]);
% ...to constrain overlaps among trials to a min of 14 TRs apart,
% (variable trial onsets will only be properly constrained by the 2nd element).
% OR add in mspec.toverlap = [1 10], e.g.
%
% Tor Wager
% April 19, 2002
% Modified April 2, 2003: changed xover, convergence, constrained
%   so that min time between trial onsets is maintained
%   also added non-identity intrinsic autocorrelation
%   and capability for d-optimality
%
% Modified May 14, 2003: changed autocorrelation function to colored noise
%
% convergence criteria:
%   all efficiency scores are the same, OR no improvement in 100 generations
%
% defaults:
%   exponential autocorrelation function based on empirical estimates
%   saturation (clipping) of the predicted response at 2.5 times the original height
%       in order to approximate nonlinear effects
%   A-optimal design efficiency computation
%
%   Also possible: 
%   Combo A-D optimal design, equally weighted
%       for A-optimal design, mspec.dopt = 0; 
%       for D-optimality, specify mspec.dopt = 1;
%   But I'm still testing the D-optimality stuff.
%
% to see a test of A vs D optimality, try:
% Vi = eye(mspec.numframes);
% for i = 1:1000, 
%    [X,paramvec] = construct_model(mspec,conditions,[]); 
%    a(i) = calcefficiency(ones(1,size(X,2) - 1),[],pinv(X),Vi,0);
%    d(i) = 1 ./ (inv(det(X'*X)) .^ (1./size(X,2)));
% end
% d = d - mean(d);
% figure; plot(a,d,'x'); xlabel('A-efficiency'),ylabel('D-efficiency')
% a1 = calcefficiency(ones(1,size(X,2) - 1),[],pinv(X),Vi,0);
% d1 = 1 ./ (inv(det(X'*X)) .^ (1./size(X,2)));
% hold on;plot(a1,d1,'rs','MarkerFaceColor','r')


t0 = clock;


% --------------------------------------------------------------------
%
% * Setup and defaults
%
% --------------------------------------------------------------------

if length(varargin) > 0, nooverlap = varargin{1};,else, nooverlap = 0;,end


for j = 1:gensize, newvec{j} = [];,end

% intrinsic autocorrelation
% --------------------------------------------------------------------
%Vi = eye(mspec.numframes);      % use this for scanner white noise

% new autocorr, from VNL experiment (n = 10)
% this is invariant to the TR as well.
f = inline('5.6034 * exp(-.93 * x)','x');
t1 = f(1:mspec.TR:10); t1 = t1 ./ t1(1);

Vi = getv('make',t1,ceil(mspec.numframes));

% nonlinearity and fitness function
% --------------------------------------------------------------------

% nonlinear clipping threshold
satthresh = 2.5;

% display results vs. random?
showres = 1;

% converge if no change in this many generations
genconverge = 10;

% optimality - default is A
if ~isfield(mspec,'dopt'), mspec.dopt = 0;, end
if mspec.dopt == 2
    disp('Computing both A- and D-optimality'),
elseif mspec.dopt == 1
    disp('Computing D-optimality'), 
else
    disp('Computing A-optimality'),
end

if nooverlap & length(nooverlap) == 1
    % we need a list of all trial onset conditions to make sure we don't create
    % illegal lists in xover
    % older, more cumbersome way based on each condition
    trialonsetinfo = get_tonsetinfo(conditions);
elseif nooverlap & length(nooverlap) > 1
    % newer way: minimum time constraint between all trial onsets
    trialonsetinfo = nooverlap(2);
else
    % no constraints
    trialonsetinfo = -Inf;
end



% model seeding - m-sequences, 10 %
% works for designs with one subpart per trial - don't know about 
% others.  will have to insert more paramvec entries
% --------------------------------------------------------------------
nmseq = round(gensize ./ 10);
%[mseqs,newvec(1:nmseq)] = genMseq(mspec.numframes,[],[],nmseq,length(conditions));

[tmp,mseqs] = genMseq(mspec.numframes,[],[],nmseq,length(conditions));  % 10 sequences

for i = 1:length(mseqs)
    % get paramvec, which constrains length to correct values
    [tmp,pvm] = construct_model(mspec,conditions,mseqs{i});
    newvec{i} = pvm;
end

% set up contrasts across conditions
% --------------------------------------------------------------------
    numframes = [];,    % ne4eded for ga3power at the end, and for contrasts
    for i = 1:length(conditions),
        for j=1:length(conditions(i).subcond), numframes(end+1) = conditions(i).subcond(j).hrfest;,
        end
    end
    
if isfield(mspec,'contrasts')
    % executes the code to apply contrast, apply_dx_contrast.m, in
    % construct_model
    disp(['Found contrast to collapse columns: ']), mspec.contrasts
end

f1 = figure; set(gcf,'Color','w'); hold on; title('Max efficiency'),xlabel('Generation')


% matrix of multiplier values for each param in first cell to add noise
tmp = cat(1,newvec{:}); for i = 1:length(tmp), tmp{i} = diff(tmp{i}); end , tmp = cat(2, tmp{:});
noisestd = .2 .* range(tmp);
disp(['Noise std if no change is: ' num2str(noisestd)])

    
% --------------------------------------------------------------------
%
% * Start generations
%
% --------------------------------------------------------------------
        
        
for i = 1:numgen
    
    t1 = clock;
    
    for j = 1:gensize
    
        % --------------------------------------------------------------------
        % * make models
        % --------------------------------------------------------------------
        % if empty (1st generation), then it constructs a model.
        [X,paramvec{j}] = construct_model(mspec,conditions,newvec{j},nooverlap);
        X = modelSaturation(X,satthresh);
        
        % tmp for Emily
        %X(:,15*4+1:15*5) = X(:,15*4+1:15*5) + X(:,15*6+1:15*7); X(:,15*6+1:15*7) = [];
        
        XX{j} = X;
        
        if i == i & j == 1
            cweights = ones(1,size(X,2) - 1);
        end
        
        if mspec.dopt == 0
            eff(j) = calcEfficiency(cweights,[],pinv(X),Vi,0);
        elseif mspec.dopt == 1
            eff(j) = det(X'*X) .^ (1./size(X,2)); % 1 ./ (inv(det(X'*X)) .^ (1./size(X,2)));  % D-optimality
        elseif mspec.dopt == 2  % for A- D-optimality
            eff2(j) = calcEfficiency(cweights,[],pinv(X),Vi,0);
            eff(j) = det(X'*X) .^ (1./size(X,2));  % D-optimality
        end
        
    end
    
    if mspec.dopt == 2  % for A- D-optimality
        zeff = (eff - mean(eff)) ./ std(eff);
        zeff2 = (eff2 - mean(eff2)) ./ std(eff2);
        zeff = zeff + zeff2;
        meffi = max(zeff);
        meffind = find(zeff == meffi); meffind = meffind(1);
        meff(i) = eff(meffind);
        meff2(i) = eff2(meffind);
        t2 = clock;
        fprintf(1,'Gen %3.0f: D-eff is %3.4f, A-eff is %3.4f, model %3.0f, time = %3.2f, elapsed = %3.0f\n',i,meff(i),meff2(i),meffind,etime(t2,t1),etime(t2,t0));
        figure(f1); subplot 121; plot(1:i,meff,'r','LineWidth',2); title('D-efficiency')
        subplot 122; plot(1:i,meff2,'b','LineWidth',2); title('A-efficiency')
        drawnow
    
    else
        [meffi, meffind] = max(eff);
        meffind = meffind(1);
        meff(i) = meffi;
        t2 = clock;
        fprintf(1,'Gen %3.0f: best eff is %3.4f, model %3.0f, time = %3.2f, elapsed = %3.0f\n',i,meffi,meffind,etime(t2,t1),etime(t2,t0));
        figure(f1); plot(1:i,meff,'LineWidth',2); drawnow
    end
    
    % -------------------------------------------------------------------- 
    % add noise if no improvement at all
    % rounding and illegal values should be handled by construct_model
    % but need to be integers >= 0 for HRF power below
    % in next iteration
    % -------------------------------------------------------------------- 
    if meffind == 1

        for j = 1:length(paramvec)
            if j == meffind
                % leave the best one alone!
            else
                for k = 1:length(paramvec{j})
                    paramvec{j}{k} = paramvec{j}{k} + round( noisestd .* randn(size(paramvec{j}{k})) );

                    paramvec{j}{k}(paramvec{j}{k} < 0) = 0;
                end
            end
        end

    end

    % --------------------------------------------------------------------
    % * save best and crossover
    % -------------------------------------------------------------------- 
    bestp = paramvec{meffind};
    
    % test power of bestp
    clear delta
    for reg = 1:length(bestp) 
        delta{reg} = zeros(mspec.numframes,1);
        delta{reg}(bestp{reg}+1) = 1;
        delta{reg} = delta{reg}(1:mspec.numframes);
    end
    P = hrf_power(mspec.TR,conditions(1).subcond.hrfest,delta,eye(length(bestp)));
    bestipow(i) = max(P.ipow);
    figure(f1); hold on; plot(bestipow,'r','LineWidth',2),drawnow
    %figure;subplot(1,2,1);imagesc(P.X2); subplot(1,2,2);imagesc(XX{meffind}); colormap gray
    
    
    % --------------------------------------------------------------------
    % * Check for convergence
    % --------------------------------------------------------------------
    if i > genconverge
        last_gens = meff(i-genconverge+1:i);
        unique_in_last_gens = length(unique(last_gens));
        % if the above is one, there has been no change in last n
        % generations
    else
        unique_in_last_gens = Inf;
    end

    if (all(eff == median(eff))) || unique_in_last_gens == 1
        isconverged = 1;
        disp(['System converged at generation ' num2str(i)])

        break
    end
        
    % save 1st generation in origvec output
    if i == 1, origvec = newvec; end
    
    [newvec] = xover(paramvec,eff,trialonsetinfo);
        %newvec{ceil(rand * length(newvec))} = bestp;   % this causes lots of the same list to be inserted randomly
     newvec{1} = bestp;
     
end

% --------------------------------------------------------------------
% * End of GA
% --------------------------------------------------------------------
    
X = construct_model(mspec,conditions,bestp,nooverlap);
      
% --------------------------------------------------------------------
% * Figure of final model matrix
% --------------------------------------------------------------------
try
    figure('Color','w'); imagesc(X); colormap(gray)
    bestcondno = cond(X);
    tmp = corrcoef(X); tmp=tmp(:); tmp(isnan(tmp) | tmp == 1) = []; maxcorr = max(tmp);
    tmp = sprintf('Best model: Cond. # = %3.2f, max pairwise corr = %3.2f',bestcondno,maxcorr);
    title(tmp), 
catch
    disp('Problem with figures')
end

% this plots lines for onsets of first event type
%hold on; plot(repmat([.5 1],length(evtonsets{1}),1)',[1+(evtonsets{1} ./ mspec.TR) 1+(evtonsets{1} ./ mspec.TR)]','r','LineWidth',2)

% --------------------------------------------------------------------
% * Text list of onset times
% --------------------------------------------------------------------
try
    disp(' '); disp('GA2.m output')
    [evtonsets,paramvec] = paramvec2onsets(mspec,conditions,bestp,nooverlap);
catch
    disp('Error building evtonsets!  Incorrect model specification in model_setup_ui.m')
    evtonsets = [];
end

fname = ['ga2_output' datestr(now,31)]; fname(fname==' ')='_';fname(end-5)='h';fname(end-2)='m';fname=fname(1:end-2);
try
    eval(['save ' fname ' mspec conditions X bestp nooverlap evtonsets paramvec'])
catch
    warning('Not all vars created?  Saving everything.')
    eval(['save ' fname])
end

if showres
    
    % --------------------------------------------------------------------
    % * Figure of final model matrix
    % --------------------------------------------------------------------
    try
        [out,P] = ga3power(bestp,mspec,conditions,eye(length(numframes)));
    catch
        disp('Could not execute power plot.')
    end
    
    
    % --------------------------------------------------------------------
    % * Figure of A vs D optimality
    % --------------------------------------------------------------------
try
    figure('Color','w');hist(diff(evtonsets{1})); xlabel('Time between trial onsets'),ylabel('Frequency'); set(gca,'FontSize',18)
    
    
    if ~exist('Vi') == 1
        disp('No Vi detected: Using white noise')
        Vi = eye(mspec.numframes);      % use this for scanner white noise
    end
    
    disp('Comparing GA results (red square) to 1000 random designs (blue x)')
     for i = 1:1000, 
        [X,paramvec] = construct_model(mspec,conditions,[],nooverlap); 
        a(i) = calcefficiency(ones(1,size(X,2) - 1),[],pinv(X),Vi,0);
        d(i) = det(X'*X) .^ (1./size(X,2));     % 1 ./ (inv(det(X'*X)) .^ (1./size(X,2)));
     end
    figure('Color','w'); plot(a,d,'x'); set(gca,'FontSize',18); xlabel('A-efficiency'),ylabel('D-efficiency')
    [X,paramvec] = construct_model(mspec,conditions,bestp,nooverlap);
    a1 = calcefficiency(ones(1,size(X,2) - 1),[],pinv(X),Vi,0);
    d1 = 1 ./ (inv(det(X'*X)) .^ (1./size(X,2)));
    hold on;plot(a1,d1,'rs','MarkerFaceColor','r')
    title('GA result vs. random designs')
catch
    disp('Problem with A-vs-D optimality figure')
end

end

return




% --------------------------------------------------------------------
% * get best half, crossover to fill in lists
% -------------------------------------------------------------------- 
function newvec = xover(paramvec,eff,t)

% save best one - now done outside xover
%b = find(eff == max(eff)); b = b(1); b = paramvec{b};


% add a little random noise to efficiency - 1% of var of efficiency
eff = eff + randn(1,length(eff)) .* .01 * var(eff);

w = find(eff > median(eff));

% we can only do crossover if not all the designs are the same
if isempty(w)
    warning('Extremely homogenous sample!')
    % add a little random noise to efficiency - 1% of var of efficiency
    eff = eff + randn(1,length(eff)) .* .01 * mean(eff);
    w = find(eff > median(eff));
end

% check all paramvec values - works if all trial onsets are fixed
%for i=1:length(eff),ttmp=cat(1,paramvec{1}{:});ttmp=ttmp(:);,tmax(i)=max(ttmp);tmin(i)=min(ttmp);,end,[tmin; tmax]

% save best half
newvec = paramvec(w);
last = size(newvec,2);  % save size for adding random variation to existing

% fill in 2nd half with crossovers of first half
while length(newvec) < length(paramvec)
    
    w = ceil(rand(1,2) * length(newvec));
    babyv = [];
    
    for i = 1:length(newvec{w(1)})
        
        babyv{i} = rcomb(newvec{w(1)}{i},newvec{w(2)}{i});
        
        if length(t) > 1    % if we have full trialonsetinfo
        % older, more cumbersome way based on each condition
        % construct a set of legal vectors - cannot be less time between events than min onset time
        
        it = 1;
        while any(diff(babyv{i}) < t(1,i))
        
            tmp = [Inf diff(babyv{i})]; tmp = t(1,i) - tmp; tmp(tmp < 0) = 0;
            babyv{i} = babyv{i} + ceil(tmp);
            
            it = it + 1;
            if it > 20, 
                %warning(['Too many iterations to get a legal vector for ' num2str(i)]), 
                babyv{i} = rcomb(newvec{w(1)}{i},newvec{w(2)}{i});
            end
            
            if it > 100
                warning(['Too many iterations to get a legal vector for cond ' num2str(i)]), 
                w = ceil(rand(1,2) * length(newvec));
                babyv{i} = rcomb(newvec{w(1)}{i},newvec{w(2)}{i});
                keyboard
            end
        end
        
        else
            % newer way constrains min time only
            % t is constrained in construct_model, so leave out here
        end
            
    end
   
    newvec{end+1} = babyv;
    
end

% add random variation to existing best half
for i = 1:length(paramvec{1})   % i indexes events

    for j = 1:last  % j indexes organisms
        
        p = sign(randn(1,length(newvec{j}{i})));
        
        if length(t) > 1    % if old way
            
        if t(3,i) == 1
    
            p(1) = 0;
            d = diff(newvec{j}{i});   % get diffs to make sure we don't move something illegally
            d = [Inf d] - t(1,i);
            p(find(d <= 1 & p < 0)) = 0; % set shift value to 0 for these points
            newvec{j}{i} = newvec{j}{i} + p;
            
            while any(diff(newvec{j}{i}) < t(1,i))
        
                tmp = [Inf diff(newvec{j}{i})]; tmp = t(1,i) - tmp; tmp(tmp < 0) = 0;
                newvec{j}{i} = newvec{j}{i} + ceil(tmp);
            
            end
            
        elseif t(3,i) < 1
            
            tmp2 = newvec{j}{i} + p;
            %disp('orig')
            %tmp2
            % ensure no values go outside proscribed range
            tmp2(tmp2 < t(1,i)) = t(1,i);
            tmp2(tmp2 > t(2,i)) = t(2,i);
            %disp('final')
            %tmp2
            
            newvec{j}{i} = tmp2;
        else
            warning('This should not happen.  t(3,:) should be 1 or -Inf')
        end
        
        else   % not old way, no constraints here (constraints in construct_model) 
            newvec{j}{i} = newvec{j}{i} + p;
            newvec{j}{i}(newvec{j}{i} < 0) = 0;
        end
            
    end
    
    % warning: may not really ensure that we have stimuli too close
end
 
%newvec{1} = b;      % re-insert best one

return


% --------------------------------------------------------------------
% * crossover
% -------------------------------------------------------------------- 
function c = rcomb(a,b)
% combines 2 vectors of numbers at a random crossover point

w = ceil(rand * (length(a) - 1)); 
c(1:w) = a(1:w);
c = [c b(w+1:end)];

return




function [t,sc] = get_tonsetinfo(conditions)
% get list of onset times to match paramvec
% this contains info needed to see if one trial is impinging on the previous one, based on the length
% so it saves a vector, t, of the minimum times in TRs that each paramvec can have
% this value is compared, for each onset vector in paramvec, to the difference among onset times (trial length)
% to make sure no trial length is less than the minimum.
%
% t has first row = min, 2nd = max, 3rd = trial onsets (1) or subpart (0)

ind = 1;
for i = 1:length(conditions)
    if size(conditions(i).onsets,2) == 2   % variable range of trial onsets
        tmp = sum(conditions(i).onsets);
        t(1,ind) = tmp(1);
        t(2,ind) = tmp(2);
        t(3,ind) = 1;
        ind = ind + 1;
    end
    
    for j = 1:length(conditions(i).subcond)
        if length(conditions(i).subcond(j).onsets) > 1 & any(diff(conditions(i).subcond(j).onsets))
            t(1:2,ind) = [min(conditions(i).subcond(j).onsets); max(conditions(i).subcond(j).onsets)];      
            t(3,ind) = 0;
            % need this for subcond, i think   %conditions(i).subcond(j).onsets(1);
            ind = ind + 1;
        end
    end
end

t(t == 0) = -Inf;

return


