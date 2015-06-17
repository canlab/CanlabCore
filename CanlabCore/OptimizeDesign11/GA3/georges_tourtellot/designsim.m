function [t,conse,allt,glm] = designsim(noisetype,dmodel,HRF,ISI,TR,noise_var,c,beta,S,varargin)
% function [t,modelse,allt,glm] = designsim(noisetype,dmodel,HRF,ISI,TR,noise_var,c,betas,S,xcfunction [if using 'myxc'],num2avgover)
%  [t,modelse,allt,glm] = designsim('myxc',M.modelatTR,HRF,ISI,TR,noise_var,c,beta,fullS,xc)
%
% noisetype:
%   '1overf' : use Luis Hernandez' 1/f model with your specified noise variance
%   'myxc'   : use your own autocorrelation function and noise variance
%
% model = matrix of regressors (cols) sampled at frequency of ISI
% noise_var = variance of 1/f noise from your scanner
% c = a contrast across the regressors
% beta = a list of betas to test
%   each beta should be a row vector, and you can have multiple rows
%
% S: smoothing matrix, or 0 for no smoothing.
% _______________________________________________________________________
% Output is a distribution of t-scores in a vector and the maximum correlation between predictors
% t-scores are in a column vector, one for each set of betas.
%
% Special mode: if 4 output arguments specified, 'plot mode'.  Produces graph and extended output.
% outputs:
%	t 		mean of t-scores for specified number of noise vectors tested, numnoise in script
%	allt	vector of all t values generated for each noise vector
%
% specifying a 4th output puts designsim.m into 'verbose' mode.
%
% Last modified   5/24/01 Tor Wager		add noise first, then smooth and filter, etc.



% if no stimuli for one condition, gives error - create extra col of zeros if this happens.
while size(beta,2) > size(dmodel,2)
 	warning('      Designsim.m : dmodel has too few columns or beta is too long! No stim. in one cond.?')    
	whos dmodel
     	dmodel = [dmodel zeros(size(dmodel,1),1)];   
end

if isempty(c)
	if nargout > 3,disp('	...designsim.m: No contrasts found.  Using first predictor for one sample T-test.'),end
	c = zeros(1,size(dmodel,2));
	c(1) = 1;
end

c = c';

done = 0;
%if nargout > 3,numnoise = 1;,else numnoise = 1;,end
% if average over noise models is given, return the average - otherwise return each individual t and se score.
if nargin > 10, numnoise = varargin{2};,else numnoise = 1;,end
if nargin > 11, NULL_flag = varargin{3};,else NULL_flag = 0;,end  %gst

% make the noise vectors
% =====================================================================================
switch noisetype
case '1overf'
   noise = make1overf(size(dmodel,1), 1/ISI) * sqrt(noise_var);
case 'myxc'

%varargin{1}    %gst
    
    for i = 1:numnoise
        noise(:,i) = noisevector(size(dmodel,1),varargin{1},noise_var)';
    end
otherwise error('unsupported noise type.')
end


   		% Create the response data by weighting all the regressors by a beta parameter
   		% adding all the regressors together
   		% and adding noise to the result  
 beta = beta';
 if nargout > 3,
    	%figure;,
    	%disp('testing with only one noise model, returning 1st beta vector only in glm.')
    	%beta = beta(:,1); noise = noise(:,1);
	c = c(:,1);	% test first contrast here - this should be changed back for normal use. this is for rundesignsimMAP
 else c = c(:,1); % tests only first contrast if you specify many iterations.
 end


% make the data vectors and fit
% =====================================================================================
for i = 1:size(beta,2)
    
	if nargout > 3
        	subplot(4,size(beta,2),3*i-2);hold on
        	plot(dmodel * beta(:,i),'r');title(['ideal model, betas = ' num2str(beta(:,i)')]);drawnow
        	clear glm
	end
 
	% make data
	% --------------------------------------------------	
	for j = 1:size(noise,2)

		try
			myfit = dmodel * beta(:,i);
		catch
			disp(['	designsim.m: dmodel and beta are not the same size - wrong beta size?'])
			whos dmodel
			beta
			error('exiting.')
		end
%        myfit = (myfit - mean(myfit));		% mean center so betas are not 'correlated' with intercept.
        myfit = (myfit - ones(size(myfit,1),1)*mean(myfit));  %gst
%        data(:,j) = myfit + noise(:,j);	%gst

        if (NULL_flag ~= 0)   %gst     NULL_flag is the myfit column vect corresponding to NULL condition
            if (NULL_flag>size(myfit,2))
    			disp(['NULL_flag (the column of the design matrix to treat as rest) is too large'])
                return;
            end
            myfit(:,NULL_flag)=0;
        end
        myfit=sum(myfit, 2);                %gst
        data(:,j) = myfit + noise(:,j)*ones(1,size(myfit,2));	  %gst
        
		% testing stuff -
		%disp(['designsim.m: beta = '])
		%beta(:,i)
		%figure;subplot(3,1,1);plot(noise(:,j))
		%subplot(3,1,2);plot(dmodel * beta(:,i))
		%subplot(3,1,3);plot(data(:,j));	
	end
	
	% smoothing
	% --------------------------------------------------
	if S, 
		dmodel = S * dmodel;
		data = S * data;
		glm.smoothed = 'yes';
	end
    
	% fit model to data and save t and se
	% --------------------------------------------------	
	for j = 1:size(noise,2)
        	[myt,myconse] = my_glm(dmodel,data(:,j),c);		% fit model
        	allt(i,j) = myt(1);
		allse(i,j) = myconse(1);
        
		% make glm structure if necessary
		% ---------------------------------------------
        	if nargout > 3
            		glm.y = data;
            		glm.X = [dmodel ones(size(dmodel,1),1)];
            		glm.betas = (glm.X \ glm.y)';
            		glm.fit = glm.X * glm.betas';
            		glm.e = glm.y - glm.fit;
            		glm.xtxi = inv(glm.X' * glm.X); 
            		glm.df = (size(glm.y,1) - size(glm.X,2));
            		glm.evar = glm.e' * glm.e / glm.df;
            		glm.se = sqrt(diag(glm.xtxi .* glm.evar))';
            		glm.t = (glm.betas ./ glm.se);
            		glm.p = tdist(glm.t,glm.df);
%            		glm.c = [c;zeros(1,size(c,2))];   %gst
            		glm.c = c;   %gst                   
            		glm.con_se = diag(sqrt(glm.evar .* (glm.c' * glm.xtxi * glm.c)))';
            		glm.con_t = glm.betas * glm.c ./ glm.con_se;
            		glm.con_p = tdist(glm.con_t,glm.df);

	    		if j == 1							% plot the first noise model.
        			subplot(4,size(beta,2),3*i-1);hold on
        			plot(data(:,j));hold on;plot(glm.fit,'r--'); title(['Data with noise added and fit']);drawnow
        			subplot(4,size(beta,2),3*i);hold on
        			plot(data(:,j) - glm.fit); title('Residuals');drawnow
 	    		end
        	end
        
	end
	% loop thru noise models.
 end	% loop thru betas.

t = mean(allt,2);
conse = mean(allse,2);


return



% rmat = abs(corrcoef(dmodel));%rmat(rmat == 1) = 0;		% zero diagonals%r = max(max(rmat));




