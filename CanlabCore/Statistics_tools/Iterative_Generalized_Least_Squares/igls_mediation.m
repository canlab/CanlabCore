function out = igls_mediation(y, x, m, varargin)
    % function [beta, betastar, iterations, elapsed_time] = igls(y, x, m, [optional arguments])
    %
    % Variance Component Estimation using IGLS/RIGLS
    %
    % y = d1 + cpx + bm + epsilon_y  where  epsilon_y ~ N(0,sigma_y*Vy)
    % 
    % m = d0 + ax + epsilon_m  where  epsilon_m ~ N(0,sigma_m*Vm)
    %
    % a ~ N(0, sigma_a), b ~ N(0, sigma_b) and cp ~ N(0, sigma_cp)
    %
    % Calculate d1, a, b, cp, sigma1, sigma2, sigma_a, sigma_b and sigma_c using Maximum Likelihood
    % methods (IGLS) and Restricted Maximum Likelihood methods (RIGLS).
    %
    % This program is based on methods which are described in the following papers:
    %
    % Goldstein, H. (1986). Multilevel mixed linear model analysis using iterative generalized least squares. Biometrika 73, 43-56.
    % Goldstein, H. (1989). Restricted unbiased iterative generalized least-squares estimation, Biometrika 76, 622-623.
    %
    % Inputs:
    %
    % y - matrix T x subjects
    % x - matrix T x subjects
    % m - matrix T x subjects
    %
    % Optional inputs:
    %
    % 'noverbose'   suppress verbose output
    % 'iter'        max number of iterations
    % 'type'        'i' for igls (default) or 'r' for rigls
    % 'eps'         epsilon for convergence : all changes < (epsilon * beta)
    % 'ar'          order of AR(p) process; default is 0 (no AR model)
    % 'within_var'  specify common variance within subjects
    % 
    % Outputs:
    %
    % beta - vector of length  5. Contains estimates of d1, d0, c, a and b.
    % betastar - vector of length 8. Contains estimates of sigma_d0^2, sigma_d1^2, sigma_a^2, sigma_c^2, sigma_b^2, sigma_ab^2, sigma_ac^2 and sigma_bc^2.
    % Sigma_y - vector of length sub. Contains an estimate of sigma_y for each subject.
    % Cov_beta - matrix 2 x 2. Contains covariance matrix for beta.
    % Cov_betastar - matrix 2 x 2. Contains covariance matrix for betastar.
    % iterations - number of iterations performed
    % elapsed_time - amount of time needed to run program
    %
    % By Martin Lindquist, April 2007
    %
    % Example: Create simulated data and test
    % ----------------------------------------------------------------
    % len = 200; sub = 20;
    % x = zeros(len,sub);
    % x(11:20,:) = 2;                   % create signal
    %
    % a = normrnd(0.6,0,sub,1);       % slope between-subjects variations
    % d0 = normrnd(2,0,sub,1);         % intercept between-subjects variations
    %
    % b = normrnd(0.3,0,sub,1);       % slope between-subjects variations
    % c = normrnd(0.5,0,sub,1);       % slope between-subjects variations
    % d1 = normrnd(3,0,sub,1);         % intercept between-subjects variations
    %
    % % Create y: Add between-subjects error (random effects) and measurement noise
    % % (within-subjects error)
    % for i=1:sub, 
    %   m(:,i) = d0(i) + a(i).*x(:,i) + normrnd(0,0.3,len,1);
    %   y(:,i) = d1(i) + c(i).*x(:,i) + b(i).*m(:,i) + normrnd(0,0.5,len,1);
    % end;
    %
    % out = igls_mediation(y, x, m);
    
    c1= clock;

    % outputs
    Phi = [];
    
    % defaults
    epsilon = 0.001;                % Convergence criterion: Min change in beta * epsilon
    num_iter = 5;
    doverbose = 1;
    arorder = 0;                    % or Zero for no AR
    type = 'i';
    within = 'common';
    
    
    for varg = 1:length(varargin)
        if ischar(varargin{varg})
            switch varargin{varg}
                % reserved keywords
                case 'verbose', doverbose = 1;
                case 'noverbose', doverbose = 0;
                case {'iterations', 'iter'}, num_iter = varargin{varg+1};
                case 'type', type = varargin{varg+1}; varargin{varg+1} = [];
                case {'epsilon', 'eps'}, epsilon = varargin{varg+1};
                case {'ar', 'arorder'} , arorder = varargin{varg+1};
                case {'within_var', 'within'} , arorder = varargin{varg+1};
                otherwise, disp(['Unknown input string option: ' varargin{varg}]);
            
            end
        end
    end
    
    % check type
    if ~(strcmp(type, 'i') || strcmp(type,'r'))
        error('Type must be ''i'' (igls) or ''r'' (rigls)');
    end
        
    [T, sub] = size(y);                              % Length of y vector (Time) x Number of subjects
    
    Tc = 2*T;                                        % length of concatenated y, m vector
    Tc2 = Tc * (Tc + 1) ./ 2;                        % num. elements in lower triangle of cov matrix
    len = sub * Tc;                                  % Total number of observations
    n_G = sub * Tc2;                                 % number of rows in var-comp est. design matrix
    
    w = [y; m];
    z = reshape(w,len,1);                            % Concatenated data
    D = zeros(len,5);                                % Design matrix
    
    one = (zeros(T,1)+1);
    null = (zeros(T,1));
    
    for s=1:sub,
         wh = ((s-1) * Tc + 1):(s * Tc);  
         D(wh,:) = [[one null x(:,s) null m(:,s)]; [null one null x(:,s) null]];      % Create combined design matrix
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Step 1: Find the OLS solution
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    beta = D \ z;                                           % Beta values
    resid = z - D * beta;                                   % Residuals

    if (strcmp(within, 'common'))                           
        V = zeros(n_G,2);                                      
    else
        V = zeros(n_G,2*sub); 
    end
    
    Sigma_y = zeros(sub,1);                                       
    Sigma_m = zeros(sub,1);                                      
    
       
    if arorder > 0
        Phi = zeros(sub,arorder);
        get_ar                                             
    else
        get_V_no_ar
    end
    
    ystar = zeros(n_G, 1);                                      % Sums of squared residuals, concatenated across Ss
    get_ystar;



    % Fit the variance parameter design matrix to ystar, est. residual variances

    G = Create_Design(x, m, V);                                % Create design matrix for variance estimation
    betastar = G \ ystar;                                      % Estimate variance components
    
    ind = zeros(size(G,2),1);
    ind(1:5) = 1;
    betastar((ind > 0) & (betastar < 0)) = 0;                  % Use max(0,betastar) to ensure nonnegative variance.
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Step 2: Iterate
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    iterations = 0;

    min_change = betastar * epsilon;

    iSig = zeros(len,len);
    Sig = zeros(len,len);
    isconverged = 0;
    
    
    while (iterations < num_iter) && ~isconverged

        for s = 1:sub
            Sig1 = ivech(G(((s-1) * Tc2 + 1):(s * Tc2),:) * betastar);
            Sig1 = Sig1 + tril(Sig1,-1)';                      
            wh = ((s-1)*Tc+1):(s*Tc);                   % which indices in Cov mtx for this subject
            Sig(wh, wh) = (Sig1);
            iSig(wh, wh) = inv(Sig(wh, wh));
        end

        tmp = D'*iSig;
        tmp2 = inv(tmp*D);
        beta = tmp2*tmp*z;                              % Beta values

        resid = z - D*beta;                             % Residuals

        betastar_old = betastar;

        %         if arorder > 0
        %             get_ar                                              % updates Phi, Sigma, and -> V (from Sigma)
        %         else
        %             get_V_no_ar
        %         end
        
        get_ystar;
   
        betastar = G \ ystar;                                           % Estimate variance components
        betastar((ind > 0) & (betastar < 0)) = 0;                       % Use max(0,betastar) to ensure nonnegative variance.

        isconverged = ~any(abs(betastar - betastar_old) > abs(min_change));
        iterations = iterations + 1;

    end   

    %     Cov_beta = inv(D'*iSigma*D);
    %     W = (Sigma - D*inv(D'*iSigma*D)*D');        % Residual inducing matrix x Sigma
    %     df_beta = (trace(W).^2)./trace(W*W);        % Satterthwaite approximation for degrees of freedom

    Cov_beta = tmp2;
    df_beta = sub - 1; %Edit later
    
%    Cov_betastar = inv(G'*G);
    
    if (strcmp(within, 'common'))
        Sigma_y = betastar(9);
        Sigma_m = betastar(10);
    else
        Sigma_y = betastar(9:(sub+9-1));
        Sigma_m = betastar((sub+9):end);
    end;
            
    betastar = betastar(1:8);                       % Remove scaling factor corresponding to within subject variance
    
    % Cov_betastar = Cov_betastar(1:8,1:8);           % Remove scaling factor corresponding to within subject variance
      
    Cb = zeros(8,8);
    GG = zeros(Tc,Tc,8);
   
    NN = zeros(T,T);
    ONE = zeros(T,T)+1;
    GG(:,:,1) = [[ONE NN]; [NN NN]];    
    GG(:,:,2) = [[NN NN]; [NN ONE]]; 
        
    for i =1:sub,
        
        XX = x(:,i)*x(:,i)';
        MX = m(:,i)*x(:,i)';
        XM = x(:,i)*m(:,i)';
        MM = m(:,i)*m(:,i)';
        wh = ((s-1)*Tc+1):(s*Tc);        
        iS = iSig(wh,wh);
        
        GG(:,:,3) = [[XX NN]; [NN NN]];
        GG(:,:,4) = [[NN NN]; [NN XX]];
        GG(:,:,5) = [[MM NN]; [NN NN]];
        GG(:,:,6) = [[NN XX]; [XX NN]];
        GG(:,:,7) = [[MX+XM NN]; [NN NN]];
        GG(:,:,8) = [[NN MX]; [XM NN]];
        
        for h=1:8,
            for l=1:8,
                Cb(h,l) = Cb(h,l) + trace(GG(:,:,h)'*iS*GG(:,:,l)*iS');
            end;
        end;
        
    end;
    
    Cov_betastar = Cb;

    betastar = betastar(1:8);                       % Remove scaling factor corresponding to within subject variance
    Cov_betastar = Cov_betastar(1:8,1:8);           % Remove scaling factor corresponding to within subject variance
        
    df_betastar = sub - length(betastar);
    
    
    c2 = clock;
    elapsed_time = etime(c2, c1);
    
        
    % save output structure

    out = struct('beta', beta, 'betastar', betastar, 'Cov_beta', Cov_beta, 'Cov_betastar', Cov_betastar, ...
        'Sigma_m', Sigma_m, 'Sigma_y', Sigma_y, 'Phi', Phi, ...
        'type', type, 'arorder', arorder, 'isconverged', isconverged, ...
        'num_obs', T, 'sub', sub, 'num_iter', num_iter, 'epsilon', epsilon, ...
        'iterations', iterations, 'elapsed_time', elapsed_time);

    % save stats
    out.t = out.beta ./ sqrt(diag(out.Cov_beta));
    out.df_beta = df_beta;
    out.p = 2 * (1 - tcdf(abs(out.t), out.df_beta));

    out.t_randvariance = out.betastar ./ sqrt(diag(out.Cov_betastar));
    out.df_betastar = df_betastar;
    out.p_randvariance = 2 * (1 - tcdf(abs(out.t_randvariance), out.df_betastar)); 
    
    % print output
    
    if doverbose
        fprintf('igls.m finished:\n');
        fprintf('Data: %3.0f observations x %3.0f subjects \n', T, sub);
        typestr = {'igls' 'rigls'};
        convergestr = {'No' 'Yes'};
        
        fprintf('Fit Type: %s\n', typestr{strcmp(type, 'r') + 1});
        fprintf('Converged: %s\n', convergestr{isconverged + 1});
        fprintf('Max iterations: %3.0f, Completed iterations: %3.0f\n', num_iter, iterations);
        fprintf('Epsilon for convergence: %3.6f\n', epsilon);
        
        fprintf('Elapsed time: %3.2f s\n', elapsed_time);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Nested functions
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function get_V_no_ar
    
        NN = zeros(T,T);
        mysig1 = vech([[eye(T) NN];[NN NN]]);
        mysig2 = vech([[NN NN];[NN eye(T)]]);

        for k = 1:sub
            if (strcmp(within, 'common'))
                wh = ((k-1) * Tc2 + 1):(k * Tc2);                             % indices in time series for this subject    
                V(wh,1) = mysig1;                                             % Create regressors for common variance
                V(wh,2) = mysig2;
            else
                wh = ((k-1) * Tc2 + 1):(k * Tc2);                             % indices in time series for this subject
                V(wh,k) = mysig1;                                             % Create 2*sub regressors otherwise               
                V(wh,(k+sub)) = mysig2;
            end
        end
    end


    function get_ar

        for k = 1:sub

            % Estimate AR parameters using the Yuke-Walker method for each subject.     
            wh = ((k-1) * Tc + 1):((k-1) * Tc + T);
            Dt = D(wh,[1 3 5]);
            res = z(wh) - (Dt * (Dt \ z(wh)));
            [a,e] = aryule(res, arorder);             % Yule-Walker
            Phi = a(2:(arorder+1));                   % Parameters of AR(p) model
            Sigma_y(k) = sqrt(e);                     % standard deviation of AR(p) model

            % Find the covariance matrix
            A = diag(ones(T,1));
            for j=1:arorder,
              A = A + diag(Phi(j)*ones((T-j),1),-j);
            end

            iA = inv(A);
        %            C1 = Sigma_y(k)^2*(iA * iA');
            C1 = (iA * iA');
            
            wh = ((k-1) * Tc + T + 1):(k * Tc);
            Dt = D(wh,[2 4]);
            res = z(wh) - (Dt * (Dt \ z(wh)));
            [a,e] = aryule(res, arorder);             % Yule-Walker
            Phi = a(2:(arorder+1));                   % Parameters of AR(p) model
            Sigma_m(k) = sqrt(e);                     % standard deviation of AR(p) model

            % Find the covariance matrix
            A = diag(ones(T,1));
            for j=1:arorder,
              A = A + diag(Phi(j)*ones((T-j),1),-j);
            end

            wh = ((k-1) * Tc2 + 1):(k * Tc2);                             % indices in time series for this subject
            iA = inv(A);
            
    %            C2 = Sigma_m(k)^2*(iA * iA');
            C2 = (iA * iA');
            
            NN = zeros(T,T);            
            mysig1 = vech([[C1 NN];[NN NN]]);
            mysig2 = vech([[NN NN];[NN C2]]);

            
            if strcmp(within, 'common')
%                V(wh) = Sigma(k)^2 * tmp;                                      % Covariance function in vech format
                V(wh,1) = vech(mysig1);                                         
                V(wh,2) = vech(mysig2);
            else
%                V(wh,k) = Sigma(k)^2 * tmp;                                    % Covariance function in vech format
                V(wh,k) = vech(mysig1);              
                V(wh,sub+k) = vech(mysig2);                      
            end

            
        end
        
    end


   function get_ystar

        if (type == 'i')           % IGLS
            for k=1:sub
                wh = ((k-1) * Tc + 1):(k * Tc);                              % indices in time series for this subject
    %                myresid = resid( wh ) - V(wh);                          % residuals for this subject with within-subject variance removed.
                myresid = resid( wh );                                       % residuals for this subject.
                tmp = vech(myresid * myresid');                              % Find vech of estimated covariance

                wh = ((k-1) * Tc2 + 1):(k * Tc2);                            % indices in time series for this subject
                ystar(wh) = tmp;
            end
        elseif (type == 'r')       % RIGLS
            for k=1:sub
                wh = ((k-1) * Tc + 1):(k*Tc);                                % indices in time series for this subject
                Dtmp = D(wh, :);
    %                rtmp = resid(wh) - V(wh);                               % residuals for this subject with within-subject variance removed.
                rtmp = resid(wh);                                            % residuals for this subject.
                rig = rtmp * rtmp' + Dtmp * inv(Dtmp' * Dtmp) * Dtmp';
                tmp = vech(rig);                                             % Find vech of estimated covariance

                wh = ((k-1) * Tc2 + 1):(k * Tc2);                            % indices in time series for this subject
                ystar(wh) = tmp;
            end
        end

    end


end % END MAIN FUNCTION


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subfunctions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = Create_Design(x, m, V)
    % function G = Create_Design(x)
    %
    % Create Design matrix for estimation of variance components
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [T, sub] = size(x);
    
    NN = zeros(T,T);
    ONE = zeros(T,T)+1;

    Tc= 2*T;
    Tc2 = Tc * (Tc+1) / 2;    
    G = zeros(sub * Tc2, 8);
  
    for i = 1:sub

        XX = x(:,i)*x(:,i)';
        MX = m(:,i)*x(:,i)';
        XM = x(:,i)*m(:,i)';
        MM = m(:,i)*m(:,i)';
        
        wh = ((i-1) * Tc2 + 1):(i * Tc2); 
        
        G(wh,1) = vech([[ONE NN]; [NN NN]]);
        G(wh,2) = vech([[NN NN]; [NN ONE]]); 
        G(wh,3) = vech([[XX NN]; [NN NN]]);
        G(wh,4) = vech([[NN NN]; [NN XX]]);
        G(wh,5) = vech([[MM NN]; [NN NN]]);
        
        G(wh,6) = vech([[NN XX]; [XX NN]]);
        G(wh,7) = vech([[MX+XM NN]; [NN NN]]);
        G(wh,8) = vech([[NN MX]; [XM NN]]);
        
    end

    G = [G V];                              % use within-subject covariance as a regressor
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V = vech(Mat)
    % function V = vech(Mat)
    %
    % Calculate vech for the matrix Mat
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    V = Mat(tril(true(size(Mat))));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mat = ivech(V)
    % function Mat = vech(V)
    %
    % Calculate the "inverse" of the vech function
    % This is much faster than matlab's squareform.m
    % It could be speeded up, probably, by operating column-wise
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    len = length(V);
    dim = -0.5 + sqrt(0.25 + 2 * len);
    Mat = zeros(dim, dim);
    ind=1;

    for i=1:dim
        for j=i:dim
            Mat(j,i) = V(ind);
            ind = ind+1;
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
