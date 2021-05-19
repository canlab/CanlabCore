function [S,KL,KH] = use_spm_filter(TR,dims,LChoice,HChoice,HParam,varargin)
% :Usage:
% ::
%
%     function [S,KL,KH] = use_spm_filter(TR,dim of filter,LChoice,HChoice,HP filter in s,[LP Gauss len in s])
%
% :Inputs:
%
%   **K{s}.LChoice:**
%        Low-pass  filtering {'hrf' 'Gaussian' 'none'}
%   **K{s}.LParam:**
%        Gaussian parameter in seconds
%   **K{s}.HChoice:**
%        High-pass filtering {'specify' 'none'}
% ..
%    05/22/01 Tor Wager
% ..

K{1}.RT = TR; 
K{1}.LChoice = LChoice;
K{1}.HChoice = HChoice;
K{1}.HParam = HParam; 
K{1}.row = ones(dims, 1);

if length(varargin) > 0
        K{1}.LParam = varargin{1};
end

    
KL = []; KH = [];

spmS = spm_filter('set', K);

S = eye(length(K{1}.row));

if ~strcmp(HChoice,'none')
    
	KH = full(spmS{1}.KH);
	S = S - KH * pinv(KH);
    
end

if ~strcmp(LChoice,'none')
	KL = full(spmS{1}.KL);
	S = KL * S;         	% lowpass * highpass; hp = I - res forming mtx of S.KH
end

return





function [vargout] = spm_filter(Action,K,Y)
% filter routine
% FORMAT [K] = spm_filter('set',K)
% FORMAT [Y] = spm_filter('apply',K,Y)
%
% Action    - 'set'   fills in filter structure K
% Action    - 'apply' applies K to Y = K*Y
% K         - filter convolution matrix or:
% K{s}      - cell of structs containing session-specific specifications
%
% K{s}.RT       - repeat time in seconds
% K{s}.row      - row of Y constituting session s
% K{s}.LChoice  - Low-pass  filtering {'hrf' 'Gaussian' 'none'}
% K{s}.LParam   - Gaussian parameter in seconds
% K{s}.HChoice  - High-pass filtering {'specify' 'none'}
% K{s}.HParam   - cut-off period in seconds
%
% K{s}.HP       - low frequencies to be removed
% K{s}.LP       - sparse toepltz low-pass convolution matrix
% 
% Y         - data matrix
%
% K         - filter structure
% Y         - filtered data K.K*Y
%___________________________________________________________________________
%
% spm_filter implements band pass filtering in an efficient way by
% using explicitly the projector matrix form of the High pass
% component.  spm_filter also configures the filter structure in
% accord with the specification fields if required
%___________________________________________________________________________
% @(#)spm_filter.m	2.4 Karl Friston 99/08/31


% set or apply
%---------------------------------------------------------------------------
switch Action

	case 'set'
	%-------------------------------------------------------------------
	for s = 1:length(K)

		% matrix order
		%-----------------------------------------------------------
		k     = length(K{s}.row);

		% make low pass filter
		%-----------------------------------------------------------
		switch K{s}.LChoice

			case 'none'
			%---------------------------------------------------
			h       = 1;
			d       = 0;

			case 'hrf'
			%---------------------------------------------------
			h       = spm_hrf(K{s}.RT);
			h       = [h; zeros(size(h))];
			g       = abs(fft(h));
			h       = real(ifft(g));
			h       = fftshift(h)';
			n       = length(h);
			d       = [1:n] - n/2 - 1;

			case 'Gaussian'
			%---------------------------------------------------
			sigma   = K{s}.LParam/K{s}.RT;
			h       = round(4*sigma);
			h       = exp(-[-h:h].^2/(2*sigma^2));
			n       = length(h);
			d       = [1:n] - (n + 1)/2;
			if      n == 1, h = 1; end

			otherwise
			%---------------------------------------------------
			error('Low pass Filter option unknown');
			return

		end

		% create and normalize low pass filter
		%-----------------------------------------------------------
		K{s}.KL = spdiags(ones(k,1)*h,d,k,k);
		K{s}.KL = spdiags(1./sum(K{s}.KL')',0,k,k)*K{s}.KL;


		% make high pass filter
		%-----------------------------------------------------------
		switch K{s}.HChoice

			case 'none'
			%---------------------------------------------------
			K{s}.KH = [];

			case 'specify'
			%---------------------------------------------------
			n       = fix(2*(k*K{s}.RT)/K{s}.HParam + 1);
			X       = spm_dctmtx(k,n);
			K{s}.KH = sparse(X(:,[2:n]));

			otherwise
			%---------------------------------------------------
			error('High pass Filter option unknown');
			return

		end

	end

	% return structure
	%-------------------------------------------------------------------
	vargout = K;


	case 'apply'
	%-------------------------------------------------------------------
	if iscell(K)


		% ensure requisite feild are present
		%-----------------------------------------------------------
		if ~isfield(K{1},'KL')
			K = spm_filter('set',K);
		end

		for s = 1:length(K)

			% select data
			%---------------------------------------------------
			y = Y(K{s}.row,:);

			% apply low pass filter
			%---------------------------------------------------
			y = K{s}.KL*y;

			% apply high pass filter
			%---------------------------------------------------
			if ~isempty(K{s}.KH)
				y = y - K{s}.KH*(K{s}.KH'*y);
			end

			% reset filtered data in Y
			%---------------------------------------------------
			Y(K{s}.row,:) = y;

		end

	% K is simply a convolution matrix
	%-------------------------------------------------------------------
	else
		Y = K*Y;
	end

	% return filtered data
	%-------------------------------------------------------------------
	vargout   = Y;


	otherwise
	%-------------------------------------------------------------------
	warning('Filter option unknown');


end



function C = spm_dctmtx(N,K,n,f)
% Creates basis functions for Discrete Cosine Transform.
% FORMAT C = spm_dctmtx(N,K,n)
%     OR C = spm_dctmtx(N,K)
%     OR D = spm_dctmtx(N,K,n,'diff')
%     OR D = spm_dctmtx(N,K,'diff')
% N - dimension
% K - order
% n - optional points to sample
%____________________________________________________________________________
% spm_dctmtx creates a matrix for the first few basis functions of a one
% dimensional discrete cosine transform.
% With the 'diff' argument, spm_dctmtx produces the derivatives of the
% DCT.
%
% See:    Fundamentals of Digital Image Processing (p 150-154).
%         Anil K. Jain 1989.
%____________________________________________________________________________
% @(#)spm_dctmtx.m	1.3 John Ashburner MRCCU/FIL 96/08/14

d = 0;

if (nargin == 2)
	n = (0:(N-1))';
	if (nargin == 3)
		d = 1;
	end
elseif (nargin == 3)
	if (strcmp(n,'diff'))
		d = 1;
		n = (0:(N-1))';
	else
		n = n(:);
	end
elseif (nargin == 4)
	n = n(:);
	if (strcmp(f,'diff'))
		d = 1;
	else
		error('Incorrect Usage');
	end
else
	error('Incorrect Usage');
end

C = zeros(size(n,1),K);

if (d == 0)
	C(:,1)=ones(size(n,1),1)/sqrt(N);
	for k=2:K
		C(:,k) = sqrt(2/N)*cos(pi*(2*n+1)*(k-1)/(2*N));
	end
else
	for k=2:K
		C(:,k) = -2^(1/2)*(1/N)^(1/2)*sin(1/2*pi*(2*n*k-2*n+k-1)/N)*pi*(k-1)/N;
	end
end


