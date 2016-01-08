function H = make_conv_mtx(sz,sampres)
% Constructs the matrix (H) for a linear convolution
% With the canonical SPM hrf
% such that Hx = conv(x,hrf)
%
% :Usage:
% ::
%
%     function H = make_conv_mtx(sz,sampres)
%
% :Inputs:
%
%   **sz:**
%        size of output matrix (elements)
%
%   **sampres:**
%        spm_hrf sampling resolution (~ TR), OR
%
%        if a vector, a custom HRF
%        sampled at the appropriate frequency.
%
% ..
% Tor Wager
% ..

if length(sampres) == 1
	ho = spm_hrf(sampres)';
else
	ho = sampres;
	if length(ho) == size(ho,1), ho = ho';, end
end

for i = 1:sz

	 h = [zeros(1,i-1) ho]; hz = sz-length(h);

	if hz > 0, 
		h = [h zeros(1,hz)];, 
	else, 
		h = h(1:sz);, 
	end

	H(i,:) = h;
end
