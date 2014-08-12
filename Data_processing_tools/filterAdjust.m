function [y,O,X,S] = filterAdjust(O)
% [y,O,X,S] = filterAdjust(OPTIONS)
%
%	O.y = signal
%  O.HP = high pass freq. cutoff
%  O.TR = sampling rate in s
%
%  O.doHP = [0 or 1] - do HP filter (spm), default is 0
%  O.doLP = [0, 1, 2] - do LP filter (spm), default is 0
%           2 = Gaussian filter with TR*2 s length
%
% O.firstimg - sets values for first image in each run to mean of the
% remaining values.  Good for removing first-image artifacts present in
% some scanner sequences.  Default is 1.
%
% O.cyclecorrection - checks for unimodal (normally distributed) data
% within each session, because some bimodal data that cycles between two
% mean scanner values has been observed in some data.  if a high proportion
% of outliers are found in non-normal data, subtracts mean of higher mode
% to adjust data. Sorry-not clear.  Check the code. Default is 0
%
% O.cyclecorrection2 - removes large transitions from data, as in
% detransition.m.  Default is 0.  Artifacts may be acquisition or
% motion-correction/resampling related.
%
% We do this after filtering, but if there's trouble, we re-do the
% scanadjust and filtering, because 'cycling' can affect these.
%
%  O.scanadjust = [0 or 1] - adjust to scan means, default is 0
%           - If O.X is entered, assumes this is the session mean matrix,
%           instead of recomputing.
%  O.percent = [0 or 1] - adjust to percent change from baseline, default is 0
%
%	O.filtertype = filter style, default is 'none'
%			'spm'		= use spm's filtering
%                       - if O.S is entered, uses this instead of
%                       recomputing
%			'fourier' = Doug's fourier filter
%           'fouriernotch' : Omit frequencies between HP(1) and HP(2)
%			'cheby'	= chebyshev
%			'Luis'		= Luis' custom filter
%			'none'		= no filtering (or leave field out).
%
% O.HP
%           for SPM, the filter cutoff in s
%           for fourier, the HP value or the [HP LP] values
%           -notches out everything slower than HP and faster than LP, in s
%           (1/s = Hz).  
%
%	O.nruns = number of runs (scanadjust), default is 1
%  O.adjustmatrix = custom adjustment matrix to regress out (e.g., movement params)
%	O.plot [0 or 1] 	= plots intermediate products, default is 0
%   O.verbose [0 or 1]   = verbose output
%   O.trimts [0 or std]  = trim overall timseries values to std, 0 for no trimming  
%   O.lindetrend         = specify linear detrending of timeseries.
%                       . occurs after adjustment and filtering and
%                       windsorizing
%                       . detrending option -> what to enter in this field:
%                       .   no detrending  -> empty, missing field, or 0
%                       .   detrending every n elements -> single number (n)
%                       .   piecewise linear detrend -> ROW vector of breakpoints
%                           (do not specify 1 as the start of the 1st segment.)
%
% by Tor Wager, 09/24/01
% modified 10/15/02 to add trial baseline adjustment and linear detrending.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% out = voistat('adjusty',y,HChoice,TR,nruns, [opt] adjustmatrix); 
% disp('*	voistat adjusty')

if ~isfield(O,'verbose'), O.verbose = 0;,end
if ~isfield(O,'trimts'), O.trimts = 0;,end
if ~isfield(O,'lindetrend'), O.lindetrend = 0;, end
if ~isfield(O,'percent'), O.percent = 0;, end

y = O.y;
TR = O.TR;
npoints = size(y,1);

if size(y,2) > size(y,1), y = y';, end

nruns = 1;
if isfield(O,'scanadjust'), if O.scanadjust, nruns = O.nruns;,end,else, O.scanadjust = 0;,end
if isfield(O,'adjustmatrix'), docustomadjust = 1;,adjustmatrix = O.adjustmatrix;,else, docustomadjust = 0;,adjustmatrix = [];,end

scanlen = npoints./nruns;
if scanlen ~= round(scanlen),warning(['Scan length is not an integer! scanlen = ' num2str(scanlen) ', rounded = ' num2str(round(scanlen))]),end
scanlen = round(scanlen);

wholeyavg = mean(y(~isnan(y)));
X = [];

if ~isfield(O,'doHP'), O.doHP = 0;,HChoice = [];,else, HChoice = O.HP;end
if O.doHP,if isempty(HChoice),error('filterAdjust: specify HP filter length in input options!'),end,end
if ~isfield(O,'doLP'), O.doLP = 0;,end
if ~isfield(O,'filtertype'), O.filtertype = 'none';,end
if ~isfield(O,'plot'), O.plot = 0;,end
if ~isfield(O,'cyclecorrection'), O.cyclecorrection = 0;,end
if ~isfield(O,'cyclecorrection2'), O.cyclecorrection2 = 0;,end
if ~isfield(O,'firstimg'), O.firstimg = 1;,end

% plot raw
if O.plot, f1 = figure('Color','w'); hold on; set(gca,'FontSize',16); plot(scale(y),'k');,title('Raw timeseries (standardized for display)'); pause(1);,end

if ((O.doHP | O.doLP) & strcmp(O.filtertype,'spm'))
	% MAKE FILTER - use_spm_filter.m
	% --------------------------------------------------------
	% Ran out of memory doing this with whole ts, so trying to make more efficient...do it scan by scan.
   if O.doHP | O.doLP
       if isfield(O,'S')
           if O.verbose, disp(	'		...voistat adjusty: Using input S filter: filter shape plot is final S matrix only'), end
           S = O.S; KL = O.S; KH = O.S;
       else
           
        
           if O.verbose, disp(	'		...voistat adjusty: making S filter'), end
        
        
           if O.doHP & O.doLP
            
               if O.doLP == 2
                
                   [S,KL,KH] = use_spm_filter(TR,scanlen,'Gaussian','specify',HChoice,TR*2);
            
               else
                
                   [S,KL,KH] = use_spm_filter(TR,scanlen,'hrf','specify',HChoice);
            
               end
        
        
           elseif O.doHP
         
               [S,KL,KH] = use_spm_filter(TR,scanlen,'none','specify',HChoice);
        
           elseif O.doLP == 2
         
               [S,KL,KH] = use_spm_filter(TR,scanlen,'Gaussian','none',[],TR*2);
        
           elseif O.doLP
         
               [S,KL,KH] = use_spm_filter(TR,scanlen,'hrf','none',[]);
        
           else error('Script bug.  Check filterAdjust.m')
        
           end
       
       end  % if S
      
   
   else warning('No filtering specified, though spm is filtertype.'), 
      S = [];
   end

	if O.plot, figure;
   		subplot(1,3,1);imagesc(S); title('Smoothing matrix (spm) *used')
   		subplot(1,3,2);imagesc(KL); title('Low-pass filter(spm)')
      subplot(1,3,3);imagesc(KH); title('High-pass filter(spm)')
      pause(1)
      close
    end
   
end

if (O.doHP | O.doLP) & strcmp(O.filtertype,'none')
    error('For filtering to be used, enter ''spm'' ''cheby'' ''doug'', etc. in O.filtertype.')
end

    
if O.cyclecorrection2
    if O.verbose, disp('		...removing large transitions from data (presumably artifacts).'),end
	
    for startimg = 1:scanlen:npoints
        tmp = y(startimg:startimg+scanlen-1);
        tmp = detransition(tmp,O.plot);
        y(startimg:startimg+scanlen-1) = tmp;
    end
    
    if O.plot, figure(f1); hold off; plot(scale(y),'k'); hold on; title('Raw (detransitioned, std for display)');,pause(1);end
end

   
if O.scanadjust
	% adjust scan means
	% --------------------------------------------------------
	if O.verbose, disp('		...adjusting scan means to session mean.'),end

    if isfield(O,'X')
        X = O.X;
    else
	    index = 1;
	    for startimg = 1:scanlen:npoints
		    X(startimg:startimg+scanlen-1,index) = 1;
		    index = index + 1;
	    end
    end
    

	if size(X,1) > size(y,1)
		disp(['	* * * voistat.m adjusty: WARNING: X is longer than y.  truncating.'])
		X = X(1:length(y),:);
	end

   %y = y - (X * (X \ y));	% y-xB to adjust to mean.
   y = y - X * pinv(X) * y;
   
   if O.plot, figure;
      imagesc(X); title('Model matrix for effects to regress out.')
      pause(1),close
      
      figure(f1); hold off; plot(scale(y),'k'); hold on; title('Session means removed, std for display)');,pause(1);
   end
   
end

if O.firstimg
    if O.verbose, disp('		...Setting values of first 2 images in each run to timeseries mean value.'),end
    tmp = 1:scanlen:length(y);
    tmp = [tmp 2:scanlen:length(y)];
    
    % first adjust to session means - rough estimate
    tmpn = 1:length(y); tmpn(tmp) = [];
    y(tmp) = mean(y(tmpn));
    
    % then adjust to expected 
    tmpx = (1:length(y))';
    [P] = polyfit(tmpx,y,3);
    yy = polyval(P,tmpx);
    y(tmp) = yy(tmp);
    
    if O.plot, figure(f1); plot(scale(y),'m'); hold on; 
        plot(scale(yy));
        title('Raw (1-2 removed, std for display)');,
    end
end

switch O.filtertype

case 'spm'
% filter each scan
% --------------------------------------------------------
if ~isempty(S)
   filty = [];
	if O.verbose, disp('		...applying to each session: high-pass filtering.'),end

	% whos y

	for startimg = 1:scanlen:npoints
		try
			sessy = y(startimg:startimg+scanlen-1,:);
		catch
			startimg
			whos y
			error(['	voistat adjusty: ERROR.  y is too short, npoints is wrong, or numscans is wrong? Check GUI values.'])
		end

		try
			sessy = S * sessy;
		catch
			try
				sessy = (S * sessy')';
			catch
				disp('		can''t multiply S and sessy: matrix dimensions not correct. Try transposing y?')
				whos S
				whos sessy
				disp(['Scanlength = ' num2str(scanlen)])
				disp(['total time pts = ' num2str(npoints)])
				disp(['Num Sessions = ' num2str(nruns)])
				error('Exiting now...')
			end
		end

		filty = [filty; sessy];
   end

% whos filty
   
   if O.plot
      figure(f1);hold off; plot(scale(y),'b','LineWidth',1);, hold on; plot(scale(filty),'r','LineWidth',1);
      legend({'Before filtering' 'After filtering'})

      %pause(1),close
   end
   
   y = filty;
   
end   
   
   % Old custom filter stuff.
       %Wn = .35;
       %[B,A] = cheby2(3,40,Wn);
       %y = filter(B,A,y); 

% Doug's custom fourier LP filter. 
%len = length(y) /2;
%ind = [1:len] ./ len; 
%ff = 1./(1+exp((ind-.30)./0.05));
%ff2 = [ff ff(len:-1:1)];
%y = real(ifft(fft(y).*ff2'));



case 'cheby'
if ~isempty(HChoice)
	nyquist = TR/2;	% in seconds, not frequencies TR/nyquist = .5
	Wn = [TR/HChoice .35];
	if O.verbose, disp([		'voistat adjusty: Chebyshev filter freqency is ' num2str(Wn)]),end
	for startimg = 1:scanlen:npoints
		sessy = y(startimg:startimg+scanlen-1);
 		[B,A] = cheby2(1,20,Wn);
		sessy = filter(B,A,y); 
	end
end

case 'fourier'
if ~isempty(HChoice)
	filty = []; nyquist = (1/TR)/2; lowerlim = 1/HChoice(1);
	hzfreqs = (1:scanlen) ./ (scanlen * TR);
	if length(HChoice) == 1
        ftopass = hzfreqs >= lowerlim & hzfreqs <= nyquist;
    else
        ftopass = hzfreqs >= lowerlim & hzfreqs <= 1./(HChoice(2));
    end
	if O.verbose, disp(['		...applying to each session: notch filtering: ' num2str(lowerlim) ' to ' num2str(nyquist) ' Hz.']),end
	for startimg = 1:scanlen:npoints
		myfft = fft(y(startimg:startimg+scanlen-1));
		mypass = zeros(1,scanlen);
		mypass(ftopass) = myfft(ftopass);
		% plot(abs(myfft));hold on;plot(abs(mypass),'rx');	
		sessy = real(ifft(mypass))';
		filty = [filty; sessy];
	end
	y = filty;
    %hold on; plot(y,'r');
end

case 'fouriernotch'
if ~isempty(HChoice)
	filty = []; nyquist = (1/TR)/2; lowerlim = 1/HChoice(1);
	hzfreqs = (1:scanlen) ./ (scanlen * TR);
    ftopass = ones(size(hzfreqs));
    ftopass(hzfreqs > lowerlim & hzfreqs < 1./(HChoice(2))) = 0;
    ftopass = find(ftopass);
    
	if O.verbose, disp(['		...applying to each session: notch filtering: ' num2str(lowerlim) ' to ' num2str(nyquist) ' Hz.']),end
	for startimg = 1:scanlen:npoints
		myfft = fft(y(startimg:startimg+scanlen-1));
		mypass = zeros(1,scanlen);
		mypass(ftopass) = myfft(ftopass);
		%figure; plot(abs(myfft));hold on;plot(abs(mypass),'rx');	
		sessy = real(ifft(mypass))';
		filty = [filty; sessy];
	end
	y = filty;
    %hold on; plot(y,'r');
end

case 'bandpass'
if ~isempty(HChoice)
	filty = []; nyquist = (1/TR)/2; lowerlim = 1/HChoice(1);
	hzfreqs = (1:scanlen) ./ (scanlen * TR);
    ftopass = zeros(size(hzfreqs));
    ftopass(hzfreqs > lowerlim & hzfreqs < 1./(HChoice(2))) = 1;
    ftopass = find(ftopass);
    
	if O.verbose, disp(['		...applying to each session: band pass filtering: ' num2str(lowerlim) ' to ' num2str(1./(HChoice(2))) ' Hz.']),end
	for startimg = 1:scanlen:npoints
		myfft = fft(y(startimg:startimg+scanlen-1));
		mypass = zeros(1,scanlen);
		mypass(ftopass) = myfft(ftopass);
		%figure; plot(abs(myfft));hold on;plot(abs(mypass),'rx');	
		sessy = real(ifft(mypass))';
		filty = [filty; sessy];
	end
	y = filty;
    %hold on; plot(y,'r');
end


case 'Luis'
   if O.plot
      y = luisFilter(y,O.TR,.35,'v');pause(1),close,pause(1),close
   else
      y = luisFilter(y,O.TR,.35);
   end
   
   
case 'none'
   
end % end switch



if O.cyclecorrection
    % check for bimodal 'cycling' of mean BOLD signal within sessions, and
    % correct if necessary, before removing scan means, etc.
    % if we find problems, then re-do session-centering
	% --------------------------------------------------------
    if O.verbose, disp('		...checking for normally distributed data (to avoid bimodal cycling observed in some data).'),end
	cycor = []; ind = 1; classes = [];
    for startimg = 1:scanlen:npoints
		   tmp = y(startimg:startimg+scanlen-1);

           [tmp,IDX] = mixture_model(tmp,O.plot);
                
           if max(IDX) > 1 & ~isempty(S)
               tmp = S * tmp;
               if O.plot, subplot(1,3,1),hold on;plot(tmp,'g'),end
           end
           
           y(startimg:startimg+scanlen-1) = tmp - mean(tmp);
           classes(ind) = max(IDX);
           
            ind = ind+1; 
    end
    cycor = classes > 1;
    
        if O.verbose, 
            
            fprintf(1,'		...%3.0f Sessions were not normally distributed.\n',sum(cycor)),
            if any(cycor), fprintf(1,['\t\t\tClasses by session: ' num2str(classes) '\n']),end
        else
            if any(cycor), fprintf(1,['\t\t\tClasses by session: ' num2str(classes) '\n']),end
        end
end

    
% do overall trimming here
% --------------------------------------------------------
if O.trimts
    [y,ntrimmed] = trimts(y,O.trimts,[],1);
    if O.verbose,fprintf(1,['\t\t...Windsorized ' num2str(ntrimmed) ' points from overall timeseries to ' num2str(O.trimts) ' std. deviations.\n']), end
    
    if O.plot
       figure(f1); plot(scale(y),'k--');, legend({'Before filtering' 'After filtering' 'After filt+trimming'})
    end
    %tpts = abs(y) > O.trimts .* std(y);
    %y(tpts) = NaN;
end



if O.lindetrend
     % linear detrending at specified knot points
     
    if O.verbose, fprintf(1,['\t\t...piecewise linear detrending requested with ' num2str(length(O.lindetrend)) ' breakpoints.\n']),end
    % set bp to [] for removal of one linear trend from the whole column
    if O.lindetrend == 1, O.lindetrend = [];, end    
    % if a single number, interpret as knotpoints every n elements.
    
    if length(O.lindetrend) == 1, O.lindetrend = O.lindetrend:O.lindetrend:length(y);, end  
    
    if O.plot, 
        figure; plot(y);, 
        set(gcf,'Position',[144         542        1306         504],'Color','w')
        disp(['         Breakpoints: ' num2str(O.lindetrend)])    
    end
    hold on;
    y = detrend(y,'linear',O.lindetrend);
    y = y + wholeyavg;  % preserve mean of the original timeseries
    
    if O.plot,
        figure(f1); hold on; plot(y,'m');
        plot([O.lindetrend' O.lindetrend']',[ones(length(O.lindetrend),1) * min(y) ones(length(O.lindetrend),1)*max(y)]','k')
        legend({'Original' 'Detrended'})
        
        mystd = std(y) * 3;
        plot([0 length(y)],[mean(y)+mystd mean(y)+mystd],'k--')
        plot([0 length(y)],[mean(y)-mystd mean(y)-mystd],'k--')
        text(length(y)-10,mean(y)+mystd,'3 std lines')
        
        drawnow
        pause(1)
        %close
    end
end
    
    
    

if O.percent
 
	if O.verbose, disp(['		...converting to % change from overall mean:  original mean is ' num2str(wholeyavg)]),end
	y = (y-mean(y))*100 / wholeyavg;		% percent change from 0.

end

if docustomadjust
		if O.verbose, disp('		...adjusting for user-entered covariates -e.g. movement params'),end
 		if sum(sum(isnan(adjustmatrix))) > 0,
			disp('		...WARNING: NaN values in user covariates!  Setting these to 0.')	
			adjustmatrix(isnan(adjustmatrix)) = 0;
		end
		try
			X = [adjustmatrix];
            X(:,end+1) = 1;     % add intercept
		catch
			warning('		voistat ''adjusty'': X and adjustmatrix are different sizes.  no covariates entered.')
			whos X
			whos adjustmatrix
        end
        O.custombeta = pinv(X) * y;
        y = y - X * O.custombeta;
 end
    


% figure; plot(abs(fft(y)),'ro');set(gca,'YLim',[0 1000])


return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
