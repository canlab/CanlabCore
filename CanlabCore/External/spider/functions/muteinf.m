function info = muteinf(A, Y, method)
	n = size(A,1);	
	Z = [A Y];	
	if(n/10 > 20)
		nbins = 20;
	else
		nbins = max(floor(n/10),10);
	end;
	pA = hist(A, nbins);
	pA = pA ./ n;
	
	i = find(pA == 0); % that's a hack!
	pA(i) = 0.00001;
	if(strcmp(method,'regression'))
		pY = hist(Y, nbins);
		pY = pY ./ n;
		j = find(pY == 0);
		pY(j) = 0.00001;	
					
		rx = abs(max(A) - min(A)) / nbins;
		ry = abs(max(Y) - min(Y)) / nbins;	
		yl = min(Y);
		p = zeros(nbins, nbins);	
		for i = 1:nbins
			xl = min(A);
			for j = 1:nbins
				%disp(['intervals [' num2str(xl) ',' num2str(xl+rx) '|' num2str(yl) ',' num2str(yl+ry) ']'])
				interval = (xl <= Z(:,1)) & (yl <= Z(:,2));
				if(j < nbins)
					interval = interval & (Z(:,1) < xl + rx);
				end;
				if(i < nbins)
					interval = interval & (Z(:,2) < yl + ry);
				end;			
				%find(interval)			
				p(i,j) = length(find(interval));

				if p(i,j) == 0 % hack!
					p(i,j) = 0.00001;
				end

				xl = xl + rx;
			end;
			yl = yl + ry;
		end;	
		HA = -sum(pA .* log(pA));
		HY = -sum(pY .* log(pY));											
		pA = repmat(pA,nbins,1);
		pY = repmat(pY',1,nbins);	
	else
		od = size(Y,2);
		cl = od;
		if(od == 1)
			pY = [length(find(Y==+1)) length(find(Y==-1))] / n;			
			cl = 2;
		else
			pY = zeros(1,od);
			for i=1:od
				pY(i) = length(find(Y==+1));
			end;
			pY = pY / n;			
		end;
		p = zeros(cl,nbins);
		rx = abs(max(A) - min(A)) / nbins;
		for i = 1:cl
			xl = min(A);
			for j = 1:nbins				
				if(i == 2) & (od == 1)
					interval = (xl <= Z(:,1)) & (Z(:,2) == -1);	
				else
					interval = (xl <= Z(:,1)) & (Z(:,i+1) == +1);
				end;
				if(j < nbins)
					interval = interval & (Z(:,1) < xl + rx);
				end;				
				%find(interval)			
				p(i,j) = length(find(interval));

				if p(i,j) == 0 % hack!
					p(i,j) = 0.00001;
				end

				xl = xl + rx;	
			end;
		end;
		HA = -sum(pA .* log(pA));
		HY = -sum(pY .* log(pY));
		pA = repmat(pA,cl,1);
		pY = repmat(pY',1,nbins);
	end;
	p = p ./ n;
	
	
	info = sum(sum(p .* log(p ./ (pA .* pY))));	
	info = 2 * info ./ (HA + HY);