function compare_slice(ovlP,sxyz,sz,VOL,varargin)
% :Usage:
% ::
%
%    compare_slice(ovlP,sxyz,sz,VOL)
%
% :Inputs:
%
%   **ovlP:**
%        is name of overlay image file (anatomical)
%
%   **sxzy:**
%        is a cell array, where each cell is a list of
%   VOXEL coordinates in 3 x n matrix
%
%   **sz:**
%        is a cell array, where each cell contains the z values
%        corresponding to sxyz.  This is used to create 
%        pseudocolor on the plots.
%
%   **VOL:**
%        is a structure containing the field M,
%        where M is the SPM99-style mapping matrix from voxel to mm 
%        space. This is used to map sxyz values onto the overlay image.
%
%        VOL must also contain field dim, which has the voxel dims
%        of the results image for all sxyz.
%        All sxyz must have the same dimensions!
%
% ..
%    tor wager, 2/16/03
% ..

if length(varargin) > 0, ptsize = varargin{1};, else, ptsize = 5;, end

docolbar = 1;

% -------------------------------------------------------------
% * Compare subjects on a slice
% -------------------------------------------------------------
bgV = spm_vol(ovlP); bgV = bgV(1);
v = spm_read_vols(bgV);
bgV.M = bgV.mat;


ss = 1;
while ~isempty(ss)
    
	ss = input(['Pick a slice (1:' num2str(VOL.dim(3)) ')']);
	
    if isempty(ss),break,end
    
	figure('Color','w'); 
    rc = 1; rc2 = 1; wh = 1;
    while rc * rc2 < length(sxyz),  %add +1 to include colorbar on this fig
        if wh == 1, rc2 = rc2+1;, else, rc = rc+1;, end
        wh = ~wh;
    end

	ssmm = voxel2mm([0 0 ss]',VOL.M);
	sso = mm2voxel(ssmm,bgV); sso = sso(3);

    % -------------------------------------------------------------
	% define color maps - biscale hot/cool
    % -------------------------------------------------------------

	% color map - hot
	% --------------------------------------------
	h1 = (0:1/99:1)';
	h2 = ones(size(h1)); 
	h3 = zeros(size(h1));
	h = [h1 h3 h3; h2 h1 h3; h2 h2 h1];
	h(1:50,:) = []; % take only red values to start
	% in new matlab: h = colormap(hot(300));

	% color map - winter
	% --------------------------------------------
	h1 = (0:1/249:1)';
	h2 = (1:-1/(249*2):.5)';
	h3 = zeros(size(h1));
	hc = [h3 h1 h2];

    % -------------------------------------------------------------
	% determine overall z-score range
    % -------------------------------------------------------------
    
	zrange = cat(2,sz{:}); 
	tmp = zrange(zrange > 0);
	tmpc = zrange(zrange < 0);

	if ~isempty(tmp)
		zrange = [min(tmp) max(tmp)];
		zh = zrange(1):(zrange(2)-zrange(1))./249:zrange(2);
		zh = round(zh*100);
	end

	if ~isempty(tmpc)
		zrangec = [min(tmpc) max(tmpc)];
		zhc = zrangec(1):(zrangec(2)-zrangec(1))./249:zrangec(2);
		zhc = round(zhc*100);
	end

    % -------------------------------------------------------------
	% loop through sets of input coordinates
    % -------------------------------------------------------------
    
	for i = 1:length(sxyz)
        
        % -------------------------------------------------------------
	    % create image of the slice
        % -------------------------------------------------------------
		subplot(rc,rc2,i), imagesc(v(:,:,sso)'); hold on; colormap gray
		set(gca,'YDir','normal'),axis image, axis off

        % -------------------------------------------------------------
	    % select xyz coordinates in slice and transform
        % -------------------------------------------------------------
		myxyz = sxyz{i}(:,sxyz{i}(3,:) == ss);
		myz = sz{i}(sxyz{i}(3,:) == ss);

		myxyzmm = voxel2mm(myxyz,VOL.M);
		myxyz = mm2voxel(myxyzmm,bgV,1)';
        
        
        % -------------------------------------------------------------
	    % find color for each xyz and plot
        % -------------------------------------------------------------        
		clear h2,clear wh
		for j = 1:length(myz)
			if myz(j) > 0, docool = 0; else, docool = 1;, end

			if docool,
				tmp = find(round(myz(j)*100) == zhc);
				if isempty(tmp), 
					tmp = find((zhc-round(myz(j)*100)).^2 == min((zhc-round(myz(j)*100)).^2));
				end
			else
				tmp = find(round(myz(j)*100) == zh);
				if isempty(tmp), 
					tmp = find((zh-round(myz(j)*100)).^2 == min((zh-round(myz(j)*100)).^2));
				end
			end

			wh(j) = tmp(1);

			if docool
				h2(j) = plot(myxyz(1,j),myxyz(2,j),'Color',hc(wh(j),:),'MarkerSize',ptsize,'MarkerFaceColor',hc(wh(j),:));
			else
				h2(j) = plot(myxyz(1,j),myxyz(2,j),'Color',h(wh(j),:),'MarkerSize',ptsize,'MarkerFaceColor',h(wh(j),:));
			end

		end
		if exist('h2') == 1, set(h2,'Marker','square'),end

	drawnow
	end

    % -------------------------------------------------------------
	% color scale bar - we must create by hand
    % -------------------------------------------------------------
    
    %cc = colormap(gray); cc(1:10,:) = repmat([1 1 1],10,1);    % [0 0 .3] for dark blue
    %colormap(cc)

    if docolbar 
	% only on the 1st time thru
    	% weird bug in this; make sep fig
	%subplot(rc,rc2,length(sxyz)+1),hold on

	zrange = cat(2,sz{:}); 
	tmp = zrange(zrange > 0);
	tmpc = zrange(zrange < 0);

    	if ~isempty(tmp)
    		figure('Color','w'); hold on;
		zh2 = zh./100;
    		%axis([0 .3 zh2(1) zh2(end)]),hold on
		%for i = 1:size(h,1), plot([0 1],[zh2(i) zh2(i)],'Color',h(i,:));, end
		%set(gca,'XTickLabel',''); % ylabel('Z-score')
		%h3 = get(gcf,'Position');
    		%set(gcf,'Position',[h3(1:2) h3(3)*.3 h3(4)*.5])
            
            for i = 2:size(h,1), fill([zh2(i-1) zh2(i-1) zh2(i) zh2(i)],[0 1 1 0],h(i,:),'EdgeColor','none');, end
		    set(gca,'YTickLabel',''); % ylabel('Z-score')
            xlabel('Z-score','FontSize',14)
            
    		docolbar = 0;
    	end


    	if ~isempty(tmpc)
    		figure('Color','w'); hold on;
		zh2 = zhc./100;
    		axis([0 .3 zh2(1) zh2(end)]),hold on
		for i = 1:size(h,1), plot([0 1],[zh2(i) zh2(i)],'Color',hc(i,:));, end
		set(gca,'XTickLabel',''); % ylabel('Z-score')
		h3 = get(gcf,'Position');
    		set(gcf,'Position',[h3(1:2) h3(3)*.3 h3(4)*.5])
    		docolbar = 0;
    	end

    
    end	% make colorbar figures

end	% end loop
