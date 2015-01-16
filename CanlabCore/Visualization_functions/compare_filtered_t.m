function compare_filtered_t(anatP,varargin)
% function compare_filtered_t(anatP,P1,P2, etc...)
% tor wager
%
% example:
% compare_filtered_t([],'rob_tmap_filtered_0001.img','rob_tmap_filtered_0002.img')
%
% Threshold spm T images and display them together in SPM orthviews
% threshold_spm_t(.005,22,0,'pos')
% compare_filtered_t([],'rfx0009/spmT_filtered_0002.img','rfx0011/spmT_filtered_0002.img', ...
% 'rfx0013/spmT_filtered_0002.img','rfx0015/spmT_filtered_0002.img','rfx0017/spmT_filtered_0002.img')

if isempty(anatP), 
    
    anatP = which('scalped_single_subj_T1.img');, 
    if isempty(anatP), anatP = spm_get(1,'img','Choose anatomical overlay image');,end
    
end

anatP = repmat(anatP,length(varargin),1);, 

spm_check_registration(anatP)

for i = 1:length(varargin)

    V = spm_vol(varargin{i});
    v = spm_read_vols(V);

    vv{i} = v;
    
    wh = find(abs(v) > 0);
    [x,y,z] = ind2sub(size(v),wh);
    Z = v(wh);
    XYZ = [x y z]';
    spm_orthviews('AddBlobs',i,XYZ,Z,V.mat)
    sxyz{i} = XYZ;
    sz{i} = Z';
end

go = 1;
if go
    
if length(varargin) > 1
    vv{1}(isnan(vv{1})) = 0; vv{2}(isnan(vv{2})) = 0;
    vvmax = squeeze(sum(sum(abs(vv{1} - vv{2}))));
    vvmax = find(vvmax == max(vvmax));
    mm = voxel2mm([0 0 vvmax]',V.mat); mm = mm(3);
    
        disp(['Maximally different slice btwn imgs 1 and 2 is slice ' num2str(vvmax) ' mm = ' num2str(mm)])
    spm_orthviews('Reposition',[0 0 mm])
    spm_orthviews('Xhairs','off')
    
    for i = 1:length(varargin)
        
        
    end
    V.M = V.mat;
    
    ZZ = cat(2,sz{:});  [tmp,x] = hist(ZZ,50);
	figure('Color','w'), hold on
	mycol = {'r' 'b' 'g' 'r--' 'b--' 'g--'};

    for i= 1:length(sz)
	zhist{i} = hist(sz{i},x);
	plot(x,zhist{i},mycol{i},'LineWidth',2),hold on;
	myleg{i} = varargin{i};
    end
    
    xlabel('t-score','FontSize',14),ylabel('Frequency','FontSize',14),legend(myleg)
    
    compare_slice(anatP,sxyz,sz,V,5)
    

 
 end
end
 
return
