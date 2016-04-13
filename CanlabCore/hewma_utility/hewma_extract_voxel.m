function [dat,stats] = hewma_extract_voxel(EXPT,coords)
% Extracts data for one voxel and runs hewma2 on it.
%
% :Usage:
% ::
%
%     function [dat,stats] = hewma_extract_voxel(EXPT,coords)
%
% Uses plot option in hewma2 to create plots.

dat = [];
stats = [];

persistent timeseriesfile
persistent lamfile
persistent files
persistent vfiles

zsl = num2str(coords(3));

mm_coords = spm_orthviews('Pos')';
fprintf(1, 'mm coordinates are (x, y, z) = (%3.0f, %3.0f, %3.0f)\t\t Voxel coords: (%3.0f, %3.0f, %3.0f)\n', mm_coords, coords);

if isempty(files) || isempty(vfiles)
    mydir = spm_get(-1,'*','Select Directory above individual subject directories',pwd);  
end

if isempty(files)
    
    % EXPT = getfunctnames2(EXPT,['z_slice' zsl '.mat'],'tmp');
    EXPT = getfunctnames2(EXPT,['z_slice' zsl '.mat'], 'tmp', [], mydir);
    
    files = str2mat(EXPT.tmp{:});
end

if isempty(vfiles)
    %EXPT = getfunctnames2(EXPT,['var_slice' zsl '.mat'],'tmp');
    EXPT = getfunctnames2(EXPT,['var_slice' zsl '.mat'], 'tmp', [], mydir);
    
    vfiles = str2mat(EXPT.tmp{:});
end

if isempty(files) || all(files(:) == ' ')
    disp('Cannot find z_slice files for subjects.  Check that path names are correct in EXPT structure.');
    return
end

if isempty(vfiles) || all(vfiles(:) == ' ')
    disp('Cannot find var_slice files for subjects.  Check that path names are correct in EXPT structure.');
    return
end

if isempty(timeseriesfile)
    
    timeseriesfile = spm_get(1,'*mat','Select hewma_timeseries.mat containing img dims',pwd);
  
end

  load(timeseriesfile, 'xdim', 'ydim');
% % try
% %     load hewma0001/hewma_timeseries xdim ydim
% % catch
% %     if ~(exist('hewma_timeseries') == 2)
% %             file = spm_get(1,'*mat','Select hewma_timeseries.mat',pwd);
% %             [dd,ff,ee] = fileparts(file);
% %             cd(dd)
% %             load hewma_timeseries xdim ydim
% %     end
% % end

if isempty(lamfile)
    lamfile = spm_get(1,'*mat','Select lam.mat file containing lambda parameter for one subject',pwd);  
end

load(lamfile, 'lam');
 
% % try
% %     load([EXPT.subjects{1} filesep 'lam'])    % load lambda param, assume same for all subjects
% % catch
% %     error('Cannot find lam.mat.  Start in dir above individual ewma results dirs.')
% % end


fprintf(1,'Loading data, subject ');

for i = 1:length(EXPT.subjects)
    fprintf(1,'%3.0f ',i);
    
    % load EWMA stat
    load(deblank(files(i,:)));
    zdat = full(zdat);
    T = size(zdat,2);
    zdat = reshape(zdat,xdim,ydim,T);


    try
        dat(i,:) = squeeze(zdat(coords(1),coords(2),:));
    catch
        disp('Warning!  data for this subject does not match others in size.  ')
        dat(i,:) = squeeze(zdat(coords(1),coords(2),1:size(dat,2)));
    end
    
    % load variance
    load(deblank(vfiles(i,:))); 
    vardat = full(vardat);
    vardat = reshape(vardat,xdim,ydim,T);

    try
        vdat(i,:) = squeeze(vardat(coords(1),coords(2),:));
    catch
        vdat(i,:) = squeeze(vardat(coords(1),coords(2),1:size(vdat,2)));
    end
    
end
fprintf(1,'\n');

if isfield(EXPT,'cov');
    mycov = EXPT.cov;
end
if ~isempty(mycov), mycov = mycov(:,1); end
        
[p,tm,Zcor,sbmean,Zpop,tvals,sb,stats] = hewma2(dat, vdat, lam,1,1,mycov);

%figure;plot(dat');
%hold on; plot(mean(dat),'k','LineWidth',2);
%title('EWMA statistics with mean');


return

