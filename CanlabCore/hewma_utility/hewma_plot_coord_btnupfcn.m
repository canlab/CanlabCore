function [dat,stats,mycov] = hewma_plot_coord_btnupfcn()
% function for loading and plotting hewma data from a voxel.  must have
% declared globals f,f2,VOL,EXPT
%
% :Usage:
% ::
%
%     [dat,stats,mycov] = hewma_plot_coord_btnupfcn
%
% see hewma_timeseries_plot, the shell function.

global f
global f2
global VOL
global EXPT

% get coordinate
coord = spm_orthviews('Pos');
coord = mm2voxel(coord',VOL);
mm = round(spm_orthviews('Pos')');



% extract data and plot
[dat,stats] = hewma_extract_voxel(EXPT,coord);

mycov = [];
% plot groups, if info there
if isfield(EXPT,'cov')
    
    f = findobj('Tag', 'hewma_grp_witherrbars');
    if isempty(f)
        f = tor_fig;
        set(f, 'Tag', 'hewma_grp_witherrbars');
    else
        figure(f)
        cla
    end

    set(gca,'FontSize',16)
    mycov = EXPT.cov(:,1);
    h = timeseries_btwngroups_plot(dat,mycov,60,1);
    
     hold on; 
     low = stats.cp_ind(mycov<0);
     hi = stats.cp_ind(mycov>0);
     plot([low; low],zeros(2,length(low)),'ys','MarkerFaceColor','b');
     plot([hi;hi],zeros(2,length(hi)),'ys','MarkerFaceColor','r');
end

return
