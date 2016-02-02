function spm_orthviews_name_axis(name, axisnum)
% :Usage:
% ::
%
%    spm_orthviews_name_axis(name, axis#)
%
% put names on spm_orthviews axes

    global st

    axes(st.vols{axisnum}.ax{3}.ax)

    name(name == '_') = ' ';

    title(name,'FontSize',24,'Color','k')



end
